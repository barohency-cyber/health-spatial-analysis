#!/usr/bin/env python3
"""空间计量入门：OLS、空间滞后模型（SLM）和空间误差模型（SEM）。

模型：
    OLS: y = Xβ + ε
    SLM: y = ρWy + Xβ + ε
    SEM: y = Xβ + u,  u = λWu + ε

估计：最大似然；ρ/λ 用网格 + 黄金分割搜索，β/σ² 在给定 ρ/λ 下用 OLS。
标准误：在最优参数处用数值 Hessian 近似，仅作为教学示例；正式论文建议用
spreg / GeoDa 复核，并报告权重敏感性。

输入：
    output/accessibility.csv   可达性结果
    data/outcomes.csv          模拟就医利用结果与协变量

输出：
    output/spatial_regression_coefficients.csv
    output/spatial_regression_fits.csv
    output/spatial_regression_fitted.csv
    output/spatial_regression_summary.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ACCESS = ROOT / "output" / "accessibility.csv"
DEFAULT_OUTCOMES = ROOT / "data" / "outcomes.csv"
DEFAULT_OUTPUT = ROOT / "output"
DEFAULT_PREFIX = "spatial_regression"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spatial_autocorr import (  # noqa: E402
    Observation,
    _moran_i,
    build_knn_weights,
)


def _normal_p(z_value: float) -> float:
    return math.erfc(abs(z_value) / math.sqrt(2))


def _chi2_1_p(statistic: float) -> float:
    if statistic <= 0:
        return 1.0
    return math.erfc(math.sqrt(statistic / 2))


def read_access(path: Path) -> list[Observation]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except FileNotFoundError as exc:
        raise ValueError(f"找不到文件：{path}；请先运行 python code/two_step_fca.py") from exc
    if not rows:
        raise ValueError(f"{path} 没有数据行")
    observations = []
    for line_no, row in enumerate(rows, start=2):
        try:
            observations.append(
                Observation(
                    community_id=row["community_id"],
                    name=row["name"],
                    lon=float(row["lon"]),
                    lat=float(row["lat"]),
                    population=float(row["population"]),
                    access_per_1000=float(row["access_per_1000"]),
                )
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{path} 第 {line_no} 行无法解析") from exc
    return observations


def read_outcomes(path: Path) -> dict[str, dict[str, float]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except FileNotFoundError as exc:
        raise ValueError(
            f"找不到文件：{path}；请先运行 python code/make_outcome_data.py"
        ) from exc
    required = {"community_id", "outpatient_utilization", "elderly_share", "income_index"}
    if rows and not required.issubset(rows[0]):
        missing = sorted(required - set(rows[0]))
        raise ValueError(f"{path} 缺少字段：{', '.join(missing)}")
    result = {}
    for row in rows:
        community_id = (row.get("community_id") or "").strip()
        if not community_id:
            continue
        result[community_id] = {
            "outpatient_utilization": float(row["outpatient_utilization"]),
            "elderly_share": float(row["elderly_share"]),
            "income_index": float(row["income_index"]),
        }
    return result


def build_data(
    observations: Sequence[Observation], outcomes: dict[str, dict[str, float]]
) -> tuple[np.ndarray, np.ndarray, list[str], list[Observation]]:
    y_values = []
    rows = []
    matched = []
    for observation in observations:
        outcome = outcomes.get(observation.community_id)
        if outcome is None:
            continue
        y_values.append(outcome["outpatient_utilization"])
        rows.append(
            [
                1.0,
                observation.access_per_1000,
                outcome["elderly_share"],
                outcome["income_index"],
            ]
        )
        matched.append(observation)
    if len(matched) < 10:
        raise ValueError(
            "匹配后的社区太少；请确认 accessibility.csv 与 outcomes.csv 的 community_id 一致"
        )
    return (
        np.asarray(y_values, dtype=float),
        np.asarray(rows, dtype=float),
        ["const", "access_per_1000", "elderly_share", "income_index"],
        matched,
    )


def _ols(y: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    residual = y - x @ beta
    rss = float(residual @ residual)
    sigma2 = rss / len(y)
    log_likelihood = -0.5 * len(y) * (math.log(2 * math.pi) + math.log(sigma2) + 1)
    return beta, residual, sigma2, log_likelihood


def _fit_stats(
    y: np.ndarray,
    x: np.ndarray,
    beta: np.ndarray,
    log_likelihood: float,
    sigma2: float,
    parameter_count: int,
    yhat: np.ndarray | None = None,
) -> dict[str, float]:
    n = len(y)
    if yhat is None:
        yhat = x @ beta
    rss = float((y - yhat) @ (y - yhat))
    tss = float((y - y.mean()) @ (y - y.mean()))
    return {
        "log_likelihood": log_likelihood,
        "aic": -2 * log_likelihood + 2 * parameter_count,
        "bic": -2 * log_likelihood + parameter_count * math.log(n),
        "sigma2": sigma2,
        "pseudo_r2": 1 - rss / tss if tss > 0 else float("nan"),
    }


def _bounds(w: np.ndarray) -> tuple[float, float]:
    eigenvalues = np.real(np.linalg.eigvals(w))
    positive = eigenvalues[eigenvalues > 1e-8]
    negative = eigenvalues[eigenvalues < -1e-8]
    upper = 0.999 if not len(positive) else min(0.999, 1 / positive.max())
    lower = -0.999 if not len(negative) else max(-0.999, 1 / negative.min())
    if lower >= upper - 0.01:
        lower, upper = -0.9, 0.9
    return lower + 1e-4, upper - 1e-4


def _golden_max(
    function, lower: float, upper: float, tolerance: float = 1e-7
) -> tuple[float, float]:
    ratio = (math.sqrt(5) - 1) / 2
    c = upper - ratio * (upper - lower)
    d = lower + ratio * (upper - lower)
    fc, fd = function(c), function(d)
    for _ in range(80):
        if upper - lower < tolerance:
            break
        if fc < fd:
            lower, c, fc = c, d, fd
            d = lower + ratio * (upper - lower)
            fd = function(d)
        else:
            upper, d, fd = d, c, fc
            c = upper - ratio * (upper - lower)
            fc = function(c)
    best = (lower + upper) / 2
    return best, function(best)


def _profile_spatial(
    y: np.ndarray,
    x: np.ndarray,
    w: np.ndarray,
    parameter: str,
) -> tuple[float, dict[str, object]]:
    n = len(y)
    identity = np.eye(n)

    def profile(value: float) -> float:
        a = identity - value * w
        sign, logdet = np.linalg.slogdet(a)
        if sign <= 0:
            return -math.inf
        if parameter == "rho":
            transformed_y, transformed_x = a @ y, x
        else:
            transformed_y, transformed_x = a @ y, a @ x
        beta, *_ = np.linalg.lstsq(transformed_x, transformed_y, rcond=None)
        residual = transformed_y - transformed_x @ beta
        sigma2 = float(residual @ residual) / n
        if sigma2 <= 0:
            return -math.inf
        return -0.5 * n * (math.log(2 * math.pi) + math.log(sigma2) + 1) + logdet

    lower, upper = _bounds(w)
    grid = np.linspace(lower, upper, 61)
    values = [profile(float(value)) for value in grid]
    best_index = int(np.argmax(values))
    left = float(grid[max(0, best_index - 1)])
    right = float(grid[min(len(grid) - 1, best_index + 1)])
    value, _ = _golden_max(profile, left, right)

    a = identity - value * w
    _, logdet = np.linalg.slogdet(a)
    if parameter == "rho":
        transformed_y, transformed_x = a @ y, x
    else:
        transformed_y, transformed_x = a @ y, a @ x
    beta, *_ = np.linalg.lstsq(transformed_x, transformed_y, rcond=None)
    residual = transformed_y - transformed_x @ beta
    sigma2 = float(residual @ residual) / n
    log_likelihood = -0.5 * n * (math.log(2 * math.pi) + math.log(sigma2) + 1) + logdet
    return value, {
        "beta": beta,
        "sigma2": sigma2,
        "log_likelihood": log_likelihood,
        "log_det": float(logdet),
        "residual": residual,
    }


def _numerical_hessian(function, parameters: np.ndarray, relative_step: float = 1e-4) -> np.ndarray:
    size = len(parameters)
    steps = relative_step * np.maximum(1.0, np.abs(parameters))
    hessian = np.zeros((size, size), dtype=float)
    base = function(parameters)
    for i in range(size):
        step_i = steps[i]
        plus_i = parameters.copy()
        minus_i = parameters.copy()
        plus_i[i] += step_i
        minus_i[i] -= step_i
        hessian[i, i] = (function(plus_i) - 2 * base + function(minus_i)) / (step_i ** 2)
        for j in range(i):
            step_j = steps[j]
            plus_plus = parameters.copy()
            plus_minus = parameters.copy()
            minus_plus = parameters.copy()
            minus_minus = parameters.copy()
            plus_plus[i] += step_i
            plus_plus[j] += step_j
            plus_minus[i] += step_i
            plus_minus[j] -= step_j
            minus_plus[i] -= step_i
            minus_plus[j] += step_j
            minus_minus[i] -= step_i
            minus_minus[j] -= step_j
            value = (
                function(plus_plus)
                - function(plus_minus)
                - function(minus_plus)
                + function(minus_minus)
            ) / (4 * step_i * step_j)
            hessian[i, j] = value
            hessian[j, i] = value
    return hessian


def _standard_errors(
    full_loglike,
    beta: np.ndarray,
    spatial_parameter: float,
    sigma2: float,
) -> tuple[np.ndarray, float | None]:
    parameters = np.concatenate([beta, [spatial_parameter], [math.log(sigma2)]])
    hessian = _numerical_hessian(full_loglike, parameters)
    try:
        covariance = np.linalg.pinv(-hessian)
    except np.linalg.LinAlgError:
        return np.full(len(beta), float("nan")), None
    variances = np.maximum(np.diag(covariance), 0)
    beta_se = np.sqrt(variances[: len(beta)])
    parameter_se = math.sqrt(variances[len(beta)]) if variances[len(beta)] > 0 else None
    return beta_se, parameter_se


def fit_ols(
    y: np.ndarray, x: np.ndarray, w: np.ndarray, terms: Sequence[str]
) -> dict[str, object]:
    n = len(y)
    beta, residual, sigma2, log_likelihood = _ols(y, x)
    rss = float(residual @ residual)
    covariance = sigma2 * np.linalg.pinv(x.T @ x)
    standard_errors = np.sqrt(np.maximum(np.diag(covariance), 0))
    z_values = beta / standard_errors
    p_values = [_normal_p(float(value)) for value in z_values]
    moran_i, _, moran_z = _moran_i(residual.tolist(), w.tolist())
    fits = _fit_stats(y, x, beta, log_likelihood, sigma2, x.shape[1] + 1)
    fits.update(
        {
            "model": "OLS",
            "spatial_parameter": None,
            "residual_moran_i": moran_i,
            "residual_moran_p": _normal_p(moran_z),
        }
    )
    coefficients = []
    for index, term in enumerate(terms):
        coefficients.append(
            {
                "model": "OLS",
                "term": term,
                "estimate": float(beta[index]),
                "std_error": float(standard_errors[index]),
                "z": float(z_values[index]),
                "p_value": p_values[index],
                "direct_effect": float(beta[index]),
                "indirect_effect": 0.0,
                "total_effect": float(beta[index]),
            }
        )
    return {
        "model": "OLS",
        "coefficients": coefficients,
        "fits": fits,
        "beta": beta,
        "residual": residual,
        "sigma2": sigma2,
        "log_likelihood": log_likelihood,
        "spatial_parameter": None,
        "yhat": x @ beta,
        "structural_residual": residual,
        "rss": rss,
    }


def fit_slm(
    y: np.ndarray, x: np.ndarray, w: np.ndarray, terms: Sequence[str]
) -> dict[str, object]:
    n = len(y)
    rho, profile = _profile_spatial(y, x, w, "rho")
    beta = np.asarray(profile["beta"])
    sigma2 = float(profile["sigma2"])
    log_likelihood = float(profile["log_likelihood"])
    identity = np.eye(n)
    a = identity - rho * w

    def full_loglike(parameters: np.ndarray) -> float:
        current_beta = parameters[: len(beta)]
        current_rho = parameters[len(beta)]
        current_log_sigma2 = parameters[len(beta) + 1]
        current_a = identity - current_rho * w
        sign, logdet = np.linalg.slogdet(current_a)
        if sign <= 0:
            return -1e100
        current_sigma2 = math.exp(current_log_sigma2)
        residual = current_a @ y - x @ current_beta
        return (
            -0.5 * n * (math.log(2 * math.pi) + current_log_sigma2)
            - float(residual @ residual) / (2 * current_sigma2)
            + logdet
        )

    beta_se, rho_se = _standard_errors(full_loglike, beta, rho, sigma2)
    z_values = beta / beta_se
    p_values = [_normal_p(float(value)) for value in z_values]
    yhat = np.linalg.solve(a, x @ beta)
    structural_residual = a @ y - x @ beta
    moran_i, _, moran_z = _moran_i(structural_residual.tolist(), w.tolist())
    fits = _fit_stats(y, x, beta, log_likelihood, sigma2, x.shape[1] + 2, yhat=yhat)
    fits.update(
        {
            "model": "SLM",
            "spatial_parameter": rho,
            "residual_moran_i": moran_i,
            "residual_moran_p": _normal_p(moran_z),
        }
    )
    coefficients = []
    multiplier = rho / (1 - rho) if abs(1 - rho) > 1e-9 else float("nan")
    for index, term in enumerate(terms):
        coefficients.append(
            {
                "model": "SLM",
                "term": term,
                "estimate": float(beta[index]),
                "std_error": float(beta_se[index]),
                "z": float(z_values[index]),
                "p_value": p_values[index],
                "direct_effect": float(beta[index]),
                "indirect_effect": float(beta[index] * multiplier),
                "total_effect": float(beta[index] * (1 + multiplier)),
            }
        )
    coefficients.append(
        {
            "model": "SLM",
            "term": "rho",
            "estimate": float(rho),
            "std_error": rho_se,
            "z": float(rho / rho_se) if rho_se else None,
            "p_value": _normal_p(float(rho / rho_se)) if rho_se else None,
            "direct_effect": None,
            "indirect_effect": None,
            "total_effect": None,
        }
    )
    return {
        "model": "SLM",
        "coefficients": coefficients,
        "fits": fits,
        "beta": beta,
        "residual": structural_residual,
        "sigma2": sigma2,
        "log_likelihood": log_likelihood,
        "spatial_parameter": rho,
        "yhat": yhat,
        "structural_residual": structural_residual,
        "rss": float((y - yhat) @ (y - yhat)),
    }


def fit_sem(
    y: np.ndarray, x: np.ndarray, w: np.ndarray, terms: Sequence[str]
) -> dict[str, object]:
    n = len(y)
    lam, profile = _profile_spatial(y, x, w, "lambda")
    beta = np.asarray(profile["beta"])
    sigma2 = float(profile["sigma2"])
    log_likelihood = float(profile["log_likelihood"])
    identity = np.eye(n)

    def full_loglike(parameters: np.ndarray) -> float:
        current_beta = parameters[: len(beta)]
        current_lambda = parameters[len(beta)]
        current_log_sigma2 = parameters[len(beta) + 1]
        current_a = identity - current_lambda * w
        sign, logdet = np.linalg.slogdet(current_a)
        if sign <= 0:
            return -1e100
        current_sigma2 = math.exp(current_log_sigma2)
        residual = current_a @ (y - x @ current_beta)
        return (
            -0.5 * n * (math.log(2 * math.pi) + current_log_sigma2)
            - float(residual @ residual) / (2 * current_sigma2)
            + logdet
        )

    beta_se, lambda_se = _standard_errors(full_loglike, beta, lam, sigma2)
    z_values = beta / beta_se
    p_values = [_normal_p(float(value)) for value in z_values]
    a = identity - lam * w
    yhat = x @ beta
    structural_residual = a @ (y - yhat)
    moran_i, _, moran_z = _moran_i(structural_residual.tolist(), w.tolist())
    fits = _fit_stats(y, x, beta, log_likelihood, sigma2, x.shape[1] + 2, yhat=yhat)
    fits.update(
        {
            "model": "SEM",
            "spatial_parameter": lam,
            "residual_moran_i": moran_i,
            "residual_moran_p": _normal_p(moran_z),
        }
    )
    coefficients = []
    for index, term in enumerate(terms):
        coefficients.append(
            {
                "model": "SEM",
                "term": term,
                "estimate": float(beta[index]),
                "std_error": float(beta_se[index]),
                "z": float(z_values[index]),
                "p_value": p_values[index],
                "direct_effect": float(beta[index]),
                "indirect_effect": 0.0,
                "total_effect": float(beta[index]),
            }
        )
    coefficients.append(
        {
            "model": "SEM",
            "term": "lambda",
            "estimate": float(lam),
            "std_error": lambda_se,
            "z": float(lam / lambda_se) if lambda_se else None,
            "p_value": _normal_p(float(lam / lambda_se)) if lambda_se else None,
            "direct_effect": None,
            "indirect_effect": None,
            "total_effect": None,
        }
    )
    return {
        "model": "SEM",
        "coefficients": coefficients,
        "fits": fits,
        "beta": beta,
        "residual": structural_residual,
        "sigma2": sigma2,
        "log_likelihood": log_likelihood,
        "spatial_parameter": lam,
        "yhat": yhat,
        "structural_residual": structural_residual,
        "rss": float((y - yhat) @ (y - yhat)),
    }


def write_outputs(
    observations: Sequence[Observation],
    y: np.ndarray,
    results: Sequence[dict[str, object]],
    terms: Sequence[str],
    summary: dict[str, object],
    output_dir: Path,
    prefix: str,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    coeff_path = output_dir / f"{prefix}_coefficients.csv"
    fits_path = output_dir / f"{prefix}_fits.csv"
    fitted_path = output_dir / f"{prefix}_fitted.csv"
    summary_path = output_dir / f"{prefix}_summary.json"

    with coeff_path.open("w", encoding="utf-8-sig", newline="") as handle:
        fields = [
            "model", "term", "estimate", "std_error", "z", "p_value",
            "direct_effect", "indirect_effect", "total_effect",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for result in results:
            for row in result["coefficients"]:
                writer.writerow(row)

    with fits_path.open("w", encoding="utf-8-sig", newline="") as handle:
        fields = [
            "model", "log_likelihood", "aic", "bic", "sigma2", "pseudo_r2",
            "spatial_parameter", "residual_moran_i", "residual_moran_p",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for result in results:
            writer.writerow(result["fits"])

    with fitted_path.open("w", encoding="utf-8-sig", newline="") as handle:
        fields = ["community_id", "name", "observed"]
        for result in results:
            model = str(result["model"])
            fields.extend([f"{model}_fitted", f"{model}_residual"])
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for index, observation in enumerate(observations):
            row = {
                "community_id": observation.community_id,
                "name": observation.name,
                "observed": f"{y[index]:.6f}",
            }
            for result in results:
                model = str(result["model"])
                yhat = np.asarray(result["yhat"])
                row[f"{model}_fitted"] = f"{yhat[index]:.6f}"
                row[f"{model}_residual"] = f"{y[index] - yhat[index]:.6f}"
            writer.writerow(row)

    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, allow_nan=False)
    return {
        "coefficients": coeff_path,
        "fits": fits_path,
        "fitted": fitted_path,
        "summary": summary_path,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="OLS / 空间滞后（SLM）/ 空间误差（SEM）教学示例",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--access", type=Path, default=DEFAULT_ACCESS, help="可达性 CSV")
    parser.add_argument("--outcomes", type=Path, default=DEFAULT_OUTCOMES, help="结果变量 CSV")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT, help="输出目录")
    parser.add_argument("--prefix", default=DEFAULT_PREFIX, help="输出文件前缀")
    parser.add_argument("--k", type=int, default=6, help="k 近邻空间权重")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        observations = read_access(args.access)
        outcomes = read_outcomes(args.outcomes)
        y, x, terms, matched = build_data(observations, outcomes)
        weights_list, _ = build_knn_weights(matched, args.k)
        w = np.asarray(weights_list, dtype=float)
        if not np.isfinite(w).all():
            raise ValueError("空间权重包含非法数值")
    except ValueError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 2

    ols = fit_ols(y, x, w, terms)
    slm = fit_slm(y, x, w, terms)
    sem = fit_sem(y, x, w, terms)
    results = [ols, slm, sem]

    lr_tests = {}
    for result in (slm, sem):
        lr = 2 * (float(result["log_likelihood"]) - float(ols["log_likelihood"]))
        lr_tests[f"{result['model']}_vs_OLS"] = {
            "lr_statistic": lr,
            "df": 1,
            "p_value": _chi2_1_p(lr),
        }

    summary = {
        "method": "OLS / SLM / SEM (maximum likelihood)",
        "dependent_variable": "outpatient_utilization",
        "covariates": terms[1:],
        "n": len(y),
        "k": args.k,
        "weights": {"type": "k nearest neighbors", "row_standardized": True},
        "fits": [result["fits"] for result in results],
        "lr_tests": lr_tests,
        "note": (
            "标准误由数值 Hessian 近似，ρ/λ 的检验在边界上应谨慎解读；"
            "这是教学实现，正式论文建议用 spreg 或 GeoDa 复核，并做不同 k 的敏感性分析。"
        ),
    }
    paths = write_outputs(matched, y, results, terms, summary, args.output_dir, args.prefix)

    print("模型比较：")
    for result in results:
        fits = result["fits"]
        print(
            f"  {fits['model']:4s}  logL={fits['log_likelihood']:.2f}  "
            f"AIC={fits['aic']:.2f}  pseudo-R²={fits['pseudo_r2']:.3f}  "
            f"resid.Moran I={fits['residual_moran_i']:.3f}"
        )
    for name, test in lr_tests.items():
        print(f"  LR {name}: {test['lr_statistic']:.2f}, p={test['p_value']:.3f}")
    print("核心系数：")
    for result in results:
        for row in result["coefficients"]:
            if row["term"] in ("access_per_1000", "rho", "lambda"):
                estimate = row["estimate"]
                p_value = row["p_value"]
                print(
                    f"  {result['model']:4s} {row['term']:16s} "
                    f"{float(estimate): .4f}  p={float(p_value):.3f}"
                    if p_value is not None
                    else f"  {result['model']:4s} {row['term']:16s} {float(estimate): .4f}"
                )
    print("输出文件：")
    for name, path in paths.items():
        print(f"  {name:12s} {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
