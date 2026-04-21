#!/usr/bin/env python3
"""
Run SHAP analysis for both variance-decomposition condition datasets using a
RandomForest regressor + TreeExplainer.

Two passes are produced per condition:
  1) full feature set (keeps layout/context features)
  2) decision-only feature set (drops row/col/category/layout proxies)

Outputs are preserved side-by-side via filename prefixing:
  - shap_*.{csv,pdf,txt}              (full feature set)
  - shap_decision_only_*.{csv,pdf,txt} (decision-only feature set)
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import shap
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_absolute_error
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path("SlabDesignFactors/results/sensitivity_variance_decomposition")
CONDITIONS = [
    ROOT / "yesdeflection_yesslabmin",
    ROOT / "yesdeflection_noslabmin",
]

# Keep runtime manageable while preserving stable global rankings.
MAX_EXPLAIN_ROWS = 250


def run_variant(
    condition_dir: Path,
    df: pd.DataFrame,
    prefix: str,
    drop_columns: list[str],
) -> tuple[float, float]:
    y = df["total_ec"].to_numpy()
    x = df.drop(columns=["total_ec"] + drop_columns, errors="ignore")
    x_enc = pd.get_dummies(x, drop_first=False).astype(np.float64)

    if x_enc.shape[1] == 0:
        raise ValueError(f"No features left for variant '{prefix}' in {condition_dir}")

    # A stable, moderate-capacity tree model that works with TreeExplainer.
    model = RandomForestRegressor(
        n_estimators=200,
        max_depth=None,
        min_samples_leaf=4,
        n_jobs=-1,
        random_state=42,
    )
    model.fit(x_enc, y)
    y_hat = model.predict(x_enc)

    # Subsample rows for SHAP computation speed.
    n = len(x_enc)
    explain_n = min(MAX_EXPLAIN_ROWS, n)
    explain_idx = np.random.default_rng(42).choice(n, size=explain_n, replace=False)
    x_explain = x_enc.iloc[explain_idx].copy()

    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(x_explain)
    if isinstance(shap_values, list):
        shap_values = shap_values[0]

    shap_df = pd.DataFrame(shap_values, columns=x_explain.columns)
    shap_df.to_csv(condition_dir / f"{prefix}_values.csv", index=False)

    importance = (
        pd.DataFrame(
            {
                "feature": x_explain.columns,
                "mean_abs_shap": np.abs(shap_values).mean(axis=0),
            }
        )
        .sort_values("mean_abs_shap", ascending=False)
        .reset_index(drop=True)
    )
    importance.to_csv(condition_dir / f"{prefix}_importance.csv", index=False)

    plt.figure()
    shap.summary_plot(
        shap_values,
        x_explain,
        plot_type="bar",
        show=False,
        max_display=25,
    )
    plt.tight_layout()
    plt.savefig(condition_dir / f"{prefix}_summary_bar.pdf")
    plt.close()

    plt.figure()
    shap.summary_plot(
        shap_values,
        x_explain,
        show=False,
        max_display=25,
    )
    plt.tight_layout()
    plt.savefig(condition_dir / f"{prefix}_summary_beeswarm.pdf")
    plt.close()

    r2 = r2_score(y, y_hat)
    mae = mean_absolute_error(y, y_hat)
    with open(condition_dir / f"{prefix}_fit_summary.txt", "w", encoding="utf-8") as f:
        f.write(f"variant={prefix}\n")
        f.write(f"rows_total={n}\n")
        f.write(f"rows_explained={explain_n}\n")
        f.write(f"features_onehot={x_enc.shape[1]}\n")
        f.write(f"dropped_columns={','.join(drop_columns) if drop_columns else '(none)'}\n")
        f.write("model=RandomForestRegressor\n")
        f.write(f"train_r2={r2:.6f}\n")
        f.write(f"train_mae={mae:.6f}\n")
        f.write("note=SHAP computed with TreeExplainer\n")

    return r2, mae


def run_condition(condition_dir: Path) -> None:
    csv_path = condition_dir / "master_design_matrix.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing master matrix: {csv_path}")

    df = pd.read_csv(csv_path)
    # Keep existing "full" output names so prior workflow still works.
    r2_full, mae_full = run_variant(
        condition_dir,
        df,
        prefix="shap",
        drop_columns=[],
    )

    # Decision-only: remove layout proxies so ranking reflects design decisions.
    # `layout` and `name` are both dropped defensively in case either appears.
    r2_dec, mae_dec = run_variant(
        condition_dir,
        df,
        prefix="shap_decision_only",
        drop_columns=["row", "col", "category", "layout", "name"],
    )

    print(f"[ok] {condition_dir}")
    print(f"     full          train R2={r2_full:.4f}, MAE={mae_full:.3f}")
    print(f"     decision-only train R2={r2_dec:.4f}, MAE={mae_dec:.3f}")


def main() -> None:
    for condition_dir in CONDITIONS:
        run_condition(condition_dir)
    print("[done] SHAP artifacts written for both conditions.")


if __name__ == "__main__":
    main()
