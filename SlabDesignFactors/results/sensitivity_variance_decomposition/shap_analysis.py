#!/usr/bin/env python3
'''SHAP analysis on the master design matrix exported by
`sensitivity_variance_decomposition.jl`.

Run separately for each condition from that condition's subfolder, e.g.

    cd SlabDesignFactors/results/sensitivity_variance_decomposition/yesdeflection_yesslabmin
    python ../shap_analysis.py master_design_matrix.csv

Requires `lightgbm`, `shap`, `pandas`, `matplotlib`. Outputs:
  - shap_summary_bar.pdf       global importance (mean |SHAP|)
  - shap_summary_beeswarm.pdf  per-sample SHAP with direction of effect
  - shap_values.csv            raw SHAP value matrix (rows × features)
'''
import sys, os
import pandas as pd
import numpy as np
import lightgbm as lgb
import shap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv_path = sys.argv[1] if len(sys.argv) > 1 else "master_design_matrix.csv"
outdir   = os.path.dirname(os.path.abspath(csv_path)) or "."

df = pd.read_csv(csv_path)
y  = df["total_ec"].values
X  = df.drop(columns=["total_ec"])

# one-hot encode categorical columns; LightGBM can also handle them natively
# via `categorical_feature`, but one-hot keeps SHAP interpretation simple.
X_enc = pd.get_dummies(X, drop_first=False)

model = lgb.LGBMRegressor(
    n_estimators=800,
    learning_rate=0.05,
    num_leaves=63,
    min_data_in_leaf=8,
    feature_fraction=0.9,
    bagging_fraction=0.9,
    bagging_freq=3,
    random_state=42,
)
model.fit(X_enc, y)

print(f"Train R²: {model.score(X_enc, y):.3f}")

explainer = shap.TreeExplainer(model)
shap_vals = explainer.shap_values(X_enc)

# Save raw SHAP matrix
shap_df = pd.DataFrame(shap_vals, columns=X_enc.columns)
shap_df.to_csv(os.path.join(outdir, "shap_values.csv"), index=False)

# Bar plot — global importance
plt.figure()
shap.summary_plot(shap_vals, X_enc, plot_type="bar", show=False, max_display=25)
plt.tight_layout()
plt.savefig(os.path.join(outdir, "shap_summary_bar.pdf"))
plt.close()

# Beeswarm — per-sample SHAP distribution
plt.figure()
shap.summary_plot(shap_vals, X_enc, show=False, max_display=25)
plt.tight_layout()
plt.savefig(os.path.join(outdir, "shap_summary_beeswarm.pdf"))
plt.close()

print("✓ SHAP outputs written to", outdir)
