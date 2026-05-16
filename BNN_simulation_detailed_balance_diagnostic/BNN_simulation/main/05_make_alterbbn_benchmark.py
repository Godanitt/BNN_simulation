#!/usr/bin/env python3
"""Build AlterBBN and model benchmark tables/figures.

This script is intentionally independent of pynucastro. It consumes the CSVs
produced by the pynucastro notebooks plus the AlterBBN numbers copied from
`./stand_cosmo.x 5`.  If the optional detailed-balance test notebook has added
a `pynucastro_db` row to `bbn_final_abundances_all_models.csv`, that row is
included automatically in all final benchmark tables and bar plots.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
DATA_DIR = BASE / "data"
FIG_DIR = BASE / "figures"
TAB_DIR = BASE / "tables"
for d in (DATA_DIR, FIG_DIR, TAB_DIR):
    d.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# 1. AlterBBN output copied from: ./stand_cosmo.x 5
# -----------------------------------------------------------------------------
ALTER_ROWS = [
    {"observable": "Yp_He4", "display": r"$Y_p$", "alterbbn_column": "Yp", "low": 2.521e-01, "cent": 2.525e-01, "high": 2.525e-01, "value": 2.525e-01, "sigma": 3.280e-04, "source": "AlterBBN stand_cosmo.x 5", "notes": "Primordial helium-4 mass fraction."},
    {"observable": "D/H", "display": r"D/H", "alterbbn_column": "H2/H", "low": 2.766e-05, "cent": 2.674e-05, "high": 2.558e-05, "value": 2.674e-05, "sigma": 5.462e-07, "source": "AlterBBN stand_cosmo.x 5", "notes": "Deuterium-to-hydrogen number ratio."},
    {"observable": "He3/H", "display": r"$^3$He/H", "alterbbn_column": "He3/H", "low": 1.056e-05, "cent": 1.067e-05, "high": 1.068e-05, "value": 1.067e-05, "sigma": 1.730e-07, "source": "AlterBBN stand_cosmo.x 5", "notes": "Helium-3-to-hydrogen number ratio."},
    {"observable": "Li7/H", "display": r"$^7$Li/H", "alterbbn_column": "Li7/H", "low": 4.582e-10, "cent": 4.972e-10, "high": 5.497e-10, "value": 4.972e-10, "sigma": 3.616e-11, "source": "AlterBBN stand_cosmo.x 5", "notes": "Final mass-7 lithium abundance reported by AlterBBN. This is the quantity compared with observed Li7/H."},
    {"observable": "Li6/H", "display": r"$^6$Li/H", "alterbbn_column": "Li6/H", "low": 1.903e-15, "cent": 1.211e-14, "high": 3.843e-14, "value": 1.211e-14, "sigma": 1.219e-14, "source": "AlterBBN stand_cosmo.x 5", "notes": "Lithium-6-to-hydrogen number ratio. Not used in the main observational comparison."},
    {"observable": "Be7/H", "display": r"$^7$Be/H", "alterbbn_column": "Be7/H", "low": 4.268e-10, "cent": 4.652e-10, "high": 5.179e-10, "value": 4.652e-10, "sigma": 3.548e-11, "source": "AlterBBN stand_cosmo.x 5", "notes": "Beryllium-7-to-hydrogen number ratio at the end of the BBN calculation."},
]

alter = pd.DataFrame(ALTER_ROWS)
alter.to_csv(DATA_DIR / "alterbbn_standard_results.csv", index=False)
alter.to_csv(TAB_DIR / "table_alterbbn_standard_results.csv", index=False)

corr_labels = ["Yp_He4", "D/H", "He3/H", "Li7/H", "Li6/H", "Be7/H"]
corr_values = np.array([
    [1.000000, 0.031046, -0.016455, -0.000513, -0.067454, -0.005372],
    [0.031046, 1.000000, -0.136597, -0.518948, 0.241945, -0.539147],
    [-0.016455, -0.136597, 1.000000, 0.106901, -0.028359, 0.106707],
    [-0.000513, -0.518948, 0.106901, 1.000000, -0.070280, 0.995547],
    [-0.067454, 0.241945, -0.028359, -0.070280, 1.000000, -0.067654],
    [-0.005372, -0.539147, 0.106707, 0.995547, -0.067654, 1.000000],
])
corr = pd.DataFrame(corr_values, columns=corr_labels)
corr.insert(0, "observable", corr_labels)
corr.to_csv(DATA_DIR / "alterbbn_correlation_matrix.csv", index=False)

# -----------------------------------------------------------------------------
# 2. Read model final values and observational references.
# -----------------------------------------------------------------------------
finals_path = DATA_DIR / "bbn_final_abundances_all_models.csv"
obs_path = DATA_DIR / "observational_abundances.csv"
if not finals_path.exists():
    raise FileNotFoundError(f"Missing {finals_path}. Run notebook 02 first.")
if not obs_path.exists():
    raise FileNotFoundError(f"Missing {obs_path}. Run notebook 00 first.")

finals = pd.read_csv(finals_path)
obs = pd.read_csv(obs_path)

MODEL_ORDER = ["pynucastro", "pynucastro_db"]
def model_sort_key(m: str) -> tuple[int, str]:
    return (MODEL_ORDER.index(m) if m in MODEL_ORDER else 100, m)

finals = finals.sort_values("model", key=lambda s: s.map(model_sort_key)).reset_index(drop=True)
pyn_map = {
    "D/H": "D/H",
    "He3/H": "He3/H",
    "Yp_He4": "X_He4",
    # For model networks, compare post-BBN mass-7 as Li7 + Be7.
    "Li7/H": "Li7_plus_Be7_over_H",
}
ordered_obs = ["Yp_He4", "D/H", "He3/H", "Li7/H"]

compare_rows: list[dict] = []
for observable in ordered_obs:
    for _, row in finals.iterrows():
        col = pyn_map[observable]
        if col not in row or pd.isna(row[col]):
            continue
        label = str(row["model"])
        compare_rows.append({
            "observable": observable,
            "kind": label,
            "value": float(row[col]),
            "sigma": np.nan,
            "source": f"{label} light network",
            "notes": f"column={col}",
        })
    a = alter.loc[alter["observable"] == observable].iloc[0]
    compare_rows.append({
        "observable": observable,
        "kind": "AlterBBN",
        "value": float(a["value"]),
        "sigma": float(a["sigma"]),
        "source": "AlterBBN stand_cosmo.x 5",
        "notes": str(a["notes"]),
    })
    matching_obs = obs[(obs["observable"] == observable) | (obs["model_column"] == pyn_map[observable])]
    if len(matching_obs) > 0:
        ob = matching_obs.iloc[0]
        compare_rows.append({
            "observable": observable,
            "kind": "observed",
            "value": float(ob["value"]),
            "sigma": float(ob["sigma"]),
            "source": str(ob["source_key"]),
            "notes": "observational reference used in the memory",
        })

compare = pd.DataFrame(compare_rows)
compare.to_csv(DATA_DIR / "comparison_pynucastro_alterbbn_observed.csv", index=False)
compare.to_csv(DATA_DIR / "comparison_models_vs_observed.csv", index=False)

# Compact final benchmark with one row per method and one column per observable.
combined_rows = []
for _, row in finals.iterrows():
    label = str(row["model"])
    combined_rows.append({
        "model": label,
        "Yp_He4": float(row["X_He4"]) if "X_He4" in row and pd.notna(row["X_He4"]) else np.nan,
        "D/H": float(row["D/H"]) if "D/H" in row and pd.notna(row["D/H"]) else np.nan,
        "He3/H": float(row["He3/H"]) if "He3/H" in row and pd.notna(row["He3/H"]) else np.nan,
        "Li7/H": float(row["Li7_plus_Be7_over_H"]) if "Li7_plus_Be7_over_H" in row and pd.notna(row["Li7_plus_Be7_over_H"]) else np.nan,
        "Li6/H": float(row["Li6/H"]) if "Li6/H" in row and pd.notna(row["Li6/H"]) else np.nan,
        "Be7/H": float(row["Be7/H"]) if "Be7/H" in row and pd.notna(row["Be7/H"]) else np.nan,
        "source": f"{label} light network",
    })

combined_rows.append({
    "model": "AlterBBN",
    "Yp_He4": float(alter.loc[alter.observable == "Yp_He4", "value"].iloc[0]),
    "D/H": float(alter.loc[alter.observable == "D/H", "value"].iloc[0]),
    "He3/H": float(alter.loc[alter.observable == "He3/H", "value"].iloc[0]),
    "Li7/H": float(alter.loc[alter.observable == "Li7/H", "value"].iloc[0]),
    "Li6/H": float(alter.loc[alter.observable == "Li6/H", "value"].iloc[0]),
    "Be7/H": float(alter.loc[alter.observable == "Be7/H", "value"].iloc[0]),
    "source": "AlterBBN stand_cosmo.x 5",
})
combined = pd.DataFrame(combined_rows)
combined.to_csv(DATA_DIR / "bbn_final_benchmark_combined.csv", index=False)

# -----------------------------------------------------------------------------
# 3. Formatting for Typst/LaTeX tables.
# -----------------------------------------------------------------------------
def fmt_sci(x: float, sig: int = 3) -> str:
    if pd.isna(x):
        return ""
    x = float(x)
    if x == 0.0:
        return "0"
    return f"{x:.{sig}e}"


def fmt_value(observable: str, x: float) -> str:
    if pd.isna(x):
        return ""
    if observable == "Yp_He4":
        return f"{float(x):.5f}"
    return fmt_sci(float(x), 3)

comparison_pretty = compare.copy()
comparison_pretty["value"] = [fmt_value(o, v) for o, v in zip(compare["observable"], compare["value"])]
comparison_pretty["sigma"] = [fmt_value(o, v) for o, v in zip(compare["observable"], compare["sigma"])]
comparison_pretty.to_csv(TAB_DIR / "table_comparison_models_vs_observed.csv", index=False)
comparison_pretty.to_latex(TAB_DIR / "table_comparison_models_vs_observed.tex", index=False, escape=False)

alter_pretty = alter[["observable", "value", "sigma", "low", "cent", "high", "source"]].copy()
for c in ["value", "sigma", "low", "cent", "high"]:
    alter_pretty[c] = [fmt_value(o, v) for o, v in zip(alter_pretty["observable"], alter_pretty[c])]
alter_pretty.to_csv(TAB_DIR / "table_alterbbn_standard_results.csv", index=False)
alter_pretty.to_latex(TAB_DIR / "table_alterbbn_standard_results.tex", index=False, escape=False)

combined_pretty = combined.copy()
for c in ["Yp_He4", "D/H", "He3/H", "Li7/H", "Li6/H", "Be7/H"]:
    combined_pretty[c] = [fmt_value(c, v) for v in combined[c]]
combined_pretty.to_csv(TAB_DIR / "table_bbn_final_benchmark_combined.csv", index=False)
combined_pretty.to_latex(TAB_DIR / "table_bbn_final_benchmark_combined.tex", index=False, escape=False)

# -----------------------------------------------------------------------------
# 4. Figures.
# -----------------------------------------------------------------------------
try:
    import scienceplots  # noqa: F401
    plt.style.use(["science", "grid"])
except Exception:
    pass

# Stable labels for plots.
def pretty_kind(k: str) -> str:
    return {
        "pynucastro": "pynucastro",
        "pynucastro_db": "pynucastro DB",
        "AlterBBN": "AlterBBN",
        "observed": "observed",
    }.get(k, k)

all_kinds = list(dict.fromkeys(compare["kind"].tolist()))
ratio_obs = ["D/H", "He3/H", "Li7/H"]
ratio_labels = [r"D/H", r"$^3$He/H", r"$^7$Li/H"]

def _plot_grouped_bars(obs_names, labels, outname, ylabel, title, yscale=None):
    sub = compare[compare["observable"].isin(obs_names)].copy()
    kinds = [k for k in all_kinds if k in set(sub["kind"])]
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    x = np.arange(len(obs_names))
    width = min(0.18, 0.78 / max(1, len(kinds)))
    offsets = (np.arange(len(kinds)) - (len(kinds) - 1) / 2.0) * width
    for j, kind in enumerate(kinds):
        vals, errs = [], []
        for lab in obs_names:
            rows = sub[(sub["observable"] == lab) & (sub["kind"] == kind)]
            if len(rows) == 0:
                vals.append(np.nan); errs.append(0.0)
            else:
                row = rows.iloc[0]
                vals.append(row["value"])
                errs.append(row["sigma"] if np.isfinite(row["sigma"]) else 0.0)
        ax.bar(x + offsets[j], vals, width=width, label=pretty_kind(kind), yerr=errs, capsize=2 if kind in {"AlterBBN", "observed"} else 0)
    if yscale:
        ax.set_yscale(yscale)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / outname)
    plt.close(fig)

_plot_grouped_bars(
    ratio_obs,
    ratio_labels,
    "fig_final_ratios_bar_models_vs_observed.pdf",
    "abundance ratio",
    "Final abundance ratios",
    yscale="log",
)
_plot_grouped_bars(
    ["Yp_He4"],
    [r"$Y_p$"],
    "fig_final_he4_bar_models_vs_observed.pdf",
    r"$X(^4\mathrm{He})$",
    r"Final helium-4 mass fraction",
    yscale=None,
)

# Combined benchmark figure with ratios only.
fig, ax = plt.subplots(figsize=(7.4, 4.8))
sub = compare[compare["observable"].isin(ratio_obs)].copy()
kinds = [k for k in all_kinds if k in set(sub["kind"])]
x = np.arange(len(ratio_obs))
width = min(0.18, 0.78 / max(1, len(kinds)))
offsets = (np.arange(len(kinds)) - (len(kinds) - 1) / 2.0) * width
for j, kind in enumerate(kinds):
    vals, errs = [], []
    for obs_name in ratio_obs:
        rows = sub[(sub["observable"] == obs_name) & (sub["kind"] == kind)]
        if len(rows) == 0:
            vals.append(np.nan); errs.append(0.0)
        else:
            row = rows.iloc[0]
            vals.append(row["value"])
            errs.append(row["sigma"] if np.isfinite(row["sigma"]) else 0.0)
    ax.bar(x + offsets[j], vals, width=width, label=pretty_kind(kind), yerr=errs, capsize=2 if kind in {"AlterBBN", "observed"} else 0)
ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels(ratio_labels)
ax.set_ylabel("number ratio")
ax.set_title("Final BBN benchmark")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_combined_benchmark_pynucastro_alterbbn_observed.pdf")
plt.close(fig)

# Evolution figures: include optional DB model if it exists in bbn_evolution_all_models.csv.
evo_path = DATA_DIR / "bbn_evolution_all_models.csv"
obs_map = {row["model_column"]: row for _, row in obs.iterrows()}
if evo_path.exists():
    evo = pd.read_csv(evo_path)
    plot_quantities = [
        ("D/H", r"D/H"),
        ("He3/H", r"$^3$He/H"),
        ("Li7_plus_Be7_over_H", r"$(^7$Li+$^7$Be)/H"),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for model_name, g in evo.groupby("model", sort=False):
        for col, label in plot_quantities:
            if col not in g.columns:
                continue
            linestyle = "-" if model_name == "pynucastro" else "--"
            ax.loglog(g["t_s"], g[col].clip(lower=1e-40), linestyle=linestyle, label=f"{pretty_kind(model_name)} {label}")
    for col, label in plot_quantities:
        if col in obs_map:
            val = obs_map[col]["value"]
            sig = obs_map[col]["sigma"]
            ax.axhspan(max(val-sig, 1e-40), val+sig, alpha=0.10)
            ax.axhline(val, lw=0.8, alpha=0.6)
    # AlterBBN points at t=1000s for visual reference.
    for obs_name, label in [("D/H", "AlterBBN D/H"), ("He3/H", r"AlterBBN $^3$He/H"), ("Li7/H", r"AlterBBN $^7$Li/H")]:
        val = alter.loc[alter.observable == obs_name, "value"].iloc[0]
        ax.scatter([1000.0], [val], marker="s", s=28, label=label)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("abundance ratio")
    ax.set_title("Key BBN abundance ratios")
    ax.legend(fontsize=6.5, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_key_ratios_evolution_vs_observed.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    for model_name, g in evo.groupby("model", sort=False):
        if "X_He4" in g.columns:
            ax.semilogx(g["t_s"], g["X_He4"], label=pretty_kind(model_name))
    he = obs[obs["model_column"] == "X_He4"].iloc[0]
    ax.axhspan(he["value"]-he["sigma"], he["value"]+he["sigma"], alpha=0.18, label="observed band")
    ax.axhline(he["value"], lw=1.0)
    ax.axhline(alter.loc[alter.observable == "Yp_He4", "value"].iloc[0], linestyle=":", lw=1.0, label="AlterBBN")
    ax.set_xlabel("time [s]")
    ax.set_ylabel(r"$X(^4\mathrm{He})$")
    ax.set_title(r"Helium-4 mass fraction")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_he4_evolution_vs_observed.pdf")
    plt.close(fig)

# Lithium standalone bar/evolution diagnostic.
fig, ax = plt.subplots(figsize=(7.2, 4.8))
sub_li = compare[compare["observable"] == "Li7/H"].copy()
li_kinds = [k for k in all_kinds if k in set(sub_li["kind"])]
x = np.arange(len(li_kinds))
vals = [sub_li[sub_li.kind == k]["value"].iloc[0] for k in li_kinds]
errs = [sub_li[sub_li.kind == k]["sigma"].iloc[0] if np.isfinite(sub_li[sub_li.kind == k]["sigma"].iloc[0]) else 0.0 for k in li_kinds]
ax.bar(x, vals, yerr=errs, capsize=2)
ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels([pretty_kind(k) for k in li_kinds], rotation=15)
ax.set_ylabel(r"$^7$Li/H or mass-7/H")
ax.set_title(r"Lithium-7 problem diagnostic")
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_lithium7_problem_comparison.pdf")
plt.close(fig)

print("Wrote AlterBBN/model benchmark tables and figures.")
print(combined.to_string(index=False))
