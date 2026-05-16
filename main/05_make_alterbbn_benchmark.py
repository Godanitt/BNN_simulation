#!/usr/bin/env python3
"""Build the AlterBBN benchmark tables and combined comparison figures.

This script is intentionally independent of pynucastro.  It consumes the CSVs
produced by the pynucastro notebooks plus the AlterBBN numbers copied from
`./stand_cosmo.x 5`, and writes the combined benchmark used by the Typst memory.
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
    {
        "observable": "Yp_He4",
        "display": r"$Y_p$",
        "alterbbn_column": "Yp",
        "low": 2.521e-01,
        "cent": 2.525e-01,
        "high": 2.525e-01,
        "value": 2.525e-01,
        "sigma": 3.280e-04,
        "source": "AlterBBN stand_cosmo.x 5",
        "notes": "Primordial helium-4 mass fraction.",
    },
    {
        "observable": "D/H",
        "display": r"D/H",
        "alterbbn_column": "H2/H",
        "low": 2.766e-05,
        "cent": 2.674e-05,
        "high": 2.558e-05,
        "value": 2.674e-05,
        "sigma": 5.462e-07,
        "source": "AlterBBN stand_cosmo.x 5",
        "notes": "Deuterium-to-hydrogen number ratio.",
    },
    {
        "observable": "He3/H",
        "display": r"$^3$He/H",
        "alterbbn_column": "He3/H",
        "low": 1.056e-05,
        "cent": 1.067e-05,
        "high": 1.068e-05,
        "value": 1.067e-05,
        "sigma": 1.730e-07,
        "source": "AlterBBN stand_cosmo.x 5",
        "notes": "Helium-3-to-hydrogen number ratio.",
    },
    {
        "observable": "Li7/H",
        "display": r"$^7$Li/H",
        "alterbbn_column": "Li7/H",
        "low": 4.582e-10,
        "cent": 4.972e-10,
        "high": 5.497e-10,
        "value": 4.972e-10,
        "sigma": 3.616e-11,
        "source": "AlterBBN stand_cosmo.x 5",
        "notes": "Final mass-7 lithium abundance reported by AlterBBN. This is the quantity compared with observed Li7/H.",
    },
    {
        "observable": "Li6/H",
        "display": r"$^6$Li/H",
        "alterbbn_column": "Li6/H",
        "low": 1.903e-15,
        "cent": 1.211e-14,
        "high": 3.843e-14,
        "value": 1.211e-14,
        "sigma": 1.219e-14,
        "source": "AlterBBN stand_cosmo.x 5",
        "notes": "Lithium-6-to-hydrogen number ratio. Not used in the main observational comparison.",
    },
    {
        "observable": "Be7/H",
        "display": r"$^7$Be/H",
        "alterbbn_column": "Be7/H",
        "low": 4.268e-10,
        "cent": 4.652e-10,
        "high": 5.179e-10,
        "value": 4.652e-10,
        "sigma": 3.548e-11,
        "source": "AlterBBN stand_cosmo.x 5",
        "notes": "Beryllium-7-to-hydrogen number ratio at the end of the BBN calculation.",
    },
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
# 2. Read pynucastro final values and observational references.
# -----------------------------------------------------------------------------
finals_path = DATA_DIR / "bbn_final_abundances_all_models.csv"
obs_path = DATA_DIR / "observational_abundances.csv"
if not finals_path.exists():
    raise FileNotFoundError(f"Missing {finals_path}. Run notebook 02 first.")
if not obs_path.exists():
    raise FileNotFoundError(f"Missing {obs_path}. Run notebook 00 first.")

finals = pd.read_csv(finals_path)
obs = pd.read_csv(obs_path)
pyn = finals.iloc[-1]

# Mapping between the common observable and the column used by the pynucastro table.
# For Li7, the fair comparison is the post-BBN mass-7 abundance.  In the pynucastro
# table it is Li7 + Be7, while AlterBBN already reports the final Li7/H abundance.
pyn_map = {
    "D/H": "D/H",
    "He3/H": "He3/H",
    "Yp_He4": "X_He4",
    "Li7/H": "Li7_plus_Be7_over_H",
}

ordered_obs = ["Yp_He4", "D/H", "He3/H", "Li7/H"]
compare_rows = []
for observable in ordered_obs:
    compare_rows.append({
        "observable": observable,
        "kind": "pynucastro",
        "value": float(pyn[pyn_map[observable]]),
        "sigma": np.nan,
        "source": "pynucastro light network",
        "notes": f"column={pyn_map[observable]}",
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
# Backwards-compatible name read by the memory.
compare.to_csv(DATA_DIR / "comparison_models_vs_observed.csv", index=False)

# Compact final benchmark with one row per method and one column per observable.
combined = pd.DataFrame([
    {
        "model": "pynucastro",
        "Yp_He4": float(pyn["X_He4"]),
        "D/H": float(pyn["D/H"]),
        "He3/H": float(pyn["He3/H"]),
        "Li7/H": float(pyn["Li7_plus_Be7_over_H"]),
        "Li6/H": np.nan,
        "Be7/H": float(pyn["Be7/H"]),
        "source": "pynucastro light network",
    },
    {
        "model": "AlterBBN",
        "Yp_He4": float(alter.loc[alter.observable == "Yp_He4", "value"].iloc[0]),
        "D/H": float(alter.loc[alter.observable == "D/H", "value"].iloc[0]),
        "He3/H": float(alter.loc[alter.observable == "He3/H", "value"].iloc[0]),
        "Li7/H": float(alter.loc[alter.observable == "Li7/H", "value"].iloc[0]),
        "Li6/H": float(alter.loc[alter.observable == "Li6/H", "value"].iloc[0]),
        "Be7/H": float(alter.loc[alter.observable == "Be7/H", "value"].iloc[0]),
        "source": "AlterBBN stand_cosmo.x 5",
    },
])
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

# Helper used for bar placement.
def _bar_positions(n_groups: int, n_kinds: int, width: float = 0.25):
    x = np.arange(n_groups)
    offsets = (np.arange(n_kinds) - (n_kinds - 1) / 2.0) * width
    return x, offsets

# Final ratios: D, He3, Li7 in log scale.
ratio_obs = ["D/H", "He3/H", "Li7/H"]
ratio_labels = [r"D/H", r"$^3$He/H", r"$^7$Li/H"]
kinds = ["pynucastro", "AlterBBN", "observed"]
fig, ax = plt.subplots(figsize=(7.4, 4.8))
x, offsets = _bar_positions(len(ratio_obs), len(kinds), width=0.24)
for j, kind in enumerate(kinds):
    vals = []
    errs = []
    for ob in ratio_obs:
        sub = compare[(compare.observable == ob) & (compare.kind == kind)]
        if len(sub) == 0:
            vals.append(np.nan)
            errs.append(0.0)
        else:
            vals.append(float(sub.value.iloc[0]))
            err = sub.sigma.iloc[0]
            errs.append(0.0 if pd.isna(err) else float(err))
    ax.bar(x + offsets[j], vals, width=0.23, label=kind, yerr=errs, capsize=2)
ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels(ratio_labels)
ax.set_ylabel("abundance ratio")
ax.set_title("Final abundance ratios")
ax.legend(frameon=True)
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_final_ratios_bar_models_vs_observed.pdf")
plt.close(fig)

# Helium-4 separate.
fig, ax = plt.subplots(figsize=(5.8, 4.2))
vals, errs = [], []
for kind in kinds:
    sub = compare[(compare.observable == "Yp_He4") & (compare.kind == kind)]
    vals.append(float(sub.value.iloc[0]))
    err = sub.sigma.iloc[0]
    errs.append(0.0 if pd.isna(err) else float(err))
ax.bar(kinds, vals, yerr=errs, capsize=3)
ax.set_ylabel(r"$Y_p$ / $X(^4\mathrm{He})$")
ax.set_ylim(min(vals) - 0.006, max(vals) + 0.006)
ax.set_title("Final helium-4 mass fraction")
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_final_he4_bar_models_vs_observed.pdf")
plt.close(fig)

# Evolution ratios with AlterBBN and observations as horizontal references.
evol = pd.read_csv(DATA_DIR / "bbn_evolution_t0_50s.csv")
fig, ax = plt.subplots(figsize=(7.3, 4.8))
curve_map = {
    "D/H": ("D/H", r"pynucastro D/H"),
    "He3/H": ("He3/H", r"pynucastro $^3$He/H"),
    "Li7/H": ("Li7_plus_Be7_over_H", r"pynucastro mass-7/H"),
}
for ob, (col, lab) in curve_map.items():
    ax.loglog(evol["t_s"], evol[col].clip(lower=1e-40), label=lab)
    a = alter[alter.observable == ob].iloc[0]
    ax.axhline(a["value"], linestyle="--", linewidth=1.0, alpha=0.9, label=f"AlterBBN {ob}")
    obrow = compare[(compare.observable == ob) & (compare.kind == "observed")]
    if len(obrow) > 0:
        val = float(obrow.value.iloc[0]); sig = float(obrow.sigma.iloc[0])
        if val - sig > 0:
            ax.axhspan(val - sig, val + sig, alpha=0.10)
ax.set_xlabel("time [s]")
ax.set_ylabel("abundance ratio")
ax.set_title("pynucastro evolution vs AlterBBN and observations")
ax.legend(fontsize=7, ncol=2, frameon=True)
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_key_ratios_evolution_vs_observed.pdf")
plt.close(fig)

# He4 evolution with AlterBBN and observed band.
fig, ax = plt.subplots(figsize=(7.0, 4.5))
ax.semilogx(evol["t_s"], evol["X_He4"], label="pynucastro")
yp = float(alter.loc[alter.observable == "Yp_He4", "value"].iloc[0])
yp_sig = float(alter.loc[alter.observable == "Yp_He4", "sigma"].iloc[0])
ax.axhline(yp, linestyle="--", linewidth=1.2, label="AlterBBN")
ax.axhspan(yp - yp_sig, yp + yp_sig, alpha=0.10)
obs_yp = compare[(compare.observable == "Yp_He4") & (compare.kind == "observed")]
if len(obs_yp) > 0:
    val = float(obs_yp.value.iloc[0]); sig = float(obs_yp.sigma.iloc[0])
    ax.axhspan(val - sig, val + sig, alpha=0.12, label="observed")
ax.set_xlabel("time [s]")
ax.set_ylabel(r"$X(^4\mathrm{He})$")
ax.set_title("Helium-4 mass fraction")
ax.legend(frameon=True)
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_he4_evolution_vs_observed.pdf")
plt.close(fig)

# Lithium diagnostic.
fig, ax = plt.subplots(figsize=(6.6, 4.5))
li_kinds = ["pynucastro", "AlterBBN", "observed"]
li_vals, li_errs = [], []
for kind in li_kinds:
    sub = compare[(compare.observable == "Li7/H") & (compare.kind == kind)]
    li_vals.append(float(sub.value.iloc[0]))
    err = sub.sigma.iloc[0]
    li_errs.append(0.0 if pd.isna(err) else float(err))
ax.bar(li_kinds, li_vals, yerr=li_errs, capsize=3)
ax.set_yscale("log")
ax.set_ylabel(r"$^7$Li/H")
ax.set_title("Mass-7 / lithium diagnostic")
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_lithium7_problem_comparison.pdf")
plt.close(fig)

# Compact benchmark figure: one panel-equivalent figure but no subplots.
fig, ax = plt.subplots(figsize=(7.0, 4.4))
benchmark = compare[compare.observable.isin(["D/H", "He3/H", "Li7/H"])]
for kind in kinds:
    sub = benchmark[benchmark.kind == kind]
    vals = [float(sub[sub.observable == ob].value.iloc[0]) for ob in ratio_obs]
    ax.plot(ratio_labels, vals, marker="o", label=kind)
ax.set_yscale("log")
ax.set_ylabel("abundance ratio")
ax.set_title("Combined BBN benchmark")
ax.legend(frameon=True)
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_combined_benchmark_pynucastro_alterbbn_observed.pdf")
plt.close(fig)

print("Wrote AlterBBN benchmark and combined comparison outputs.")
print(DATA_DIR / "alterbbn_standard_results.csv")
print(DATA_DIR / "comparison_pynucastro_alterbbn_observed.csv")
