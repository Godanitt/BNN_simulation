#!/usr/bin/env python3
"""Build and test an alternative pynucastro network with reverse rates.

The goal is diagnostic rather than cosmetic: the baseline BBN network burns
D/H far below the AlterBBN/observed scale.  The first physical suspect is the
incomplete treatment of the deuterium bottleneck, especially reverse channels
such as d -> n + p and the reverse rates associated with p(n,gamma)d,
d(p,gamma)He3, etc.

This script therefore creates a second network, `bbn_network_db.py`, in which
reverse rates are added by detailed balance whenever `pynucastro.DerivedRate`
can construct them.  If a local pynucastro version cannot derive some rates,
the script keeps a detailed diagnostic and falls back to ReacLib reverse rates
for coverage.  It then integrates the new network on exactly the same imposed
thermodynamic trajectory as the baseline model.

Outputs
-------
data/bbn_network_db_rates.csv
    All rates in the detailed-balance candidate network.
data/bbn_network_reverse_diagnostic.csv
    Coverage of the physically important reverse channels.
data/pynucastro_db_status.csv
    One-row status summary.
data/bbn_evolution_pynucastro_db.csv
    Evolution from 50 s to 1000 s, if integration succeeds.
data/bbn_final_abundances_pynucastro_db.csv
    Final row for the DB network, if integration succeeds.

If integration succeeds, the script also appends the DB model to:

data/bbn_evolution_all_models.csv
data/bbn_final_abundances_all_models.csv
data/bbn_final_compact.csv
data/bbn_snapshots_all_models.csv
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import importlib.util
import json
import math
import sys
import traceback
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

BASE = Path(__file__).resolve().parent
DATA_DIR = BASE / "data"
TAB_DIR = BASE / "tables"
FIG_DIR = BASE / "figures"
for d in (DATA_DIR, TAB_DIR, FIG_DIR):
    d.mkdir(parents=True, exist_ok=True)

NUCLEI_LIGHT = ["n", "p", "d", "t", "he3", "he4", "li6", "li7", "be7"]
# Important inverse channels for the deuterium bottleneck and mass-7 flow.
EXPECTED_REVERSE_CHANNELS = [
    ("d_to_n_p", ["d"], ["n", "p"], "photodisintegration counterpart of p(n,gamma)d"),
    ("he3_to_p_d", ["he3"], ["p", "d"], "reverse of d(p,gamma)he3"),
    ("t_to_n_d", ["t"], ["n", "d"], "reverse of d(n,gamma)t"),
    ("he4_to_p_t", ["he4"], ["p", "t"], "reverse of t(p,gamma)he4"),
    ("he4_to_n_he3", ["he4"], ["n", "he3"], "reverse of he3(n,gamma)he4"),
    ("be7_to_he4_he3", ["be7"], ["he4", "he3"], "reverse of he3(alpha,gamma)be7"),
    ("li7_to_he4_t", ["li7"], ["he4", "t"], "reverse of t(alpha,gamma)li7"),
]

ALIAS_TO_CANON = {
    "neutron": "n", "n": "n",
    "h1": "p", "p": "p", "proton": "p",
    "h2": "d", "d": "d",
    "h3": "t", "t": "t",
    "he3": "he3", "he4": "he4",
    "li6": "li6", "li7": "li7", "be7": "be7",
}


def canon_name(obj) -> str:
    s = str(obj).strip().lower().replace(" ", "")
    return ALIAS_TO_CANON.get(s, s)


def species_key(items: Iterable[object]) -> tuple[str, ...]:
    return tuple(sorted(canon_name(x) for x in items))


def rate_rows(rates: list[object]) -> pd.DataFrame:
    rows = []
    for i, r in enumerate(rates):
        rows.append({
            "index": i,
            "rate": str(r),
            "reactants": " + ".join(str(x) for x in getattr(r, "reactants", [])),
            "products": " + ".join(str(x) for x in getattr(r, "products", [])),
            "reactants_canon": " + ".join(species_key(getattr(r, "reactants", []))),
            "products_canon": " + ".join(species_key(getattr(r, "products", []))),
            "Q_MeV": getattr(r, "Q", np.nan),
            "fname": getattr(r, "fname", ""),
            "label": getattr(r, "label", ""),
            "weak": bool(getattr(r, "weak", False)),
            "derived_from_inverse": bool(getattr(r, "derived_from_inverse", False)),
            "rate_class": type(r).__name__,
        })
    return pd.DataFrame(rows)


def write_status(status: str, message: str, **extra) -> None:
    row = {"status": status, "message": message}
    row.update(extra)
    pd.DataFrame([row]).to_csv(DATA_DIR / "pynucastro_db_status.csv", index=False)


def export_not_run_diagnostics(reason: str) -> None:
    rows = []
    for name, reactants, products, note in EXPECTED_REVERSE_CHANNELS:
        rows.append({
            "channel": name,
            "reactants_expected": " + ".join(reactants),
            "products_expected": " + ".join(products),
            "found": False,
            "matched_rate": "",
            "Q_MeV": "",
            "rate_class": "",
            "derived_from_inverse": "",
            "note": note,
            "diagnostic": reason,
        })
    pd.DataFrame(rows).to_csv(DATA_DIR / "bbn_network_reverse_diagnostic.csv", index=False)


def import_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def unique_rates(rates: Iterable[object]) -> list[object]:
    """Remove exact duplicates while preserving order.

    We key primarily on reactants/products/class/source label.  This is loose
    enough to avoid repeated entries but conservative enough not to collapse
    physically distinct weak/reaclib rates with the same endpoints.
    """
    out = []
    seen = set()
    for r in rates:
        key = (
            species_key(getattr(r, "reactants", [])),
            species_key(getattr(r, "products", [])),
            getattr(r, "fname", ""),
            getattr(r, "label", ""),
            type(r).__name__,
            bool(getattr(r, "weak", False)),
            bool(getattr(r, "derived_from_inverse", False)),
        )
        if key not in seen:
            out.append(r)
            seen.add(key)
    return out


def build_db_network() -> tuple[object, list[object], pd.DataFrame, pd.DataFrame]:
    import pynucastro as pyna

    rl = pyna.ReacLibLibrary()

    # Baseline forward library.  This mirrors the philosophy of the original
    # notebook, but we deliberately include li6 if the local database supports it.
    lib_forward = rl.linking_nuclei(NUCLEI_LIGHT, with_reverse=False)
    forward_rates = list(getattr(lib_forward, "rates", []))

    # ReacLib library with reverse entries.  We keep it as fallback and for
    # comparison; some entries are already marked as derived_from_inverse.
    lib_with_reaclib_reverse = rl.linking_nuclei(NUCLEI_LIGHT, with_reverse=True)
    reaclib_reverse_rates = list(getattr(lib_with_reaclib_reverse, "rates", []))

    # Try explicit detailed-balance rates from the forward reactions.  This is
    # the most diagnostic part: if d -> n+p appears here, the bottleneck test is
    # meaningful.  DerivedRate rejects weak, tabular or already-derived rates;
    # every failure is stored in a CSV rather than silently ignored.
    derived_rates = []
    derived_attempt_rows = []
    DerivedRate = getattr(pyna, "DerivedRate", None)
    if DerivedRate is None:
        try:
            from pynucastro.rates import DerivedRate  # type: ignore
        except Exception:
            from pynucastro.rates.derived_rate import DerivedRate  # type: ignore

    for r in forward_rates:
        try:
            dr = DerivedRate(r, use_pf=False)
            # Keep only if all species remain inside our light network.
            all_nucs = [canon_name(x) for x in getattr(dr, "reactants", []) + getattr(dr, "products", [])]
            if all(x in set(NUCLEI_LIGHT) for x in all_nucs):
                derived_rates.append(dr)
                ok = True
                msg = "derived"
            else:
                ok = False
                msg = "outside selected nuclei"
        except Exception as exc:
            ok = False
            msg = f"{type(exc).__name__}: {exc}"
        derived_attempt_rows.append({
            "source_rate": str(r),
            "source_reactants": " + ".join(str(x) for x in getattr(r, "reactants", [])),
            "source_products": " + ".join(str(x) for x in getattr(r, "products", [])),
            "source_Q_MeV": getattr(r, "Q", np.nan),
            "derived_ok": ok,
            "message": msg,
        })

    pd.DataFrame(derived_attempt_rows).to_csv(DATA_DIR / "bbn_network_db_derived_attempts.csv", index=False)

    # Prefer explicit DerivedRate products, but include ReacLib reverse entries as
    # fallback for channels that DerivedRate cannot build.  This maximizes coverage
    # while the diagnostic tells us which rates are truly DB-derived.
    all_rates = unique_rates(forward_rates + derived_rates + reaclib_reverse_rates)

    # Construct PythonNetwork robustly across pynucastro versions.
    try:
        net = pyna.PythonNetwork(rates=all_rates)
    except TypeError:
        lib = pyna.Library(rates=all_rates)
        net = pyna.PythonNetwork(libraries=lib)

    net.write_network(str(BASE / "bbn_network_db.py"))

    rates_df = rate_rows(list(net.rates))
    rates_df.to_csv(DATA_DIR / "bbn_network_db_rates.csv", index=False)

    diag_rows = []
    for channel, reactants, products, note in EXPECTED_REVERSE_CHANNELS:
        rk = tuple(sorted(reactants))
        pk = tuple(sorted(products))
        matched = rates_df[
            (rates_df["reactants_canon"].apply(lambda x: tuple(x.split(" + ")) if x else tuple()) == rk)
            & (rates_df["products_canon"].apply(lambda x: tuple(x.split(" + ")) if x else tuple()) == pk)
        ]
        if len(matched) > 0:
            m = matched.iloc[0]
            diag_rows.append({
                "channel": channel,
                "reactants_expected": " + ".join(reactants),
                "products_expected": " + ".join(products),
                "found": True,
                "matched_rate": m["rate"],
                "Q_MeV": m["Q_MeV"],
                "rate_class": m["rate_class"],
                "derived_from_inverse": m["derived_from_inverse"],
                "note": note,
                "diagnostic": "present in bbn_network_db.py",
            })
        else:
            diag_rows.append({
                "channel": channel,
                "reactants_expected": " + ".join(reactants),
                "products_expected": " + ".join(products),
                "found": False,
                "matched_rate": "",
                "Q_MeV": "",
                "rate_class": "",
                "derived_from_inverse": "",
                "note": note,
                "diagnostic": "missing; bottleneck likely still incomplete",
            })

    diag_df = pd.DataFrame(diag_rows)
    diag_df.to_csv(DATA_DIR / "bbn_network_reverse_diagnostic.csv", index=False)
    return net, list(net.rates), rates_df, diag_df


def make_initial_Y(net, np_ratio: float = 1.0 / 7.0) -> np.ndarray:
    name_to_i = {canon_name(name): i for i, name in enumerate(net.names)}
    Y0 = np.zeros(net.nnuc)
    Yn = np_ratio / (1.0 + np_ratio)
    Yp = 1.0 / (1.0 + np_ratio)
    Y0[name_to_i["n"]] = Yn
    Y0[name_to_i["p"]] = Yp
    return Y0


def nuc_col(net, prefix: str, key: str) -> str:
    for name in net.names:
        if canon_name(name) == key:
            return f"{prefix}_{str(name).replace(' ', '')}"
    raise KeyError(f"Could not find {key} in {net.names}")


def integrate_db_network() -> pd.DataFrame:
    from thermo_model import make_thermo_functions, ETA10_PLANCK_APPROX

    net = import_module(BASE / "bbn_network_db.py", "bbn_network_db")
    snapshots = pd.read_csv(DATA_DIR / "bbn_assignment_snapshots.csv")
    thermo = make_thermo_functions(snapshots, eta10_ref=ETA10_PLANCK_APPROX)

    def rhs_var(t, Y):
        Y_safe = np.maximum(np.asarray(Y, dtype=float), 0.0)
        return np.asarray(net.rhs(t, Y_safe, float(thermo.rho_of_t(t)), float(thermo.T_of_t(t))), dtype=float)

    t0 = 50.0
    t_end = 1000.0
    t_eval = np.unique(np.concatenate([np.geomspace(t0, t_end, 700), np.array([50.0, 200.0, 1000.0])]))
    Y0 = make_initial_Y(net, np_ratio=1.0 / 7.0)

    sol = solve_ivp(
        rhs_var,
        (t0, t_end),
        Y0,
        method="BDF",
        t_eval=t_eval,
        rtol=1.0e-8,
        atol=1.0e-30,
        max_step=2.0,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    Y_hist = np.maximum(sol.y.T, 0.0)
    X_hist = Y_hist * net.A[np.newaxis, :]
    df = pd.DataFrame({
        "model": "pynucastro_db",
        "t0_s": t0,
        "n_over_p_initial": 1.0 / 7.0,
        "p_over_n_initial": 7.0,
        "t_s": sol.t,
        "T9": thermo.T9_of_t(sol.t),
        "T_K": thermo.T_of_t(sol.t),
        "rho_g_cm3": thermo.rho_of_t(sol.t),
        "eta": thermo.eta_of_t(sol.t),
        "eta10": 1.0e10 * thermo.eta_of_t(sol.t),
        "baryon_sum": np.sum(X_hist, axis=1),
    })
    for i, name in enumerate(net.names):
        clean = str(name).replace(" ", "")
        df[f"Y_{clean}"] = Y_hist[:, i]
        df[f"X_{clean}"] = X_hist[:, i]

    Yp = df[nuc_col(net, "Y", "p")].replace(0, np.nan)
    def ratio(key: str):
        try:
            return df[nuc_col(net, "Y", key)] / Yp
        except KeyError:
            return np.nan

    df["D/H"] = ratio("d")
    df["T/H"] = ratio("t")
    df["He3/H"] = ratio("he3")
    df["Li6/H"] = ratio("li6")
    df["Li7/H"] = ratio("li7")
    df["Be7/H"] = ratio("be7")
    df["Li7_plus_Be7_over_H"] = df["Li7/H"].fillna(0.0) + df["Be7/H"].fillna(0.0)
    df["X_He4"] = df[nuc_col(net, "X", "he4")]
    df["X_H1"] = df[nuc_col(net, "X", "p")]
    return df


def snapshot_table(df: pd.DataFrame, snapshot_times=(50.0, 200.0, 1000.0)) -> pd.DataFrame:
    rows = []
    for t in snapshot_times:
        if t < df["t_s"].min() or t > df["t_s"].max():
            continue
        row = {"model": df["model"].iloc[0], "t0_s": df["t0_s"].iloc[0], "t_s": t}
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        for col in numeric_cols:
            row[col] = float(np.interp(t, df["t_s"], df[col]))
        rows.append(row)
    return pd.DataFrame(rows)


def merge_db_outputs(df_db: pd.DataFrame) -> None:
    df_db.to_csv(DATA_DIR / "bbn_evolution_pynucastro_db.csv", index=False)
    final_db = df_db.sort_values("t_s").tail(1).copy()
    final_db.to_csv(DATA_DIR / "bbn_final_abundances_pynucastro_db.csv", index=False)

    evo_path = DATA_DIR / "bbn_evolution_all_models.csv"
    if evo_path.exists():
        evo = pd.read_csv(evo_path)
        evo = evo[evo["model"] != "pynucastro_db"]
        evo = pd.concat([evo, df_db], ignore_index=True, sort=False)
    else:
        evo = df_db.copy()
    evo.to_csv(evo_path, index=False)

    final_path = DATA_DIR / "bbn_final_abundances_all_models.csv"
    if final_path.exists():
        finals = pd.read_csv(final_path)
        finals = finals[finals["model"] != "pynucastro_db"]
        finals = pd.concat([finals, final_db], ignore_index=True, sort=False)
    else:
        finals = final_db.copy()
    finals.to_csv(final_path, index=False)

    snap_path = DATA_DIR / "bbn_snapshots_all_models.csv"
    snap_db = snapshot_table(df_db)
    if snap_path.exists():
        snaps = pd.read_csv(snap_path)
        snaps = snaps[snaps["model"] != "pynucastro_db"]
        snaps = pd.concat([snaps, snap_db], ignore_index=True, sort=False)
    else:
        snaps = snap_db.copy()
    snaps.to_csv(snap_path, index=False)

    obs_cols = ["model", "t0_s", "t_s", "T9", "rho_g_cm3", "eta10", "D/H", "He3/H", "X_He4", "Li7/H", "Be7/H", "Li7_plus_Be7_over_H", "baryon_sum"]
    compact = finals[[c for c in obs_cols if c in finals.columns]].copy()
    compact.to_csv(DATA_DIR / "bbn_final_compact.csv", index=False)

    # Tables read by Typst.
    compact.to_csv(TAB_DIR / "table_final_models.csv", index=False)
    compact.to_latex(TAB_DIR / "table_final_models.tex", index=False, escape=False)
    snaps_cols = [c for c in obs_cols if c in snaps.columns]
    snaps[snaps_cols].to_csv(TAB_DIR / "table_snapshots_models.csv", index=False)
    snaps[snaps_cols].to_latex(TAB_DIR / "table_snapshots_models.tex", index=False, escape=False)


def main() -> int:
    try:
        import pynucastro as pyna  # noqa: F401
    except Exception as exc:
        msg = f"pynucastro is not importable in this environment: {exc}"
        export_not_run_diagnostics(msg)
        write_status("not_run", msg)
        print(msg)
        return 0

    try:
        net, rates, rates_df, diag_df = build_db_network()
    except Exception as exc:
        tb = traceback.format_exc()
        msg = f"Could not build detailed-balance network: {type(exc).__name__}: {exc}"
        export_not_run_diagnostics(msg)
        write_status("build_failed", msg, traceback=tb)
        print(tb)
        return 1

    missing = diag_df[~diag_df["found"].astype(bool)]
    print(f"Detailed-balance candidate network written with {len(rates)} rates.")
    print(f"Important reverse channels found: {len(diag_df) - len(missing)}/{len(diag_df)}")

    try:
        df_db = integrate_db_network()
        merge_db_outputs(df_db)
        final = df_db.sort_values("t_s").tail(1).iloc[0]
        write_status(
            "success",
            "pynucastro_db network built and integrated successfully",
            n_rates=len(rates),
            reverse_channels_found=int(len(diag_df) - len(missing)),
            reverse_channels_total=int(len(diag_df)),
            final_D_over_H=float(final["D/H"]),
            final_He3_over_H=float(final["He3/H"]),
            final_Yp=float(final["X_He4"]),
            final_Li7_over_H=float(final["Li7_plus_Be7_over_H"]),
        )
        print("Integration OK. Final diagnostic:")
        print(pd.read_csv(DATA_DIR / "pynucastro_db_status.csv").to_string(index=False))
    except Exception as exc:
        tb = traceback.format_exc()
        msg = f"Network was built, but integration failed: {type(exc).__name__}: {exc}"
        write_status(
            "integration_failed",
            msg,
            n_rates=len(rates),
            reverse_channels_found=int(len(diag_df) - len(missing)),
            reverse_channels_total=int(len(diag_df)),
            traceback=tb,
        )
        print(tb)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
