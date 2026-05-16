#!/usr/bin/env python3
"""Build and test an alternative pynucastro network with reverse rates.

Diagnostic script for the deuterium-bottleneck problem. It tries to build a
second network, ``bbn_network_db.py``, with inverse rates by detailed balance.
It always writes diagnostic CSVs and exits with code 0 so that the notebook can
show the actual status instead of stopping with a bare CalledProcessError.
"""
from __future__ import annotations

from pathlib import Path
import importlib.util
import traceback
from typing import Iterable, Any

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


def canon_name(obj: Any) -> str:
    s = str(obj).strip().lower().replace(" ", "")
    return ALIAS_TO_CANON.get(s, s)


def species_key(items: Iterable[Any]) -> tuple[str, ...]:
    return tuple(sorted(canon_name(x) for x in items))


def library_rates(lib: Any) -> list[Any]:
    """Return rates from a pynucastro Library across API versions."""
    if lib is None:
        return []
    if hasattr(lib, "get_rates"):
        return list(lib.get_rates())
    if hasattr(lib, "rates"):
        return list(lib.rates)
    return []


def network_rates(net: Any) -> list[Any]:
    if hasattr(net, "rates"):
        return list(net.rates)
    if hasattr(net, "get_rates"):
        return list(net.get_rates())
    return []


def rate_id_key(r: Any) -> tuple[Any, ...]:
    return (
        species_key(getattr(r, "reactants", [])),
        species_key(getattr(r, "products", [])),
        getattr(r, "id", None),
        getattr(r, "fname", ""),
        getattr(r, "label", ""),
        type(r).__name__,
        bool(getattr(r, "weak", False)),
        bool(getattr(r, "derived_from_inverse", False)),
    )


def structural_rate_key(r: Any) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Key used by pynucastro to decide whether two rates duplicate a channel.

    Different ReacLib entries, explicit DerivedRate calls and Library.derived_backward()
    can produce several rates with the same reactants and products. PythonNetwork
    refuses those duplicates. For this diagnostic network we keep one representative
    channel and record the discarded candidates.
    """
    return (species_key(getattr(r, "reactants", [])), species_key(getattr(r, "products", [])))


def rate_preference_score(r: Any) -> tuple[int, int, int, str]:
    """Lower is better when choosing between duplicate channels."""
    # For reverse photodisintegration channels prefer rates derived by detailed
    # balance, because that is exactly what this diagnostic is testing.
    derived = bool(getattr(r, "derived_from_inverse", False)) or type(r).__name__.lower().startswith("derived")
    weak = bool(getattr(r, "weak", False))
    text = str(r)
    return (0 if derived else 1, 0 if not weak else 1, len(text), text)


def unique_rates(rates: Iterable[Any]) -> list[Any]:
    """Exact de-duplication, preserving genuinely distinct channels."""
    out, seen = [], set()
    for r in rates:
        key = rate_id_key(r)
        if key not in seen:
            out.append(r)
            seen.add(key)
    return out


def deduplicate_channels_for_pythonnetwork(rates: Iterable[Any]) -> tuple[list[Any], pd.DataFrame]:
    """Remove channel duplicates that trigger RateDuplicationError.

    Returns the selected rates and a CSV-friendly table of discarded duplicates.
    """
    groups: dict[tuple[tuple[str, ...], tuple[str, ...]], list[Any]] = {}
    for r in unique_rates(rates):
        groups.setdefault(structural_rate_key(r), []).append(r)

    selected: list[Any] = []
    rows = []
    for key, group in groups.items():
        ordered = sorted(group, key=rate_preference_score)
        keep = ordered[0]
        selected.append(keep)
        if len(ordered) > 1:
            for cand in ordered:
                rows.append({
                    "reactants_canon": " + ".join(key[0]),
                    "products_canon": " + ".join(key[1]),
                    "kept": cand is keep,
                    "rate": str(cand),
                    "rate_class": type(cand).__name__,
                    "fname": getattr(cand, "fname", ""),
                    "label": getattr(cand, "label", ""),
                    "Q_MeV": getattr(cand, "Q", np.nan),
                    "weak": bool(getattr(cand, "weak", False)),
                    "derived_from_inverse": bool(getattr(cand, "derived_from_inverse", False)),
                })
    return selected, pd.DataFrame(rows)


def rate_rows(rates: list[Any]) -> pd.DataFrame:
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


def write_status(status: str, message: str, **extra: Any) -> None:
    row = {"status": status, "message": message}
    row.update(extra)
    pd.DataFrame([row]).to_csv(DATA_DIR / "pynucastro_db_status.csv", index=False)


def export_empty_diagnostics(reason: str) -> None:
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
    pd.DataFrame(columns=["source_rate", "derived_ok", "message"]).to_csv(
        DATA_DIR / "bbn_network_db_derived_attempts.csv", index=False
    )
    pd.DataFrame(columns=["index", "rate", "reactants", "products"]).to_csv(
        DATA_DIR / "bbn_network_db_rates.csv", index=False
    )


def import_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def get_library_class(pyna):
    if hasattr(pyna, "Library"):
        return pyna.Library
    from pynucastro.rates import Library  # type: ignore
    return Library


def get_derived_rate_class(pyna):
    if hasattr(pyna, "DerivedRate"):
        return pyna.DerivedRate
    try:
        from pynucastro.rates import DerivedRate  # type: ignore
    except Exception:
        from pynucastro.rates.derived_rate import DerivedRate  # type: ignore
    return DerivedRate


def make_network_from_rates(pyna, rates: list[Any]):
    Library = get_library_class(pyna)
    errors = []
    try:
        return pyna.PythonNetwork(rates=rates)
    except Exception as exc:
        errors.append(f"PythonNetwork(rates=...): {type(exc).__name__}: {exc}")

    try:
        lib = Library(rates=rates)
        return pyna.PythonNetwork(libraries=lib)
    except Exception as exc:
        errors.append(f"PythonNetwork(libraries=Library(...)): {type(exc).__name__}: {exc}")

    try:
        lib = Library(rates=rates)
        return pyna.PythonNetwork(libraries=[lib])
    except Exception as exc:
        errors.append(f"PythonNetwork(libraries=[Library(...)]): {type(exc).__name__}: {exc}")

    raise RuntimeError("Could not construct PythonNetwork. Tried:\n" + "\n".join(errors))


def build_db_network() -> tuple[Any, list[Any], pd.DataFrame, pd.DataFrame]:
    import pynucastro as pyna

    rl = pyna.ReacLibLibrary()
    lib_forward = rl.linking_nuclei(NUCLEI_LIGHT, with_reverse=False)
    forward_rates = library_rates(lib_forward)

    if not forward_rates:
        fallback_nuclei = ["n", "p", "d", "t", "he3", "he4", "li7", "be7"]
        lib_forward = rl.linking_nuclei(fallback_nuclei, with_reverse=False)
        forward_rates = library_rates(lib_forward)

    if not forward_rates:
        raise RuntimeError("No forward rates were selected from ReacLib. Check the pynucastro/ReacLib installation.")

    derived_rates = []
    derived_attempt_rows = []

    if hasattr(lib_forward, "derived_backward"):
        try:
            lib_derived = lib_forward.derived_backward(use_pf=False)
            derived_rates = library_rates(lib_derived)
            for dr in derived_rates:
                source = getattr(dr, "source_rate", None)
                derived_attempt_rows.append({
                    "source_rate": getattr(source, "id", str(source)),
                    "derived_rate": str(dr),
                    "derived_ok": True,
                    "message": "created by Library.derived_backward(use_pf=False)",
                })
        except Exception as exc:
            derived_attempt_rows.append({
                "source_rate": "Library.derived_backward",
                "derived_rate": "",
                "derived_ok": False,
                "message": f"{type(exc).__name__}: {exc}",
            })

    DerivedRate = get_derived_rate_class(pyna)
    existing_derived_keys = {rate_id_key(r) for r in derived_rates}
    allowed = set(NUCLEI_LIGHT)
    for r in forward_rates:
        try:
            dr = DerivedRate(source_rate=r, use_pf=False)
            all_nucs = [canon_name(x) for x in list(getattr(dr, "reactants", [])) + list(getattr(dr, "products", []))]
            if all(x in allowed for x in all_nucs):
                if rate_id_key(dr) not in existing_derived_keys:
                    derived_rates.append(dr)
                    existing_derived_keys.add(rate_id_key(dr))
                ok = True
                msg = "created by explicit DerivedRate(source_rate=..., use_pf=False)"
                dr_text = str(dr)
            else:
                ok = False
                msg = "derived rate leaves selected light-nuclide set"
                dr_text = str(dr)
        except Exception as exc:
            ok = False
            msg = f"{type(exc).__name__}: {exc}"
            dr_text = ""
        derived_attempt_rows.append({
            "source_rate": str(r),
            "source_reactants": " + ".join(str(x) for x in getattr(r, "reactants", [])),
            "source_products": " + ".join(str(x) for x in getattr(r, "products", [])),
            "source_Q_MeV": getattr(r, "Q", np.nan),
            "derived_rate": dr_text,
            "derived_ok": ok,
            "message": msg,
        })

    pd.DataFrame(derived_attempt_rows).to_csv(DATA_DIR / "bbn_network_db_derived_attempts.csv", index=False)

    # Candidate pool. We use direct forward rates plus the detailed-balance
    # backward rates. ReacLib reverse rates are only a fallback because they can
    # duplicate the DerivedRate channels and trigger RateDuplicationError.
    lib_with_reverse = rl.linking_nuclei(NUCLEI_LIGHT, with_reverse=True)
    reaclib_with_reverse_rates = library_rates(lib_with_reverse)

    candidate_rates = unique_rates(forward_rates + derived_rates)
    if not candidate_rates:
        raise RuntimeError("Candidate detailed-balance network contains zero rates.")

    # First remove exact/structural duplicates. This is the important fix for
    # the previous build failure.
    all_rates, duplicate_df = deduplicate_channels_for_pythonnetwork(candidate_rates)
    duplicate_df.to_csv(DATA_DIR / "bbn_network_db_duplicate_candidates.csv", index=False)

    try:
        net = make_network_from_rates(pyna, all_rates)
    except Exception as first_exc:
        # Fallback: include ReacLib reverse rates too, then deduplicate again.
        # This helps with older pynucastro versions where DerivedRate creation
        # silently fails for some channels.
        fallback_rates, duplicate_df2 = deduplicate_channels_for_pythonnetwork(
            forward_rates + derived_rates + reaclib_with_reverse_rates
        )
        duplicate_df2.to_csv(DATA_DIR / "bbn_network_db_duplicate_candidates.csv", index=False)
        try:
            net = make_network_from_rates(pyna, fallback_rates)
            all_rates = fallback_rates
        except Exception as second_exc:
            raise RuntimeError(
                "Could not construct PythonNetwork after structural de-duplication. "
                f"First attempt: {type(first_exc).__name__}: {first_exc}. "
                f"Fallback attempt: {type(second_exc).__name__}: {second_exc}"
            )
    net.write_network(str(BASE / "bbn_network_db.py"))

    rates = network_rates(net)
    rates_df = rate_rows(rates)
    rates_df.to_csv(DATA_DIR / "bbn_network_db_rates.csv", index=False)

    def split_key(s: Any) -> tuple[str, ...]:
        if pd.isna(s) or str(s).strip() == "":
            return tuple()
        return tuple(str(s).split(" + "))

    diag_rows = []
    for channel, reactants, products, note in EXPECTED_REVERSE_CHANNELS:
        rk = tuple(sorted(reactants))
        pk = tuple(sorted(products))
        matched = rates_df[
            (rates_df["reactants_canon"].apply(split_key) == rk)
            & (rates_df["products_canon"].apply(split_key) == pk)
        ]
        if len(matched) > 0:
            matched2 = matched[matched["derived_from_inverse"].astype(bool)]
            m = (matched2 if len(matched2) else matched).iloc[0]
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
    return net, rates, rates_df, diag_df


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
    snapshots_path = DATA_DIR / "bbn_assignment_snapshots.csv"
    if not snapshots_path.exists():
        raise FileNotFoundError(f"Missing {snapshots_path}; run notebook 00 first.")
    snapshots = pd.read_csv(snapshots_path)
    thermo = make_thermo_functions(snapshots, eta10_ref=ETA10_PLANCK_APPROX)

    def rhs_var(t, Y):
        Y_safe = np.maximum(np.asarray(Y, dtype=float), 0.0)
        return np.asarray(net.rhs(t, Y_safe, float(thermo.rho_of_t(t)), float(thermo.T_of_t(t))), dtype=float)

    t0 = 50.0
    t_end = 1000.0
    t_eval = np.unique(np.concatenate([np.geomspace(t0, t_end, 500), np.array([50.0, 200.0, 1000.0])]))
    Y0 = make_initial_Y(net, np_ratio=1.0 / 7.0)

    sol = solve_ivp(
        rhs_var,
        (t0, t_end),
        Y0,
        method="BDF",
        t_eval=t_eval,
        rtol=1.0e-8,
        atol=1.0e-30,
        max_step=5.0,
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
    compact.to_csv(TAB_DIR / "table_final_models.csv", index=False)
    compact.to_latex(TAB_DIR / "table_final_models.tex", index=False, escape=False)

    snap_cols = [c for c in obs_cols if c in snaps.columns]
    snaps[snap_cols].to_csv(TAB_DIR / "table_snapshots_models.csv", index=False)
    snaps[snap_cols].to_latex(TAB_DIR / "table_snapshots_models.tex", index=False, escape=False)


def main() -> int:
    try:
        import pynucastro as pyna  # noqa: F401
    except Exception as exc:
        msg = f"pynucastro is not importable in this environment: {exc}"
        export_empty_diagnostics(msg)
        write_status("not_run", msg)
        print(msg)
        return 0

    try:
        _net_obj, rates, _rates_df, diag_df = build_db_network()
    except Exception as exc:
        tb = traceback.format_exc()
        msg = f"Could not build detailed-balance network: {type(exc).__name__}: {exc}"
        export_empty_diagnostics(msg)
        write_status("build_failed", msg, traceback=tb)
        print(tb)
        return 0

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
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
