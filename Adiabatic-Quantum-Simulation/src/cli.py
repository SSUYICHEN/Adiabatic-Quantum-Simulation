"""
cli.py - Command-line interface for the SSH-Hubbard QA toolkit.

Usage
-----
  python -m sshh_qa fidelity      [options]   # step 0: choose T_A and steps
  python -m sshh_qa berry         [options]   # measure the Berry phase (PBC)
  python -m sshh_qa polarization  [options]   # measure the polarization (OBC)
  python -m sshh_qa plot --kind ... --data ...  # (re)draw from saved JSON

Run `python -m sshh_qa <command> -h` for the options of each command.
"""

from __future__ import annotations

import argparse
import os

from . import __version__
from . import experiments as ex
from . import plotting as pl


# --------------------------- arg helpers -----------------------------
def _floats(s):
    return [float(x) for x in str(s).replace(" ", "").split(",") if x != ""]


def _ints(s):
    return [int(x) for x in str(s).replace(" ", "").split(",") if x != ""]


DEFAULT_OUT = "QA_results"


# ------------------------------ fidelity -----------------------------
def cmd_fidelity(a):
    out_json = a.out or os.path.join(DEFAULT_OUT, "fidelity",
                                     f"fidelity_N{a.N}_U{a.U:g}_{a.reference}.json")
    print(f"[fidelity] N={a.N} t1={a.t1} t2={a.t2} Ne={a.electrons} U={a.U} "
          f"{'PBC' if a.pbc else 'OBC'} reference={a.reference}")
    result = ex.fidelity_scan(
        N_cells=a.N, t1=a.t1, t2=a.t2, number_of_electrons=a.electrons,
        U=a.U, T_A_values=a.TA, steps_list=a.steps, PBC=a.pbc,
        reference=a.reference, out_json=out_json,
    )
    print(f"  -> saved JSON: {out_json}")
    _report_fidelity(result)
    if not a.no_plot:
        img = os.path.splitext(out_json)[0] + ".png"
        pl.plot_fidelity(result, img)
        print(f"  -> saved figure: {img}")


def _report_fidelity(result):
    """Print a compact table and a suggested (T_A, steps) operating point."""
    steps = result["params"]["steps_list"]
    print("  steps : " + "  ".join(f"{s:>5d}" for s in steps))
    best = None
    for T_A, fids in sorted(result["curves"].items()):
        print(f"  T={T_A:<4g}: " + "  ".join(f"{f:5.3f}" for f in fids))
        for s, f in zip(steps, fids):
            if f >= 0.999 and (best is None or s < best[1]):
                best = (T_A, s, f)
    if best:
        print(f"  suggestion: fidelity>=0.999 first reached at "
              f"T_A={best[0]:g}, steps={best[1]} (fidelity={best[2]:.4f})")


# ------------------------------- berry -------------------------------
def cmd_berry(a):
    out_dir = a.out or os.path.join(DEFAULT_OUT, "berry",
                                    f"N{a.N}_TA{a.TA:g}_steps{a.steps}_UA{a.UA:g}")
    print(f"[berry] N={a.N} Ne={a.electrons} t1={a.t1} UA={a.UA} "
          f"deltaU={a.delta_U} TA={a.TA} steps={a.steps}")
    records = ex.berry_sweep(
        N_cells=a.N, number_of_electrons=a.electrons, t1=a.t1,
        t2_values=a.t2, U_A=a.UA, delta_U_values=a.delta_U,
        T_A=a.TA, steps=a.steps, mode=a.mode, out_dir=out_dir,
    )
    print(f"  -> {len(records)} points saved to: {out_dir}")
    for d in records:
        print(f"    dU={d['delta_U']:<7g} w={d['t2']:<5g} "
              f"|z|={d['Twist_Amplitude']:.4f} "
              f"gamma={d['Berry_Phase_pi_wrapped']:+.4f} pi")
    if not a.no_plot:
        img = os.path.join(out_dir, "BerryPhase_PhaseDiagram.png")
        pl.plot_berry_phase(records, img, t1_default=a.t1)
        print(f"  -> saved figure: {img}")


# --------------------------- polarization ----------------------------
def cmd_polarization(a):
    electrons = a.electrons if a.electrons is not None else 2 * a.N + 2
    out_dir = a.out or os.path.join(DEFAULT_OUT, "polarization",
                                    f"N{a.N}_TA{a.TA:g}_steps{a.steps}_t1{a.t1:g}")
    print(f"[polarization] N={a.N} Ne={electrons} t1={a.t1} t2={a.t2} "
          f"UA={a.UA} deltaU={a.delta_U} TA={a.TA} steps={a.steps} (OBC)")
    records = ex.polarization_sweep(
        N_cells=a.N, total_electrons=electrons, t1=a.t1, t2=a.t2,
        U_A=a.UA, delta_U_values=a.delta_U, T_A=a.TA, steps=a.steps,
        out_dir=out_dir,
    )
    print(f"  -> {len(records)} curves saved to: {out_dir}")
    for r in records:
        du = r["metadata"]["delta_U"]
        edge = r["data"]["n_A_minus_n_B"][0]
        print(f"    dU={du:<7g} edge polarization (cell 0) = {edge:+.4f}")
    if not a.no_plot:
        img = os.path.join(out_dir, "Polarization_vs_DeltaU.png")
        pl.plot_polarization(records, img)
        print(f"  -> saved figure: {img}")


# ------------------------------- plot --------------------------------
def cmd_plot(a):
    out = a.out or (os.path.splitext(a.data)[0] + ".png"
                    if a.data.endswith(".json") else
                    os.path.join(a.data, f"{a.kind}.png"))
    if a.kind == "fidelity":
        import json
        with open(a.data, "r", encoding="utf-8") as f:
            pl.plot_fidelity(json.load(f), out)
    elif a.kind == "berry":
        pl.plot_berry_phase(a.data, out, t1_default=a.t1)
    elif a.kind == "polarization":
        pl.plot_polarization(a.data, out)
    print(f"  -> saved figure: {out}")


# ------------------------------ parser -------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        prog="sshh_qa",
        description="SSH-Hubbard adiabatic quantum simulation toolkit.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--version", action="version",
                   version=f"sshh_qa {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    # fidelity
    f = sub.add_parser("fidelity", help="Trotter-convergence check (choose T_A & steps)",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    f.add_argument("--N", type=int, default=6, help="number of unit cells")
    f.add_argument("--t1", type=float, default=0.5)
    f.add_argument("--t2", type=float, default=1.5)
    f.add_argument("--electrons", type=int, default=12, help="total electrons")
    f.add_argument("--U", type=float, default=0.0,
                   help="uniform Hubbard U (use 0 for the non-interacting check)")
    f.add_argument("--TA", type=_floats, default=[1, 5, 15, 80],
                   help="comma-separated annealing times")
    f.add_argument("--steps", type=_ints,
                   default=[1, 10, 20, 30, 40, 50, 60, 80, 100, 120, 150],
                   help="comma-separated Trotter step counts")
    f.add_argument("--reference", choices=["ground", "adjacent"], default="ground",
                   help="'ground' vs psi(0) (U=0 check); 'adjacent' between steps (U!=0)")
    bnd = f.add_mutually_exclusive_group()
    bnd.add_argument("--pbc", dest="pbc", action="store_true", default=True)
    bnd.add_argument("--obc", dest="pbc", action="store_false")
    f.add_argument("--out", help="output JSON path")
    f.add_argument("--no-plot", action="store_true")
    f.set_defaults(func=cmd_fidelity)

    # berry
    b = sub.add_parser("berry", help="Berry phase vs inter-cell hopping (PBC)",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    b.add_argument("--N", type=int, default=6)
    b.add_argument("--electrons", type=int, default=12)
    b.add_argument("--t1", type=float, default=1.0)
    b.add_argument("--t2", type=_floats,
                   default=[0, 0.25, 0.5, 0.75, 0.99, 1.01, 1.25, 1.5, 1.75, 2.0],
                   help="comma-separated inter-cell hoppings w to scan")
    b.add_argument("--UA", type=float, default=0.01, help="U on the A sublattice")
    b.add_argument("--delta-U", dest="delta_U", type=_floats, default=[0.0],
                   help="comma-separated Delta U = U_B - U_A values")
    b.add_argument("--TA", type=float, default=1.0)
    b.add_argument("--steps", type=int, default=40)
    b.add_argument("--mode", choices=["up_spin", "total"], default="up_spin")
    b.add_argument("--out", help="output directory")
    b.add_argument("--no-plot", action="store_true")
    b.set_defaults(func=cmd_berry)

    # polarization
    pol = sub.add_parser("polarization", help="electron polarization vs Delta U (OBC)",
                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pol.add_argument("--N", type=int, default=6)
    pol.add_argument("--electrons", type=int, default=None,
                     help="total electrons (default = half-filling+2 = 2N+2)")
    pol.add_argument("--t1", type=float, default=0.5)
    pol.add_argument("--t2", type=float, default=1.5)
    pol.add_argument("--UA", type=float, default=1.0)
    pol.add_argument("--delta-U", dest="delta_U", type=_floats,
                     default=[0, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3],
                     help="comma-separated Delta U values")
    pol.add_argument("--TA", type=float, default=1.0)
    pol.add_argument("--steps", type=int, default=40)
    pol.add_argument("--out", help="output directory")
    pol.add_argument("--no-plot", action="store_true")
    pol.set_defaults(func=cmd_polarization)

    # plot
    pt = sub.add_parser("plot", help="(re)draw a figure from saved JSON",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pt.add_argument("--kind", required=True,
                    choices=["fidelity", "berry", "polarization"])
    pt.add_argument("--data", required=True,
                    help="JSON file (fidelity) or directory (berry/polarization)")
    pt.add_argument("--t1", type=float, default=1.0,
                    help="t1 used for the berry gap threshold (fallback)")
    pt.add_argument("--out", help="output PNG path")
    pt.set_defaults(func=cmd_plot)

    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
