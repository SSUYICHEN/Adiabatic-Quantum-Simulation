"""
aqs_braket.py - Command line front-end for running the SSHH experiments on
AWS Braket (IonQ trapped-ion QPUs).

IonQ on Braket costs REAL money. Recommended order:

  python aqs_braket.py targets
  python aqs_braket.py polarization --estimate
  python aqs_braket.py polarization --target SV1 --shots 200
  python aqs_braket.py polarization --target "Forte 1" --shots 200 --i-accept-cost
  python aqs_braket.py fetch --manifest "braket_results/polarization/manifest_*.json" --plot
"""

import argparse
import glob
import json

from aqs import braket as bk


def _floats(s):
    return [float(x) for x in str(s).replace(" ", "").split(",") if x != ""]


def _resolve(manifest_glob):
    hits = sorted(glob.glob(manifest_glob))
    if not hits:
        raise SystemExit(f"no manifest matches: {manifest_glob}")
    return hits[-1]


def cmd_targets(a):
    bk.list_targets()


def cmd_polarization(a):
    bk.submit_polarization(
        N_cells=a.N, electrons=a.electrons, v=a.v, w=a.w, U_A=a.UA,
        delta_U_values=a.delta_U, T_A=a.TA, steps=a.steps, shots=a.shots,
        binary=a.binary, target=a.target, out_dir=a.out,
        estimate_only=a.estimate, accept_cost=a.i_accept_cost)


def cmd_berry(a):
    bk.submit_berry(
        N_cells=a.N, electrons=a.electrons, v=a.v, w_values=a.w, U_A=a.UA,
        delta_U_values=a.delta_U, T_A=a.TA, steps=a.steps, shots=a.shots,
        binary=a.binary, target=a.target, out_dir=a.out,
        estimate_only=a.estimate, accept_cost=a.i_accept_cost)


def cmd_fetch(a):
    manifest = _resolve(a.manifest)
    records, out_dir = bk.fetch_results(manifest, out_dir=a.out,
                                        save_counts=a.save_counts)
    if a.plot and records:
        from aqs import plotting as pl
        with open(manifest, "r", encoding="utf-8") as f:
            exp = json.load(f)["experiment"]
        if exp == "berry":
            pl.plot_berry_phase(out_dir, f"{out_dir}/BerryPhase_Braket.png")
            print(f"  figure -> {out_dir}/BerryPhase_Braket.png")
        else:
            pl.plot_polarization(out_dir, f"{out_dir}/Polarization_Braket.png")
            try:
                pl.plot_polarization_cross_section(
                    out_dir, f"{out_dir}/Polarization_CrossSection_Braket.png", cell_index=0)
            except Exception as e:
                print(f"  (cross section skipped: {e})")
            print(f"  figure -> {out_dir}/Polarization_Braket.png")


def build_parser():
    p = argparse.ArgumentParser(
        prog="aqs_braket",
        description="Run SSHH experiments on AWS Braket (IonQ).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    sub.add_parser("targets", help="list Braket devices visible to this account"
                   ).set_defaults(func=cmd_targets)

    def common(sp, default_out):
        sp.add_argument("--N", type=int, default=6)
        sp.add_argument("--UA", type=float, default=0.01)
        sp.add_argument("--delta-U", dest="delta_U", type=_floats,
                        default=[0.0, 0.01, 0.1, 1.0])
        sp.add_argument("--TA", type=float, default=1.0)
        sp.add_argument("--steps", type=int, default=20, help="Trotter intervals L")
        sp.add_argument("--shots", type=int, default=200,
                        help="IonQ is expensive; default kept low")
        sp.add_argument("--binary", type=int, default=1, choices=[0, 1],
                        help="1=full evolution (default), 0=state preparation only")
        sp.add_argument("--target", default=bk.SIMULATOR,
                        help='"SV1" (default) / "Forte 1" / "Forte Enterprise 1"')
        sp.add_argument("--estimate", action="store_true",
                        help="print cost estimate only; submit nothing")
        sp.add_argument("--i-accept-cost", dest="i_accept_cost",
                        action="store_true",
                        help="required to submit to a QPU (acknowledges real cost)")
        sp.add_argument("--out", default=default_out)

    b = sub.add_parser("berry", formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                       help="Berry phase vs w (PBC)")
    b.add_argument("--electrons", type=int, default=12)
    b.add_argument("--v", type=float, default=1.0)
    b.add_argument("--w", type=_floats, default=[0, 0.5, 0.99, 1.01, 1.5, 2.0])
    common(b, "braket_results/berry")
    b.set_defaults(func=cmd_berry)

    pol = sub.add_parser("polarization", formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                         help="sublattice polarization vs dU (OBC)")
    pol.add_argument("--electrons", type=int, default=14)
    pol.add_argument("--v", type=float, default=0.5)
    pol.add_argument("--w", type=float, default=1.5)
    common(pol, "braket_results/polarization")
    pol.set_defaults(func=cmd_polarization)

    ft = sub.add_parser("fetch", help="retrieve results and post-process")
    ft.add_argument("--manifest", required=True, help="path or glob")
    ft.add_argument("--out", default=None)
    ft.add_argument("--plot", action="store_true")
    ft.add_argument("--save-counts", action="store_true")
    ft.set_defaults(func=cmd_fetch)

    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    args.func(args)