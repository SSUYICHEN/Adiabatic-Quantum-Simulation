"""
aqs_azure.py - Command line front-end for running the SSHH experiments on
Azure Quantum (Quantinuum H2).

Recommended order (Quantinuum shots are expensive -- validate before you spend):

  # 0. list the machines this workspace can see
  python aqs_azure.py targets

  # 1. FREE syntax check on H2-1sc (validates the circuit, charges nothing)
  python aqs_azure.py polarization --syntax-check

  # 2. estimate the HQC cost on the emulator (no charge)
  python aqs_azure.py polarization --estimate

  # 3. submit to the emulator (this DOES cost HQC)
  python aqs_azure.py polarization --shots 500

  # 4. retrieve + post-process + plot
  python aqs_azure.py fetch --manifest azure_results/polarization/manifest_*.json --plot

Defaults match the IBM runs so the three data sets are directly comparable:
  berry        : UA=0.01, dU in {0, 0.001, 0.01,0.1,1}, v=1,   w=0..2, 12 e-, PBC, L=40
  polarization : UA=0.01, dU in {0, 0.001, 0.01,0.1,1}, v=0.5, w=1.5,  14 e-, OBC, L=40
"""

import argparse
import glob
import json

from aqs import azure as az


def _floats(s):
    return [float(x) for x in str(s).replace(" ", "").split(",") if x != ""]


def _resolve(manifest_glob):
    hits = sorted(glob.glob(manifest_glob))
    if not hits:
        raise SystemExit(f"no manifest matches: {manifest_glob}")
    return hits[-1]  # newest


def cmd_targets(a):
    az.list_targets()


def cmd_polarization(a):
    az.submit_polarization(
        N_cells=a.N, electrons=a.electrons, v=a.v, w=a.w, U_A=a.UA,
        delta_U_values=a.delta_U, T_A=a.TA, steps=a.steps, shots=a.shots,
        target=a.target, out_dir=a.out,
        syntax_check=a.syntax_check, estimate_only=a.estimate)


def cmd_berry(a):
    az.submit_berry(
        N_cells=a.N, electrons=a.electrons, v=a.v, w_values=a.w, U_A=a.UA,
        delta_U_values=a.delta_U, T_A=a.TA, steps=a.steps, shots=a.shots,
        target=a.target, out_dir=a.out,
        syntax_check=a.syntax_check, estimate_only=a.estimate)


def cmd_fetch(a):
    manifest = _resolve(a.manifest)
    records, out_dir = az.fetch_results(manifest, out_dir=a.out,
                                        save_counts=a.save_counts)
    if a.plot and records:
        from aqs import plotting as pl
        with open(manifest, "r", encoding="utf-8") as f:
            exp = json.load(f)["experiment"]
        if exp == "berry":
            pl.plot_berry_phase(out_dir, f"{out_dir}/BerryPhase_Azure.png")
            print(f"  figure -> {out_dir}/BerryPhase_Azure.png")
        else:
            pl.plot_polarization(out_dir, f"{out_dir}/Polarization_Azure.png")
            try:
                pl.plot_polarization_cross_section(
                    out_dir, f"{out_dir}/Polarization_CrossSection_Azure.png", cell_index=0)
            except Exception as e:
                print(f"  (cross section skipped: {e})")
            print(f"  figure -> {out_dir}/Polarization_Azure.png")


def build_parser():
    p = argparse.ArgumentParser(
        prog="aqs_azure",
        description="Run SSHH experiments on Azure Quantum (Quantinuum H2).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    sub.add_parser("targets", help="list backends visible to this workspace"
                   ).set_defaults(func=cmd_targets)

    def common(sp, default_out):
        sp.add_argument("--N", type=int, default=6)
        sp.add_argument("--UA", type=float, default=0.01)
        sp.add_argument("--delta-U", dest="delta_U", type=_floats,
                        default=[0.0, 0.01, 0.1, 1.0])
        sp.add_argument("--TA", type=float, default=1.0)
        sp.add_argument("--steps", type=int, default=40, help="Trotter intervals L")
        sp.add_argument("--shots", type=int, default=500,
                        help="Quantinuum shots are costly; default kept low")
        sp.add_argument("--target", default=az.EMULATOR,
                        help=f"{az.EMULATOR} (default) / {az.SYNTAX_CHECKER} / {az.HARDWARE}")
        sp.add_argument("--syntax-check", action="store_true",
                        help=f"force {az.SYNTAX_CHECKER} (FREE validation)")
        sp.add_argument("--estimate", action="store_true",
                        help="estimate HQC cost only; submit nothing")
        sp.add_argument("--out", default=default_out)

    b = sub.add_parser("berry", formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                       help="Berry phase vs w (PBC)")
    b.add_argument("--electrons", type=int, default=12)
    b.add_argument("--v", type=float, default=1.0)
    b.add_argument("--w", type=_floats, default=[0, 0.5, 0.99, 1.01, 1.5, 2.0])
    common(b, "azure_results/berry")
    b.set_defaults(func=cmd_berry)

    pol = sub.add_parser("polarization", formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                         help="sublattice polarization vs dU (OBC)")
    pol.add_argument("--electrons", type=int, default=14)
    pol.add_argument("--v", type=float, default=0.5)
    pol.add_argument("--w", type=float, default=1.5)
    common(pol, "azure_results/polarization")
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