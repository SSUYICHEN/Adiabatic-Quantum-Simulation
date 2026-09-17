"""
aqs_ibm.py - Command line front-end for running the SSHH experiments on
IBM Quantum hardware.

Typical session
---------------
  # 0. See circuit sizes and predicted fidelity WITHOUT using any quota:
  python aqs_ibm.py berry        --dry-run
  python aqs_ibm.py polarization --dry-run

  # 1. Submit (records job IDs to a manifest JSON):
  python aqs_ibm.py berry        --shots 8192
  python aqs_ibm.py polarization --shots 8192

  # 2. Check on them later (queues can be hours):
  python aqs_ibm.py status --manifest ibm_results/berry/manifest_berry.json

  # 3. Retrieve + post-process + plot:
  python aqs_ibm.py fetch --manifest ibm_results/berry/manifest_berry.json --plot

Defaults follow the settings agreed for this study:
  berry        : UA=0.01, dU in {0,0.01,0.1,1}, v=1,   w=0..2, 12 e-, PBC
  polarization : UA=0.01, dU in {0,0.01,0.1,1}, v=0.5, w=1.5,  14 e-, OBC
  both         : T_A=1, L=20
"""

import argparse
import os

from aqs import ibm


def _floats(s):
    return [float(x) for x in str(s).replace(" ", "").split(",") if x != ""]


def cmd_berry(a):
    ibm.submit_berry(
        N_cells=a.N, electrons=a.electrons, v=a.v, w_values=a.w,
        U_A=a.UA, delta_U_values=a.delta_U, T_A=a.TA, steps=a.steps,
        shots=a.shots, backend_name=a.backend,
        optimization_level=a.opt, out_dir=a.out, dry_run=a.dry_run)


def cmd_polarization(a):
    ibm.submit_polarization(
        N_cells=a.N, electrons=a.electrons, v=a.v, w=a.w,
        U_A=a.UA, delta_U_values=a.delta_U, T_A=a.TA, steps=a.steps,
        shots=a.shots, backend_name=a.backend,
        optimization_level=a.opt, out_dir=a.out, dry_run=a.dry_run)


def cmd_status(a):
    ibm.job_status(a.manifest)


def cmd_fetch(a):
    records, out_dir = ibm.fetch_results(a.manifest, out_dir=a.out,
                                         save_counts=a.save_counts)
    if a.plot and records:
        from aqs import plotting as pl
        import json
        with open(a.manifest, "r", encoding="utf-8") as f:
            exp = json.load(f)["experiment"]
        if exp == "berry":
            img = os.path.join(out_dir, "BerryPhase_PhaseDiagram_IBM.png")
            pl.plot_berry_phase(out_dir, img)
            print(f"  figure -> {img}")
            try:
                w_cut = max(r["w"] for r in records)
                img2 = os.path.join(out_dir, f"BerryPhase_CrossSection_w_{w_cut:g}_IBM.png")
                pl.plot_berry_cross_section(out_dir, img2, w_cut=w_cut)
                print(f"  figure -> {img2}")
            except Exception as e:
                print(f"  (cross section skipped: {e})")
        else:
            img = os.path.join(out_dir, "Polarization_vs_DeltaU_IBM.png")
            pl.plot_polarization(out_dir, img)
            print(f"  figure -> {img}")
            try:
                img2 = os.path.join(out_dir, "Polarization_CrossSection_cell_0_IBM.png")
                pl.plot_polarization_cross_section(out_dir, img2, cell_index=0)
                print(f"  figure -> {img2}")
            except Exception as e:
                print(f"  (cross section skipped: {e})")


def build_parser():
    p = argparse.ArgumentParser(
        prog="aqs_ibm",
        description="Run SSHH topological experiments on IBM Quantum hardware.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp, default_out):
        sp.add_argument("--N", type=int, default=6)
        sp.add_argument("--UA", type=float, default=0.01)
        sp.add_argument("--delta-U", dest="delta_U", type=_floats,
                        default=[0.0, 0.001, 0.01, 0.1, 1.0])
        sp.add_argument("--TA", type=float, default=1.0)
        sp.add_argument("--steps", type=int, default=40, help="Trotter intervals L")
        sp.add_argument("--shots", type=int, default=8192)
        sp.add_argument("--backend", default=ibm.DEFAULT_BACKEND)
        sp.add_argument("--opt", type=int, default=3, help="optimization_level")
        sp.add_argument("--out", default=default_out)
        sp.add_argument("--dry-run", action="store_true",
                        help="transpile and report sizes; submit nothing")

    b = sub.add_parser("berry", help="Berry phase vs w (PBC)",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    b.add_argument("--electrons", type=int, default=12)
    b.add_argument("--v", type=float, default=1.0)
    b.add_argument("--w", type=_floats,
                   default=[0, 0.25, 0.5, 0.75, 0.99, 1.01, 1.25, 1.5, 1.75, 2.0])
    common(b, "ibm_results/berry")
    b.set_defaults(func=cmd_berry)

    pol = sub.add_parser("polarization", help="sublattice polarization vs dU (OBC)",
                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pol.add_argument("--electrons", type=int, default=14)
    pol.add_argument("--v", type=float, default=0.5)
    pol.add_argument("--w", type=float, default=1.5)
    common(pol, "ibm_results/polarization")
    pol.set_defaults(func=cmd_polarization)

    st = sub.add_parser("status", help="check job status from a manifest")
    st.add_argument("--manifest", required=True)
    st.set_defaults(func=cmd_status)

    ft = sub.add_parser("fetch", help="retrieve results and post-process")
    ft.add_argument("--manifest", required=True)
    ft.add_argument("--out", default=None)
    ft.add_argument("--plot", action="store_true")
    ft.add_argument("--save-counts", action="store_true",
                    help="embed raw counts in the JSON (large files)")
    ft.set_defaults(func=cmd_fetch)

    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    args.func(args)