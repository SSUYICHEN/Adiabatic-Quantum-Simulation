"""
aqs_ddsim.py - Command line front-end for running the SSHH experiments on the
local MQT DDSIM decision-diagram simulator.

Unlike the cloud backends this is synchronous and free: one command runs the
whole sweep and writes plots. No auth, no manifest, no fetch.

  # N=6 cross-check against your statevector / paper results
  python aqs_ddsim.py berry        --N 6 --steps 20 --plot
  python aqs_ddsim.py polarization --N 6 --steps 20 --plot

  # push to larger systems (decision diagrams may fit where statevector can't)
  python aqs_ddsim.py berry        --N 8  --steps 20 --plot
  python aqs_ddsim.py polarization --N 10 --steps 20 --plot

Defaults follow the study settings:
  berry        : UA=0.01, dU in {0,0.01,0.1,1}, v=1,   w=0..2, n=2N,   PBC
  polarization : UA=0.01, dU in {0,0.01,0.1,1}, v=0.5, w=1.5,  n=2N+2, OBC
Electron count defaults to 2N (berry) / 2N+2 (polarization) unless --electrons.
"""

import argparse

from aqs import ddsim_run as dd


def _floats(s):
    return [float(x) for x in str(s).replace(" ", "").split(",") if x != ""]


def cmd_berry(a):
    records, data_dir = dd.run_berry(
        N_cells=a.N, electrons=a.electrons, v=a.v, w_values=a.w,
        U_A=a.UA, delta_U_values=a.delta_U, T_A=a.TA, steps=a.steps,
        shots=a.shots, out_dir=a.out)
    if a.plot and records:
        from aqs import plotting as pl
        img = f"{data_dir}/BerryPhase_ddsim_N{a.N}.png"
        pl.plot_berry_phase(data_dir, img, v_default=a.v)
        print(f"  figure -> {img}")


def cmd_polarization(a):
    records, data_dir = dd.run_polarization(
        N_cells=a.N, electrons=a.electrons, v=a.v, w=a.w,
        U_A=a.UA, delta_U_values=a.delta_U, T_A=a.TA, steps=a.steps,
        shots=a.shots, out_dir=a.out)
    if a.plot and records:
        from aqs import plotting as pl
        img = f"{data_dir}/Polarization_ddsim_N{a.N}.png"
        pl.plot_polarization(data_dir, img)
        try:
            pl.plot_polarization_cross_section(
                data_dir, f"{data_dir}/Polarization_CrossSection_ddsim_N{a.N}.png",
                cell_index=0)
        except Exception as e:
            print(f"  (cross section skipped: {e})")
        print(f"  figure -> {img}")


def build_parser():
    p = argparse.ArgumentParser(
        prog="aqs_ddsim",
        description="Run SSHH experiments on the local MQT DDSIM simulator.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp, default_out):
        sp.add_argument("--N", type=int, default=6)
        sp.add_argument("--electrons", type=int, default=None,
                        help="default: 2N (berry) / 2N+2 (polarization)")
        sp.add_argument("--UA", type=float, default=0.01)
        sp.add_argument("--delta-U", dest="delta_U", type=_floats,
                        default=[0.0, 0.01, 0.1, 1.0])
        sp.add_argument("--TA", type=float, default=1.0)
        sp.add_argument("--steps", type=int, default=20, help="Trotter intervals L")
        sp.add_argument("--shots", type=int, default=8192)
        sp.add_argument("--out", default=default_out)
        sp.add_argument("--plot", action="store_true")

    b = sub.add_parser("berry", formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                       help="Berry phase vs w (PBC)")
    b.add_argument("--v", type=float, default=1.0)
    b.add_argument("--w", type=_floats,
                   default=[0, 0.25, 0.5, 0.75, 0.99, 1.01, 1.25, 1.5, 1.75, 2.0])
    common(b, "ddsim_results/berry")
    b.set_defaults(func=cmd_berry)

    pol = sub.add_parser("polarization", formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                         help="sublattice polarization vs dU (OBC)")
    pol.add_argument("--v", type=float, default=0.5)
    pol.add_argument("--w", type=float, default=1.5)
    common(pol, "ddsim_results/polarization")
    pol.set_defaults(func=cmd_polarization)

    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    args.func(args)