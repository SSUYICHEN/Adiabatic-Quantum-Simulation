"""
cli.py - Command-line interface for the AQS toolkit.

Built-in SSH-Hubbard model:
  aqs fidelity      [options]      # choose T_A and steps (L)
  aqs berry         [options]      # Berry phase vs w (PBC)
  aqs polarization  [options]      # polarization vs Delta U (OBC)

Arbitrary Hamiltonian file:
  aqs measure --hamiltonian H.json --property berry,polarization

Utilities:
  aqs plot --kind ... --data ...   # (re)draw from saved JSON
  aqs selftest [--backend cudaq]   # verify a backend against qiskit

Add --backend cudaq to any measurement command to run on an NVIDIA GPU.
Run `aqs <command> -h` for each command's options.
"""

from __future__ import annotations

import argparse
import json
import os

from . import __version__
from . import experiments as ex
from . import plotting as pl
from .backends import selftest as run_selftest

DEFAULT_OUT = "aqs_results"


def _floats(s):
    return [float(x) for x in str(s).replace(" ", "").split(",") if x != ""]


def _ints(s):
    return [int(x) for x in str(s).replace(" ", "").split(",") if x != ""]


def _strs(s):
    return [x for x in str(s).replace(" ", "").split(",") if x != ""]


# ------------------------------ fidelity -----------------------------
def cmd_fidelity(a):
    out_json = a.out or os.path.join(DEFAULT_OUT, "fidelity",
                                     f"fidelity_N{a.N}_U{a.U:g}_{a.reference}.json")
    print(f"[fidelity] N={a.N} v={a.v} w={a.w} Ne={a.electrons} U={a.U} "
          f"{'PBC' if a.pbc else 'OBC'} ref={a.reference} backend={a.backend}")
    result = ex.fidelity_scan(N_cells=a.N, v=a.v, w=a.w, number_of_electrons=a.electrons,
                              U=a.U, T_A_values=a.TA, steps_list=a.steps, PBC=a.pbc,
                              reference=a.reference, out_json=out_json, backend=a.backend)
    print(f"  -> saved JSON: {out_json}")
    steps = result["params"]["steps_list"]
    print("  steps : " + "  ".join(f"{s:>5d}" for s in steps))
    best = None
    for T_A, fids in sorted(result["curves"].items()):
        print(f"  T={T_A:<4g}: " + "  ".join(f"{f:5.3f}" for f in fids))
        for s, f in zip(steps, fids):
            if f >= 0.999 and (best is None or s < best[1]):
                best = (T_A, s, f)
    if best:
        print(f"  suggestion: fidelity>=0.999 first at T_A={best[0]:g}, steps={best[1]}")
    if not a.no_plot:
        img = os.path.splitext(out_json)[0] + ".png"
        pl.plot_fidelity(result, img)
        print(f"  -> saved figure: {img}")


# ------------------------------- berry -------------------------------
def cmd_berry(a):
    out_dir = a.out or os.path.join(DEFAULT_OUT, "berry",
                                    f"N{a.N}_TA{a.TA:g}_steps{a.steps}_UA{a.UA:g}_{a.backend}")
    print(f"[berry] N={a.N} Ne={a.electrons} v={a.v} UA={a.UA} dU={a.delta_U} "
          f"TA={a.TA} steps={a.steps} backend={a.backend}")
    records = ex.berry_sweep(N_cells=a.N, number_of_electrons=a.electrons, v=a.v,
                             w_values=a.w, U_A=a.UA, delta_U_values=a.delta_U,
                             T_A=a.TA, steps=a.steps, mode=a.mode, out_dir=out_dir,
                             backend=a.backend)
    print(f"  -> {len(records)} points saved to: {out_dir}")
    for d in records:
        print(f"    dU={d['delta_U']:<7g} w={d['w']:<5g} |z|={d['Twist_Amplitude']:.4f} "
              f"gamma={d['Berry_Phase_pi_wrapped']:+.4f} pi")
    if not a.no_plot:
        img = os.path.join(out_dir, "BerryPhase_PhaseDiagram.png")
        pl.plot_berry_phase(records, img, v_default=a.v)
        print(f"  -> saved figure: {img}")
        # cross section gamma vs Delta U at fixed w (Fig. 3 bottom); colours
        # match the corresponding points of the phase diagram above.
        w_cut = a.w_cut if a.w_cut is not None else max(a.w)
        if len(a.delta_U) > 1 and any(abs(w - w_cut) < 1e-9 for w in a.w):
            img2 = os.path.join(out_dir, f"BerryPhase_CrossSection_w_{w_cut:g}.png")
            pl.plot_berry_cross_section(records, img2, w_cut=w_cut)
            print(f"  -> saved figure: {img2}")


# --------------------------- polarization ----------------------------
def cmd_polarization(a):
    electrons = a.electrons if a.electrons is not None else 2 * a.N + 2
    out_dir = a.out or os.path.join(DEFAULT_OUT, "polarization",
                                    f"N{a.N}_TA{a.TA:g}_steps{a.steps}_v{a.v:g}_{a.backend}")
    print(f"[polarization] N={a.N} Ne={electrons} v={a.v} w={a.w} UA={a.UA} "
          f"dU={a.delta_U} TA={a.TA} steps={a.steps} backend={a.backend} (OBC)")
    records = ex.polarization_sweep(N_cells=a.N, total_electrons=electrons, v=a.v, w=a.w,
                                    U_A=a.UA, delta_U_values=a.delta_U, T_A=a.TA,
                                    steps=a.steps, out_dir=out_dir, backend=a.backend)
    print(f"  -> {len(records)} curves saved to: {out_dir}")
    for r in records:
        print(f"    dU={r['metadata']['delta_U']:<7g} "
              f"edge(cell 0)={r['data']['n_A_minus_n_B'][0]:+.4f}")
    if not a.no_plot:
        img = os.path.join(out_dir, "Polarization_vs_DeltaU.png")
        pl.plot_polarization(records, img)
        print(f"  -> saved figure: {img}")
        # cross section: edge polarization (unit cell a.cell) vs Delta U
        # (Fig. 4 bottom); colours match the curves of the figure above.
        if len(a.delta_U) > 1:
            img2 = os.path.join(out_dir, f"Polarization_CrossSection_cell_{a.cell}.png")
            pl.plot_polarization_cross_section(records, img2, cell_index=a.cell)
            print(f"  -> saved figure: {img2}")


# ------------------------------ measure (arbitrary H) ----------------
def cmd_measure(a):
    from .backends import get_backend
    from .hamiltonian import measure_hamiltonian
    props = a.property
    out_json = a.out or os.path.join(
        DEFAULT_OUT, "measure",
        os.path.splitext(os.path.basename(a.hamiltonian))[0] + f"_{a.backend}.json")
    print(f"[measure] file={a.hamiltonian} property={props} backend={a.backend}")
    result = measure_hamiltonian(a.hamiltonian, props, get_backend(a.backend))
    os.makedirs(os.path.dirname(out_json), exist_ok=True)
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=4)
    print(f"  -> saved JSON: {out_json}")
    for name, res in result["results"].items():
        if name == "berry":
            print(f"    berry: |z|={res['Twist_Amplitude']:.4f} "
                  f"gamma={res['Berry_Phase_pi_wrapped']:+.4f} pi")
        elif name == "polarization":
            print(f"    polarization per cell (n_A - n_B): {res['n_A_minus_n_B']}")
    # optional polarization figure
    if not a.no_plot and "polarization" in result["results"]:
        rec = {"metadata": {"delta_U": result["metadata"].get("delta_U", 0.0)},
               "data": result["results"]["polarization"]}
        img = os.path.splitext(out_json)[0] + "_polarization.png"
        pl.plot_polarization([rec], img)
        print(f"  -> saved figure: {img}")


# ------------------------------ selftest -----------------------------
def cmd_selftest(a):
    print(f"[selftest] backend={a.backend} ...")
    ok, fid = run_selftest(a.backend)
    status = "PASS" if ok else "FAIL"
    print(f"  statevector vs qiskit reference: fidelity={fid:.8f}  [{status}]")
    if not ok:
        raise SystemExit(1)


# ------------------------------- plot --------------------------------
def cmd_plot(a):
    out = a.out or (os.path.splitext(a.data)[0] + ".png" if a.data.endswith(".json")
                    else os.path.join(a.data, f"{a.kind}.png"))
    if a.kind == "fidelity":
        with open(a.data, "r", encoding="utf-8") as f:
            pl.plot_fidelity(json.load(f), out)
    elif a.kind == "berry":
        pl.plot_berry_phase(a.data, out, v_default=a.v)
    elif a.kind == "polarization":
        pl.plot_polarization(a.data, out)
    elif a.kind == "berry-cut":
        pl.plot_berry_cross_section(a.data, out, w_cut=a.w_cut)
    elif a.kind == "polarization-cut":
        pl.plot_polarization_cross_section(a.data, out, cell_index=a.cell)
    print(f"  -> saved figure: {out}")


# ------------------------------ parser -------------------------------
def _add_backend(p):
    p.add_argument("--backend", choices=["qiskit", "cudaq"], default="qiskit",
                   help="simulation backend (cudaq = NVIDIA GPU, Linux only)")


def build_parser():
    p = argparse.ArgumentParser(prog="aqs",
                                description="Adiabatic quantum simulation of topological properties.",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--version", action="version", version=f"aqs {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    # fidelity
    f = sub.add_parser("fidelity", help="Trotter-convergence check (choose T_A & steps)",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    f.add_argument("--N", type=int, default=6)
    f.add_argument("--v", type=float, default=0.5, help="intra-cell hopping v")
    f.add_argument("--w", type=float, default=1.5, help="inter-cell hopping w")
    f.add_argument("--electrons", type=int, default=12)
    f.add_argument("--U", type=float, default=0.0)
    f.add_argument("--TA", type=_floats, default=[1, 5, 15, 80])
    f.add_argument("--steps", type=_ints, default=[1, 10, 20, 30, 40, 50, 60, 80, 100, 120, 150])
    f.add_argument("--reference", choices=["ground", "adjacent"], default="ground")
    g = f.add_mutually_exclusive_group()
    g.add_argument("--pbc", dest="pbc", action="store_true", default=True)
    g.add_argument("--obc", dest="pbc", action="store_false")
    f.add_argument("--out"); f.add_argument("--no-plot", action="store_true")
    _add_backend(f); f.set_defaults(func=cmd_fidelity)

    # berry
    b = sub.add_parser("berry", help="Berry phase vs inter-cell hopping w (PBC)",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    b.add_argument("--N", type=int, default=6)
    b.add_argument("--electrons", type=int, default=12)
    b.add_argument("--v", type=float, default=1.0, help="intra-cell hopping v")
    b.add_argument("--w", type=_floats, default=[0, 0.25, 0.5, 0.75, 0.99, 1.01, 1.25, 1.5, 1.75, 2.0],
                   help="inter-cell hoppings w to scan (comma-separated)")
    b.add_argument("--UA", type=float, default=0.01)
    b.add_argument("--delta-U", dest="delta_U", type=_floats, default=[0.0])
    b.add_argument("--TA", type=float, default=1.0)
    b.add_argument("--steps", type=int, default=40)
    b.add_argument("--mode", choices=["up_spin", "total"], default="up_spin")
    b.add_argument("--w-cut", dest="w_cut", type=float, default=None,
                   help="w value for the gamma-vs-DeltaU cross-section figure "
                        "(default: the largest scanned w; must be in --w)")
    b.add_argument("--out"); b.add_argument("--no-plot", action="store_true")
    _add_backend(b); b.set_defaults(func=cmd_berry)

    # polarization
    pol = sub.add_parser("polarization", help="electron polarization vs Delta U (OBC)",
                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pol.add_argument("--N", type=int, default=6)
    pol.add_argument("--electrons", type=int, default=None, help="default = 2N+2")
    pol.add_argument("--v", type=float, default=0.5)
    pol.add_argument("--w", type=float, default=1.5)
    pol.add_argument("--UA", type=float, default=1.0)
    pol.add_argument("--delta-U", dest="delta_U", type=_floats,
                     default=[0, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3])
    pol.add_argument("--TA", type=float, default=1.0)
    pol.add_argument("--steps", type=int, default=40)
    pol.add_argument("--cell", type=int, default=0,
                     help="unit-cell index for the polarization-vs-DeltaU "
                          "cross-section figure (0 = A-sublattice edge)")
    pol.add_argument("--out"); pol.add_argument("--no-plot", action="store_true")
    _add_backend(pol); pol.set_defaults(func=cmd_polarization)

    # measure (arbitrary Hamiltonian)
    m = sub.add_parser("measure", help="measure topological properties of a Hamiltonian file",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    m.add_argument("--hamiltonian", required=True, help="path to a Hamiltonian JSON file")
    m.add_argument("--property", type=_strs, default=["berry"],
                   help="comma-separated: berry,polarization")
    m.add_argument("--out"); m.add_argument("--no-plot", action="store_true")
    _add_backend(m); m.set_defaults(func=cmd_measure)

    # selftest
    st = sub.add_parser("selftest", help="verify a backend's statevector against qiskit",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    _add_backend(st); st.set_defaults(func=cmd_selftest)

    # plot
    pt = sub.add_parser("plot", help="(re)draw a figure from saved JSON",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pt.add_argument("--kind", required=True,
                    choices=["fidelity", "berry", "polarization",
                             "berry-cut", "polarization-cut"])
    pt.add_argument("--data", required=True)
    pt.add_argument("--v", type=float, default=1.0, help="v for the berry gap threshold")
    pt.add_argument("--w-cut", dest="w_cut", type=float, default=1.5,
                    help="w value for --kind berry-cut")
    pt.add_argument("--cell", type=int, default=0,
                    help="unit-cell index for --kind polarization-cut")
    pt.add_argument("--out"); pt.set_defaults(func=cmd_plot)

    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
