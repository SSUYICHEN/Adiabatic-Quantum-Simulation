"""Smoke + contract tests for the experiment runners, plotting and the CLI.

These use the smallest systems that still exercise every branch, so the whole
file stays fast. The point is not to re-verify the physics (the dedicated test
modules do that) but to catch wiring breakage: JSON envelopes changing shape,
plotters rejecting the records their own runners emit, CLI arguments drifting.
"""
import json
import os
import tempfile

import numpy as np
import pytest

from aqs import experiments as ex
from aqs import plotting as pl
from aqs.cli import build_parser, main


# ----------------------------------------------------------------- runners
def test_fidelity_scan_envelope_and_json():
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "sub", "f.json")
        res = ex.fidelity_scan(N_cells=2, v=0.5, w=1.5, number_of_electrons=4,
                               U=0.0, T_A_values=(1,), steps_list=(2, 4),
                               out_json=out)
        assert res["experiment"] == "fidelity" and res["backend"] == "qiskit"
        assert res["params"]["steps_list"] == [2, 4]
        assert list(res["curves"]) == [1.0]
        assert len(res["curves"][1.0]) == 2
        assert all(0.0 <= f <= 1.0 + 1e-9 for f in res["curves"][1.0])
        assert os.path.exists(out), "out_json should create missing directories"
        with open(out) as f:
            assert json.load(f)["experiment"] == "fidelity"


def test_fidelity_scan_adjacent_reference_differs_from_ground():
    kw = dict(N_cells=2, v=0.5, w=1.5, number_of_electrons=4, U=0.3,
              T_A_values=(2,), steps_list=(2, 6))
    a = ex.fidelity_scan(reference="ground", **kw)["curves"][2.0]
    b = ex.fidelity_scan(reference="adjacent", **kw)["curves"][2.0]
    assert a != b


def test_berry_sweep_records_and_files():
    with tempfile.TemporaryDirectory() as d:
        recs = ex.berry_sweep(N_cells=2, number_of_electrons=4, v=1.0,
                              w_values=(0.5, 1.5), U_A=0.01,
                              delta_U_values=(0.0, 0.1), T_A=1.0, steps=4,
                              out_dir=d)
        assert len(recs) == 4
        for r in recs:
            # 't1'/'t2' kept as backward-compatible aliases for v/w
            assert r["t1"] == r["v"] and r["t2"] == r["w"]
            assert r["Twist_Amplitude"] == pytest.approx(
                abs(complex(r["Z_N_real"], r["Z_N_imag"])))
            assert r["delta_U"] == pytest.approx(r["U_B"] - r["U_A"])
        assert len([f for f in os.listdir(d) if f.endswith(".json")]) == 4


def test_polarization_sweep_records_and_files():
    with tempfile.TemporaryDirectory() as d:
        recs = ex.polarization_sweep(N_cells=2, total_electrons=6, v=0.5, w=1.5,
                                     U_A=1.0, delta_U_values=(0.0, 0.5),
                                     T_A=1.0, steps=4, out_dir=d)
        assert len(recs) == 2
        for r in recs:
            data, meta = r["data"], r["metadata"]
            assert meta["PBC"] is False, "polarization is an OBC diagnostic"
            assert data["unit_cell_j"] == [0, 1]
            assert data["n_A_minus_n_B"] == pytest.approx(
                [a - b for a, b in zip(data["n_A"], data["n_B"])], abs=1e-9)
        assert len(os.listdir(d)) == 2


def test_resolve_backend_accepts_an_object_or_a_name():
    from aqs.backends import get_backend
    bk = get_backend("qiskit")
    a = ex.fidelity_scan(N_cells=2, v=0.5, w=1.5, number_of_electrons=4,
                         T_A_values=(1,), steps_list=(2,), backend=bk)
    b = ex.fidelity_scan(N_cells=2, v=0.5, w=1.5, number_of_electrons=4,
                         T_A_values=(1,), steps_list=(2,), backend="qiskit")
    assert a["curves"] == b["curves"]


# ---------------------------------------------------------------- plotting
def test_plotters_accept_the_records_their_runners_emit():
    with tempfile.TemporaryDirectory() as d:
        fid = ex.fidelity_scan(N_cells=2, v=0.5, w=1.5, number_of_electrons=4,
                               T_A_values=(1,), steps_list=(2, 4))
        berry = ex.berry_sweep(N_cells=2, number_of_electrons=4, v=1.0,
                               w_values=(0.5, 1.5), delta_U_values=(0.0, 0.1),
                               T_A=1.0, steps=4)
        pol = ex.polarization_sweep(N_cells=2, total_electrons=6, v=0.5, w=1.5,
                                    U_A=1.0, delta_U_values=(0.0, 0.5),
                                    T_A=1.0, steps=4)
        for fn, arg, name in [
            (pl.plot_fidelity, fid, "f.png"),
            (pl.plot_berry_phase, berry, "b.png"),
            (pl.plot_polarization, pol, "p.png"),
            (lambda s, o: pl.plot_berry_cross_section(s, o, w_cut=1.5), berry, "bc.png"),
            (lambda s, o: pl.plot_polarization_cross_section(s, o, cell_index=0), pol, "pc.png"),
        ]:
            out = os.path.join(d, name)
            fn(arg, out)
            assert os.path.getsize(out) > 1000, name


def test_plotters_can_read_back_a_directory_of_json():
    with tempfile.TemporaryDirectory() as d:
        bd, pd = os.path.join(d, "b"), os.path.join(d, "p")
        ex.berry_sweep(N_cells=2, number_of_electrons=4, v=1.0,
                       w_values=(0.5, 1.5), delta_U_values=(0.0,),
                       T_A=1.0, steps=4, out_dir=bd)
        ex.polarization_sweep(N_cells=2, total_electrons=6, v=0.5, w=1.5,
                              U_A=1.0, delta_U_values=(0.0,), T_A=1.0,
                              steps=4, out_dir=pd)
        pl.plot_berry_phase(bd, os.path.join(d, "b.png"))
        pl.plot_polarization(pd, os.path.join(d, "p.png"))
        assert os.path.getsize(os.path.join(d, "b.png")) > 1000


def test_plotters_raise_on_empty_input():
    with tempfile.TemporaryDirectory() as d:
        with pytest.raises(ValueError, match="No Berry-phase records"):
            pl.plot_berry_phase(d, os.path.join(d, "x.png"))
        with pytest.raises(ValueError, match="No polarization records"):
            pl.plot_polarization(d, os.path.join(d, "y.png"))
        with pytest.raises(ValueError, match="No Berry-phase records"):
            pl.plot_berry_cross_section([], os.path.join(d, "z.png"), w_cut=9.0)


# --------------------------------------------------------------------- CLI
def test_parser_argument_types():
    p = build_parser()
    a = p.parse_args(["berry", "--N", "3", "--w", "0.5,1.0", "--delta-U", "0,0.1"])
    assert a.N == 3 and a.w == [0.5, 1.0] and a.delta_U == [0.0, 0.1]
    a = p.parse_args(["fidelity", "--steps", "1,2,3"])
    assert a.steps == [1, 2, 3] and a.pbc is True
    assert p.parse_args(["fidelity", "--obc"]).pbc is False
    a = p.parse_args(["measure", "--hamiltonian", "h.json",
                      "--property", "berry,polarization"])
    assert a.property == ["berry", "polarization"]


def test_parser_requires_a_subcommand():
    with pytest.raises(SystemExit):
        build_parser().parse_args([])


def _run(argv, cwd):
    old = os.getcwd()
    try:
        os.chdir(cwd)
        main(argv)
    finally:
        os.chdir(old)


def test_cli_fidelity_end_to_end():
    with tempfile.TemporaryDirectory() as d:
        _run(["fidelity", "--N", "2", "--electrons", "4", "--TA", "1",
              "--steps", "2,4"], d)
        base = os.path.join(d, "aqs_results", "fidelity")
        assert any(f.endswith(".json") for f in os.listdir(base))
        assert any(f.endswith(".png") for f in os.listdir(base))


def test_cli_berry_and_polarization_end_to_end():
    with tempfile.TemporaryDirectory() as d:
        _run(["berry", "--N", "2", "--electrons", "4", "--w", "0.5,1.5",
              "--delta-U", "0,0.1", "--steps", "4"], d)
        _run(["polarization", "--N", "2", "--electrons", "6", "--delta-U",
              "0,0.5", "--steps", "4"], d)
        for kind in ("berry", "polarization"):
            root = os.path.join(d, "aqs_results", kind)
            pngs = [f for r, _, fs in os.walk(root) for f in fs if f.endswith(".png")]
            assert pngs, kind


def test_cli_no_plot_suppresses_figures():
    with tempfile.TemporaryDirectory() as d:
        _run(["berry", "--N", "2", "--electrons", "4", "--w", "1.5",
              "--steps", "4", "--no-plot"], d)
        root = os.path.join(d, "aqs_results", "berry")
        pngs = [f for r, _, fs in os.walk(root) for f in fs if f.endswith(".png")]
        assert not pngs


def test_cli_measure_and_replot():
    src = os.path.join(os.getcwd(), "examples", "ssh_spinless_N3_topological.json")
    if not os.path.exists(src):
        pytest.skip("example file not present")
    with tempfile.TemporaryDirectory() as d:
        _run(["measure", "--hamiltonian", src, "--property",
              "berry,polarization"], d)
        out = os.path.join(d, "aqs_results", "measure")
        assert any(f.endswith(".json") for f in os.listdir(out))


def test_cli_selftest_passes():
    with tempfile.TemporaryDirectory() as d:
        _run(["selftest", "--backend", "qiskit"], d)      # raises SystemExit on failure


def test_polarization_cross_section_raises_on_empty_input():
    with tempfile.TemporaryDirectory() as d:
        with pytest.raises(ValueError, match="No polarization records"):
            pl.plot_polarization_cross_section([], os.path.join(d, "x.png"))
