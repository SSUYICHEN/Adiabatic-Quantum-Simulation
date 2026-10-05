"""Tests for `aqs plot` (re-drawing figures from saved JSON) and entry points."""
import json
import os
import subprocess
import sys
import tempfile

import pytest

from aqs import experiments as ex
from aqs.cli import main


def _run(argv, cwd):
    old = os.getcwd()
    try:
        os.chdir(cwd)
        main(argv)
    finally:
        os.chdir(old)


@pytest.fixture()
def saved(tmp_path):
    """A directory tree of saved runs for the plot subcommand to read back."""
    d = str(tmp_path)
    fj = os.path.join(d, "fid.json")
    res = ex.fidelity_scan(N_cells=2, v=0.5, w=1.5, number_of_electrons=4,
                           T_A_values=(1,), steps_list=(2, 4), out_json=fj)
    bd, pd = os.path.join(d, "berry"), os.path.join(d, "pol")
    ex.berry_sweep(N_cells=2, number_of_electrons=4, v=1.0,
                   w_values=(0.5, 1.5), delta_U_values=(0.0, 0.1),
                   T_A=1.0, steps=4, out_dir=bd)
    ex.polarization_sweep(N_cells=2, total_electrons=6, v=0.5, w=1.5, U_A=1.0,
                          delta_U_values=(0.0, 0.5), T_A=1.0, steps=4, out_dir=pd)
    return {"dir": d, "fidelity_json": fj, "berry_dir": bd, "pol_dir": pd}


def test_plot_fidelity_from_json_defaults_to_sibling_png(saved):
    _run(["plot", "--kind", "fidelity", "--data", saved["fidelity_json"]],
         saved["dir"])
    png = os.path.splitext(saved["fidelity_json"])[0] + ".png"
    assert os.path.getsize(png) > 1000


def test_plot_berry_from_directory_defaults_to_kind_png(saved):
    _run(["plot", "--kind", "berry", "--data", saved["berry_dir"]], saved["dir"])
    assert os.path.getsize(os.path.join(saved["berry_dir"], "berry.png")) > 1000


def test_plot_polarization_from_directory(saved):
    _run(["plot", "--kind", "polarization", "--data", saved["pol_dir"]],
         saved["dir"])
    assert os.path.getsize(os.path.join(saved["pol_dir"], "polarization.png")) > 1000


def test_plot_cross_sections(saved):
    out_b = os.path.join(saved["dir"], "bc.png")
    out_p = os.path.join(saved["dir"], "pc.png")
    _run(["plot", "--kind", "berry-cut", "--data", saved["berry_dir"],
          "--w-cut", "1.5", "--out", out_b], saved["dir"])
    _run(["plot", "--kind", "polarization-cut", "--data", saved["pol_dir"],
          "--cell", "0", "--out", out_p], saved["dir"])
    assert os.path.getsize(out_b) > 1000 and os.path.getsize(out_p) > 1000


def test_plot_explicit_out_path_is_honoured(saved):
    out = os.path.join(saved["dir"], "nested", "custom.png")
    _run(["plot", "--kind", "fidelity", "--data", saved["fidelity_json"],
          "--out", out], saved["dir"])
    assert os.path.getsize(out) > 1000, "_save must create missing directories"


def test_plot_rejects_an_unknown_kind(saved):
    with pytest.raises(SystemExit):
        _run(["plot", "--kind", "nope", "--data", saved["berry_dir"]],
             saved["dir"])


# ------------------------------------------------------------ entry points
def test_module_entry_point_runs():
    """python -m aqs --version must work (covers __main__.py)."""
    r = subprocess.run([sys.executable, "-m", "aqs", "--version"],
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr
    assert "aqs" in (r.stdout + r.stderr).lower()


def test_selftest_failure_exits_nonzero(monkeypatch):
    """cmd_selftest must raise SystemExit(1) when the backend disagrees."""
    import aqs.cli as cli
    monkeypatch.setattr(cli, "run_selftest", lambda name: (False, 0.5))
    with tempfile.TemporaryDirectory() as d:
        with pytest.raises(SystemExit) as e:
            _run(["selftest", "--backend", "qiskit"], d)
        assert e.value.code == 1
