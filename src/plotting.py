"""
plotting.py - Publication-quality figures for every experiment.

  * plot_fidelity        : Trotter-convergence curves (fidelity vs steps).
  * plot_berry_phase     : Berry phase vs w with a topology-aware solid/dotted
                           split at the gap-closing threshold.
  * plot_polarization    : sublattice polarization vs unit cell per Delta U.

Accept either in-memory records (from experiments.py) or a directory of JSON
files, so figures can be produced inline or from saved data.
"""

from __future__ import annotations

import glob
import json
import os

import matplotlib
matplotlib.use("Agg")  # headless-safe; figures are saved, not shown
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np


def _save(fig, out_path):
    d = os.path.dirname(out_path)
    if d and not os.path.exists(d):
        os.makedirs(d)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------
def plot_fidelity(result, out_path):
    curves = result["curves"]
    steps_list = result["params"]["steps_list"]
    items = sorted(curves.items(), key=lambda kv: float(kv[0]))
    fig, ax = plt.subplots(figsize=(8, 4))
    colors = cm.viridis(np.linspace(0, 0.9, max(len(items), 1)))
    for idx, (T_A, fids) in enumerate(items):
        ax.plot(steps_list, fids, "o-", color=colors[idx], label=fr"$T = {float(T_A):g}$")
    ax.axhline(1.0, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("Steps", fontsize=15)
    ax.set_ylabel("Fidelity", fontsize=15)
    ax.tick_params(axis="both", labelsize=12)
    ax.legend(loc="best", frameon=True, fontsize=13)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    _save(fig, out_path)
    return out_path


# ---------------------------------------------------------------------
def _w_of(d):
    return d.get("w", d.get("t2"))


def _load_berry_records(source):
    if isinstance(source, (list, tuple)):
        return list(source)
    records = []
    for fp in glob.glob(os.path.join(source, "*.json")):
        with open(fp, "r", encoding="utf-8") as f:
            data = json.load(f).get("data", {})
            if data and (_w_of(data) is not None):
                records.append(data)
    return records


def plot_berry_phase(source, out_path, v_default=1.0, fine_step=0.01):
    records = _load_berry_records(source)
    if not records:
        raise ValueError(f"No Berry-phase records found in: {source!r}")
    groups = {}
    for d in records:
        groups.setdefault(round(float(d.get("delta_U", 0.0)), 6), []).append(d)
    dus = sorted(groups.keys())
    fig, ax = plt.subplots(figsize=(8, 4))
    colors = cm.plasma(np.linspace(0, 0.9, max(len(dus), 1)))
    for idx, du in enumerate(dus):
        pts = sorted(groups[du], key=lambda d: _w_of(d))
        w = np.array([_w_of(d) for d in pts])
        z = np.array([d["Z_N_real"] + 1j * d["Z_N_imag"] for d in pts])
        raw = np.imag(np.log(z))
        g = (np.mod(raw + np.pi / 2.0, 2.0 * np.pi) - np.pi / 2.0) / np.pi
        v_val = float(pts[0].get("v", pts[0].get("t1", v_default)))
        ub_val = float(pts[0].get("U_B", 0.0))
        ax.plot([], [], "o-", color=colors[idx], markersize=4, linewidth=1, label=fr"$\Delta U={du:g}$")
        ax.plot(w, g, "o", color=colors[idx], markersize=4)
        if len(w) >= 2:
            wf = np.arange(w.min(), w.max() + 1e-3, fine_step)
            gf = np.interp(wf, w, g)
            gap = np.minimum(np.abs(v_val + wf), np.abs(v_val - wf))
            solid = ub_val < gap
            start, cur = 0, "-" if solid[0] else ":"
            for i in range(1, len(wf)):
                st = "-" if solid[i] else ":"
                if st != cur:
                    ax.plot(wf[start:i + 1], gf[start:i + 1], cur, color=colors[idx], linewidth=1)
                    cur, start = st, i
            ax.plot(wf[start:], gf[start:], cur, color=colors[idx], linewidth=1)
    ax.set_xlabel("Inter-cell hopping amplitude $w$", fontsize=15)
    ax.set_ylabel(r"Berry Phase $\gamma$ ($\pi$)", fontsize=15)
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis="both", labelsize=12)
    ax.legend(loc="upper right", fontsize=11, ncol=2)
    fig.tight_layout()
    _save(fig, out_path)
    return out_path


# ---------------------------------------------------------------------
def _load_polarization_records(source):
    if isinstance(source, (list, tuple)):
        return list(source)
    records = []
    for fp in sorted(glob.glob(os.path.join(source, "*.json"))):
        with open(fp, "r", encoding="utf-8") as f:
            rec = json.load(f)
            if "data" in rec and "n_A_minus_n_B" in rec.get("data", {}):
                records.append(rec)
    return records


def plot_polarization(source, out_path):
    records = _load_polarization_records(source)
    if not records:
        raise ValueError(f"No polarization records found in: {source!r}")
    records.sort(key=lambda r: float(r["metadata"].get("delta_U", 0.0)))
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = cm.plasma(np.linspace(0, 0.9, max(len(records), 1)))
    j_list = records[0]["data"]["unit_cell_j"]
    for idx, rec in enumerate(records):
        du = float(rec["metadata"].get("delta_U", 0.0))
        ax.plot(rec["data"]["unit_cell_j"], rec["data"]["n_A_minus_n_B"],
                "o-", color=colors[idx], markersize=8, linewidth=2, label=fr"$\Delta U={du:g}$")
    ax.axhline(0.0, color="gray", linestyle="--", alpha=0.6)
    ax.set_xlabel("Unit Cell Index $j$", fontsize=15)
    ax.set_ylabel(r"Polarization $\langle n_{A,j}\rangle-\langle n_{B,j}\rangle$", fontsize=15)
    ax.set_xticks(j_list)
    ax.tick_params(axis="both", labelsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=12)
    fig.tight_layout()
    _save(fig, out_path)
    return out_path
