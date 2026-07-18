"""
plotting.py - Publication-quality figures for every experiment.

  * plot_fidelity           : Trotter-convergence curves (fidelity vs L).
  * plot_berry_phase        : Berry phase vs w with a topology-aware
                              solid/dotted split at the gap-closing threshold.
  * plot_polarization       : sublattice polarization vs unit cell per Delta U.
  * plot_berry_cross_section        : Berry phase vs Delta U at fixed w
                                      (cf. Fig. 3 bottom of the paper).
  * plot_polarization_cross_section : edge polarization vs Delta U
                                      (cf. Fig. 4 bottom of the paper).

The two cross-section figures colour each point with the SAME colour that its
Delta U curve has in the corresponding phase-diagram figure (shared plasma
colour map via _delta_u_colors), so top and bottom panels match visually.

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
    ax.set_xlabel(r"Number of Trotter intervals $L$", fontsize=15)
    ax.set_ylabel("Fidelity", fontsize=15)
    ax.tick_params(axis="both", labelsize=12)
    ax.legend(loc="best", frameon=True, fontsize=13)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    _save(fig, out_path)
    return out_path


# ---------------------------------------------------------------------
def _delta_u_colors(delta_us):
    """One colour per Delta U value, shared by the phase-diagram figures and
    the cross-section figures so that identical Delta U -> identical colour."""
    dus = sorted({round(float(du), 6) for du in delta_us})
    colors = cm.plasma(np.linspace(0, 0.9, max(len(dus), 1)))
    return {du: colors[i] for i, du in enumerate(dus)}


def _symlog_deltaU_axis(ax, dus, linthresh=1e-3):
    ax.set_xscale("symlog", linthresh=linthresh)
    lo = min(dus); hi = max(dus)
    ax.set_xlim(min(-0.1 * linthresh, lo - 0.1 * linthresh), max(hi * 3.0, linthresh))
    ax.set_xlabel(r"Hubbard Interaction $\Delta U$", fontsize=15)
    ax.tick_params(axis="both", labelsize=12)
    ax.grid(True, which="both", alpha=0.3)


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
    cmap = _delta_u_colors(dus)
    colors = [cmap[du] for du in dus]
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
    cmap = _delta_u_colors(float(r["metadata"].get("delta_U", 0.0)) for r in records)
    j_list = records[0]["data"]["unit_cell_j"]
    for rec in records:
        du = round(float(rec["metadata"].get("delta_U", 0.0)), 6)
        ax.plot(rec["data"]["unit_cell_j"], rec["data"]["n_A_minus_n_B"],
                "o-", color=cmap[du], markersize=8, linewidth=2, label=fr"$\Delta U={du:g}$")
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


# ---------------------------------------------------------------------
# Cross-section figures (bottom panels of Fig. 3 / Fig. 4 in the paper)
# ---------------------------------------------------------------------
def plot_berry_cross_section(source, out_path, w_cut=1.5, linthresh=1e-3,
                             w_tol=1e-6):
    """Berry phase gamma(pi) vs Delta U at a fixed inter-cell hopping w=w_cut
    (cf. Fig. 3, bottom panel). Each marker is drawn with the SAME colour that
    its Delta U curve has in plot_berry_phase, so the two figures correspond
    point-for-point."""
    records = _load_berry_records(source)
    pts = [d for d in records if abs(float(_w_of(d)) - w_cut) < w_tol]
    if not pts:
        raise ValueError(f"No Berry-phase records found at w={w_cut!r} in: {source!r}")
    # colour map built from ALL Delta U present in the dataset -> matches the
    # phase-diagram legend colours exactly.
    cmap = _delta_u_colors(float(d.get("delta_U", 0.0)) for d in records)

    pts.sort(key=lambda d: float(d.get("delta_U", 0.0)))
    dus = [round(float(d.get("delta_U", 0.0)), 6) for d in pts]
    z = np.array([d["Z_N_real"] + 1j * d["Z_N_imag"] for d in pts])
    raw = np.imag(np.log(z))
    gamma = (np.mod(raw + np.pi / 2.0, 2.0 * np.pi) - np.pi / 2.0) / np.pi

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(dus, gamma, "-", color="gray", alpha=0.5, linewidth=1.5, zorder=1)
    ax.scatter(dus, gamma, c=[cmap[du] for du in dus], s=45, zorder=2,
               edgecolors="none")
    ax.axhline(0.0, color="gray", linestyle="--", alpha=0.5)
    _symlog_deltaU_axis(ax, dus, linthresh=linthresh)
    ax.set_ylabel(r"Berry Phase $\gamma$ ($\pi$)", fontsize=15)
    ax.set_title(fr"Cross section at $w = {w_cut:g}$", fontsize=13)
    fig.tight_layout()
    _save(fig, out_path)
    return out_path


def plot_polarization_cross_section(source, out_path, cell_index=0,
                                    linthresh=1e-3):
    """Sublattice polarization <n_A,j> - <n_B,j> at one unit cell (default
    j=0, the A-sublattice edge) vs Delta U (cf. Fig. 4, bottom panel). Marker
    colours match the Delta U curve colours of plot_polarization."""
    records = _load_polarization_records(source)
    if not records:
        raise ValueError(f"No polarization records found in: {source!r}")
    records.sort(key=lambda r: float(r["metadata"].get("delta_U", 0.0)))
    cmap = _delta_u_colors(float(r["metadata"].get("delta_U", 0.0))
                           for r in records)

    dus, vals = [], []
    for rec in records:
        dus.append(round(float(rec["metadata"].get("delta_U", 0.0)), 6))
        vals.append(float(rec["data"]["n_A_minus_n_B"][cell_index]))

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(dus, vals, "-", color="gray", alpha=0.5, linewidth=1.5, zorder=1)
    ax.scatter(dus, vals, c=[cmap[du] for du in dus], s=45, zorder=2,
               edgecolors="none")
    _symlog_deltaU_axis(ax, dus, linthresh=linthresh)
    ax.set_ylabel(
        fr"Polarization $\langle n_{{A,{cell_index}}}\rangle-"
        fr"\langle n_{{B,{cell_index}}}\rangle$", fontsize=15)
    ax.set_title(fr"Cross section at unit cell $j = {cell_index}$", fontsize=13)
    fig.tight_layout()
    _save(fig, out_path)
    return out_path
