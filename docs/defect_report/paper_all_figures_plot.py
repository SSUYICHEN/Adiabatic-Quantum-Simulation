"""Side-by-side 修正前 | 修正後 versions of every figure in the paper (Fig. 2, 3, 4),
drawn in the paper's own style (aqs.plotting colour maps / axes). Reads
paper_fig2.json and paper_impact.json next to this script; writes
figures/paper_fig{2,3,4}_compare.svg and figures/paper_all_compare.svg.
Run: uv run --no-sync python docs/defect_report/paper_all_figures_plot.py"""
import json, pathlib, os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.font_manager as fm
from aqs.plotting import _delta_u_colors, _symlog_deltaU_axis

_cjk = [f.name for f in fm.fontManager.ttflist if "CJK" in f.name and "Sans" in f.name]
matplotlib.rcParams["font.family"] = ([_cjk[0]] if _cjk else []) + ["DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

S = pathlib.Path(__file__).parent
FIG = S / "figures"
F2 = json.load(open(S / "paper_fig2.json"))
F34 = json.load(open(S / "paper_impact.json"))
MODE = "up_spin"   # the code default the paper's figures were produced with
C_BEFORE, C_AFTER = "#c2410c", "#1d5fb3"

def title(ax, txt, color):
    ax.set_title(txt, loc="left", fontsize=11, color=color, pad=8)

# ------------------------------------------------------------------ Fig. 2
def draw_fig2(axes, version):
    L = F2["meta"]["L"]; Ts = F2["meta"]["T"]
    colors = cm.viridis(np.linspace(0, 0.9, len(Ts)))
    for ax, (tag, ref) in zip(axes, (("U0", "ground"), ("U1", "adjacent"))):
        for c, T in zip(colors, Ts):
            y = F2[tag][version][str(T)][ref]
            ax.plot(L, y, "o-", color=c, ms=4, lw=1, label=fr"$T={T}$")
        ax.axhline(1.0, color="gray", ls="--", alpha=0.6); ax.grid(True, alpha=0.3)
        ax.set_ylabel("Fidelity"); ax.legend(fontsize=8, loc="lower right" if tag == "U0" else "center right")
    axes[1].set_xlabel("Steps $L$")

# ------------------------------------------------------------------ Fig. 3
def draw_fig3(axes, version):
    top = F34["fig3_top"]; dus = sorted({r["delta_U"] for r in top}); cmap = _delta_u_colors(dus)
    ax = axes[0]
    for du in dus:
        pts = sorted([(r["w"], r[version][MODE]["gamma_pi_wrapped"]) for r in top if r["delta_U"] == du])
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", color=cmap[round(du, 6)], ms=4, lw=1, label=fr"$\Delta U={du:g}$")
    ax.set_xlabel("Inter-cell hopping amplitude $w$"); ax.set_ylabel(r"Berry Phase $\gamma$ ($\pi$)")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=7.5, ncol=2, loc="center right")
    ax = axes[1]
    rows = sorted([r for r in F34["fig3_bottom"] if r["w"] == 2.0], key=lambda r: r["delta_U"])
    dus_b = [r["delta_U"] for r in rows]; cmap_b = _delta_u_colors(dus_b)
    g = [r[version][MODE]["gamma_pi_wrapped"] for r in rows]
    ax.plot(dus_b, g, "-", color="gray", alpha=0.5, lw=1.5, zorder=1)
    ax.scatter(dus_b, g, c=[cmap_b[round(d, 6)] for d in dus_b], s=40, zorder=2)
    ax.axhline(0, color="gray", ls="--", alpha=0.5); _symlog_deltaU_axis(ax, dus_b)
    ax.set_xlabel(r"Hubbard Interaction $\Delta U$", fontsize=10); ax.tick_params(labelsize=9)
    ax.set_ylabel(r"Berry Phase $\gamma$ ($\pi$)"); ax.set_ylim(-0.55, 0.45)

# ------------------------------------------------------------------ Fig. 4
def draw_fig4(axes, version):
    rows = sorted(F34["fig4"], key=lambda r: r["delta_U"]); dus = [r["delta_U"] for r in rows]; cmap = _delta_u_colors(dus)
    ax = axes[0]
    for r in rows:
        ax.plot(range(1, 7), r[version]["pol"], "o-", color=cmap[round(r["delta_U"], 6)], ms=4, lw=1, label=fr"$\Delta U={r['delta_U']:g}$")
    ax.axhline(0, color="gray", ls="--", alpha=0.6); ax.set_xticks(range(1, 7)); ax.grid(True, alpha=0.3)
    ax.set_xlabel("Unit Cell Index $j$"); ax.set_ylabel(r"Polarization $\langle n_{A,j}\rangle-\langle n_{B,j}\rangle$")
    ax.legend(fontsize=7.5, ncol=2, loc="upper right"); ax.set_ylim(-1.4, 1.4)
    ax = axes[1]
    v = [r[version]["pol"][0] for r in rows]
    ax.plot(dus, v, "-", color="gray", alpha=0.5, lw=1.5, zorder=1)
    ax.scatter(dus, v, c=[cmap[round(d, 6)] for d in dus], s=40, zorder=2)
    _symlog_deltaU_axis(ax, dus); ax.set_xlabel(r"Hubbard Interaction $\Delta U$", fontsize=10); ax.tick_params(labelsize=9)
    ax.set_ylabel(r"Polarization $\langle n_{A,1}\rangle-\langle n_{B,1}\rangle$"); ax.set_ylim(0.68, 1.27)

# ------------------------------------------------------------------ per-figure 2x2 compares
specs = [
    ("paper_fig2_compare", draw_fig2, "論文 Fig. 2：Trotter 收斂（上 U=0 對初態、下 U=1 相鄰 L）", ["上：U = 0", "下：U = 1"]),
    ("paper_fig3_compare", draw_fig3, "論文 Fig. 3：多體 Berry 相位（上 γ 對 w、下 w=2 時 γ 對 ΔU）", ["上：γ 對 w", "下：w = 2"]),
    ("paper_fig4_compare", draw_fig4, "論文 Fig. 4：子格點極化（上 剖面、下 第一晶胞對 ΔU）", ["上：剖面", "下：cell 1"]),
]
for name, draw, suptitle, _ in specs:
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.2))
    draw(axes[:, 0], "before"); draw(axes[:, 1], "after")
    title(axes[0, 0], "修正前（論文所用程式）", C_BEFORE); title(axes[0, 1], "修正後（三個缺陷皆修正）", C_AFTER)
    fig.suptitle(suptitle, x=0.01, ha="left", fontsize=12); fig.tight_layout()
    fig.savefig(FIG / f"{name}.svg", bbox_inches="tight")
    if os.environ.get("PNG_DIR"): fig.savefig(pathlib.Path(os.environ["PNG_DIR"]) / f"{name}.png", dpi=110, bbox_inches="tight")
    plt.close(fig); print("wrote", name)

# ------------------------------------------------------------------ one tall poster with everything
fig, axes = plt.subplots(6, 2, figsize=(12, 20))
draw_fig2(axes[0:2, 0], "before"); draw_fig2(axes[0:2, 1], "after")
draw_fig3(axes[2:4, 0], "before"); draw_fig3(axes[2:4, 1], "after")
draw_fig4(axes[4:6, 0], "before"); draw_fig4(axes[4:6, 1], "after")
for row, lab in ((0, "Fig. 2 上（U = 0）"), (1, "Fig. 2 下（U = 1）"), (2, "Fig. 3 上"), (3, "Fig. 3 下（w = 2）"), (4, "Fig. 4 上"), (5, "Fig. 4 下（cell 1）")):
    axes[row, 0].annotate(lab, (-0.22, 0.5), xycoords="axes fraction", rotation=90, va="center", ha="center", fontsize=11, weight="bold")
title(axes[0, 0], "修正前（論文所用程式）", C_BEFORE); title(axes[0, 1], "修正後（三個缺陷皆修正）", C_AFTER)
fig.suptitle("論文全部圖形：修正前 | 修正後（N=6, T=1, L=40；Fig. 2 依論文掃 T、L）", x=0.01, ha="left", fontsize=13)
fig.tight_layout(rect=(0.02, 0, 1, 0.985)); fig.savefig(FIG / "paper_all_compare.svg", bbox_inches="tight")
if os.environ.get("PNG_DIR"): fig.savefig(pathlib.Path(os.environ["PNG_DIR"]) / "paper_all_compare.png", dpi=80, bbox_inches="tight")
plt.close(fig); print("wrote paper_all_compare")

# ------------------------------------------------------------------ numbers for the text
L = F2["meta"]["L"]
worst = 0.0
for tag in ("U0", "U1"):
    for T in F2["meta"]["T"]:
        for ref in ("ground", "adjacent"):
            a = np.array(F2[tag]["before"][str(T)][ref]); b = np.array(F2[tag]["after"][str(T)][ref])
            worst = max(worst, float(np.max(np.abs(a - b))))
print(f"Fig 2: max |before - after| over all points = {worst:.2e}")
