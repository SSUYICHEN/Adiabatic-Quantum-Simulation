import json, pathlib, sys
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
_cjk = [f.name for f in fm.fontManager.ttflist if "CJK" in f.name and "Sans" in f.name]
matplotlib.rcParams["font.family"] = ([_cjk[0]] if _cjk else []) + ["DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

S = pathlib.Path(__file__).parent
R = json.load(open(S / "paper_impact.json"))
FIG = S / "figures"
ORANGE = ["#f8c9b2", "#f4a883", "#ef8757", "#eb6834", "#c94e1f", "#a03a12", "#772a0b", "#4d1a05"]
BLUE   = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]
MODE = "up_spin"   # the code's default, i.e. what the paper's figures were produced with

def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color="#e5e5e5", linewidth=0.6); ax.set_axisbelow(True); ax.tick_params(length=0)

# ---------------- Fig 3 top: gamma vs w, one curve per dU, before | after
dUs = sorted({r["delta_U"] for r in R["fig3_top"]})
fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
for ax, tag, ramp, title in ((axes[0], "before", ORANGE, "修正前（論文 Fig. 3 上圖所用程式）"), (axes[1], "after", BLUE, "修正後")):
    for i, dU in enumerate(dUs):
        pts = sorted([(r["w"], r[tag][MODE]["gamma_pi_wrapped"]) for r in R["fig3_top"] if r["delta_U"] == dU])
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "-o", ms=4, lw=1.4, color=ramp[i], label=f"ΔU = {dU:g}")
    ax.axhline(0, color="#999", lw=0.8, ls="--"); ax.axhline(1, color="#999", lw=0.8, ls="--")
    ax.set_xlabel("inter-cell hopping w"); ax.set_title(title, loc="left", fontsize=11); style(ax)
axes[0].set_ylabel("Berry phase γ (π)"); axes[1].legend(frameon=False, fontsize=8, ncol=2, loc="center right")
fig.suptitle("論文 Fig. 3 上圖重算（N=6, v=1, U_A=0.01, T=1, L=40）", x=0.01, ha="left", fontsize=11)
fig.tight_layout(); fig.savefig(FIG / "paper_fig3_top.svg"); 

# zoom on topological side
fig, ax = plt.subplots(figsize=(7, 3.6))
for i, dU in enumerate(dUs):
    for tag, ramp, ls in (("before", ORANGE, "--"), ("after", BLUE, "-")):
        pts = sorted([(r["w"], r[tag][MODE]["gamma_pi_wrapped"]) for r in R["fig3_top"] if r["delta_U"] == dU and r["w"] > 1])
        ax.plot([p[0] for p in pts], [p[1] for p in pts], ls, marker="o", ms=3.5, lw=1.2, color=ramp[i],
                label=(f"ΔU = {dU:g}" if tag == "after" else None))
ax.axhline(0, color="#999", lw=0.8)
ax.set_xlabel("w（拓樸相，w > 1）"); ax.set_ylabel("γ (π)"); ax.set_title("放大：虛線 = 修正前（橘），實線 = 修正後（藍）", loc="left", fontsize=11); style(ax)
ax.legend(frameon=False, fontsize=8, ncol=2)
fig.tight_layout(); fig.savefig(FIG / "paper_fig3_top_zoom.svg"); 

# ---------------- Fig 3 bottom: gamma vs dU at w = 2.0 (and 1.5), log-x
fig, axes = plt.subplots(1, 2, figsize=(12, 3.8), sharey=True)
for ax, w in zip(axes, (1.5, 2.0)):
    rows = sorted([r for r in R["fig3_bottom"] if r["w"] == w], key=lambda r: r["delta_U"])
    x = [max(r["delta_U"], 1e-4) for r in rows]
    ax.plot(x, [r["before"][MODE]["gamma_pi_wrapped"] for r in rows], "--o", color="#eb6834", ms=5, label="修正前（論文所用）")
    ax.plot(x, [r["after"][MODE]["gamma_pi_wrapped"] for r in rows], "-o", color="#2a78d6", ms=5, label="修正後")
    ax.plot(x, [r["after_negU"][MODE]["gamma_pi_wrapped"] for r in rows], "x", color="#1c5cab", ms=9, mew=1.5, label="修正後，但 U_A, U_B → −U_A, −U_B")
    ax.set_xscale("log"); ax.axhline(0, color="#999", lw=0.8)
    ax.set_xlabel("ΔU（0 畫在 10⁻⁴）"); ax.set_title(f"w = {w}", loc="left", fontsize=11); style(ax)
axes[0].set_ylabel("Berry phase γ (π)"); axes[0].legend(frameon=False, fontsize=8)
fig.suptitle("論文 Fig. 3 下圖重算：γ 對 ΔU（N=6, v=1, U_A=0.01）", x=0.01, ha="left", fontsize=11)
fig.tight_layout(); fig.savefig(FIG / "paper_fig3_bottom.svg"); 

# ---------------- Fig 4: polarization profile, before | after
rows4 = sorted(R["fig4"], key=lambda r: r["delta_U"])
fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
for ax, tag, ramp, title in ((axes[0], "before", ORANGE, "修正前（論文 Fig. 4 上圖所用程式）"), (axes[1], "after", BLUE, "修正後")):
    for i, r in enumerate(rows4):
        ax.plot(range(1, 7), r[tag]["pol"], "-o", ms=4, lw=1.4, color=ramp[i + 1], label=f"ΔU = {r['delta_U']:g}")
    ax.axhline(0, color="#999", lw=0.8, ls="--"); ax.set_xlabel("unit cell j"); ax.set_title(title, loc="left", fontsize=11); style(ax)
axes[0].set_ylabel("⟨n_A,j⟩ − ⟨n_B,j⟩"); axes[1].legend(frameon=False, fontsize=8, ncol=2, loc="lower left")
fig.suptitle("論文 Fig. 4 上圖重算（N=6, v=0.1, w=1, U_A=0.01, OBC, 14 電子）", x=0.01, ha="left", fontsize=11)
fig.tight_layout(); fig.savefig(FIG / "paper_fig4_top.svg"); 

fig, ax = plt.subplots(figsize=(6.5, 3.6))
x = [max(r["delta_U"], 1e-3) for r in rows4]
ax.plot(x, [r["before"]["pol"][0] for r in rows4], "--o", color="#eb6834", ms=5, label="修正前（論文所用）")
ax.plot(x, [r["after"]["pol"][0] for r in rows4], "-o", color="#2a78d6", ms=5, label="修正後")
ax.plot(x, [r["after_negU"]["pol"][0] for r in rows4], "x", color="#1c5cab", ms=9, mew=1.5, label="修正後，U → −U")
ax.set_xscale("log"); ax.set_xlabel("ΔU（0 畫在 10⁻³）"); ax.set_ylabel("cell 1 polarization")
ax.set_title("論文 Fig. 4 下圖重算：第一晶胞極化對 ΔU", loc="left", fontsize=11); style(ax); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "paper_fig4_bottom.svg"); 

# ---------------- text summary
print("=== Fig 3 top: gamma/pi at w=1.5 (topological) and w=0.5 (trivial)")
for dU in dUs:
    r15 = [r for r in R["fig3_top"] if r["delta_U"] == dU and r["w"] == 1.5][0]
    r05 = [r for r in R["fig3_top"] if r["delta_U"] == dU and r["w"] == 0.5][0]
    print(f"dU={dU:<7} w=1.5 before {r15['before'][MODE]['gamma_pi_wrapped']:+.4f} after {r15['after'][MODE]['gamma_pi_wrapped']:+.4f} | "
          f"w=0.5 before {r05['before'][MODE]['gamma_pi_wrapped']:+.4f} after {r05['after'][MODE]['gamma_pi_wrapped']:+.4f}")
print("=== Fig 3 bottom (w=2.0): before / after / after(-U)")
for r in sorted([r for r in R["fig3_bottom"] if r["w"] == 2.0], key=lambda r: r["delta_U"]):
    print(f"dU={r['delta_U']:<7} {r['before'][MODE]['gamma_pi_wrapped']:+.4f} {r['after'][MODE]['gamma_pi_wrapped']:+.4f} {r['after_negU'][MODE]['gamma_pi_wrapped']:+.4f}")
print("=== Fig 4: pol profile before / after ; after(-U)")
for r in rows4:
    b = np.round(r["before"]["pol"], 3); a = np.round(r["after"]["pol"], 3); n = np.round(r["after_negU"]["pol"], 3)
    print(f"dU={r['delta_U']:<5} before {b}\n          after  {a}\n          after(-U) {n}  max|before-after(-U)|={np.max(np.abs(np.array(r['before']['pol'])-np.array(r['after_negU']['pol']))):.1e}")
print("=== total-mode check (paper Eq.46) at w=1.5:")
for dU in (0.0, 0.1, 0.3):
    r15 = [r for r in R["fig3_top"] if r["delta_U"] == dU and r["w"] == 1.5][0]
    print(f"dU={dU} up_spin after {r15['after']['up_spin']['gamma_pi_wrapped']:+.4f}  total after {r15['after']['total']['gamma_pi_wrapped']:+.4f}  total before {r15['before']['total']['gamma_pi_wrapped']:+.4f}")
