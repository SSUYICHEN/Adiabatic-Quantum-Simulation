"""Generate the figures for docs/defect1/defect1.typ into docs/defect1/figures/.

Run from the repo root:
    uv run --no-sync --with pylatexenc python docs/defect1/make_figures.py
"""
import os, pathlib, shutil, contextlib
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import openfermion as of
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
import matplotlib.font_manager as fm
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector, Operator
from qiskit.circuit.library import UnitaryGate

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent.parent
FIG = HERE / "figures"; FIG.mkdir(exist_ok=True)

_cjk = [f.name for f in fm.fontManager.ttflist if "CJK" in f.name and "Sans" in f.name]
matplotlib.rcParams["font.family"] = ([_cjk[0]] if _cjk else []) + ["DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
BEFORE, AFTER, REF, INK, GRAY = "#eb6834", "#2a78d6", "#1baf7a", "#1a1a1a", "#8a8a8a"

def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e5e5e5", linewidth=0.6); ax.set_axisbelow(True); ax.tick_params(length=0)

def save(fig, name):
    fig.savefig(FIG / f"{name}.svg", bbox_inches="tight")
    plt.close(fig); print("wrote", name)

# ------------------------------------------------------------------ the two implementations
def legacy_givens(theta, phi):
    gv = QuantumCircuit(2, name="Givens")
    gv.rz(phi, 1); gv.rz(-phi, 0); gv.cx(0, 1); gv.cry(-2.0 * theta, 1, 0); gv.cx(0, 1)
    return gv.to_instruction()

def fixed_givens(theta, phi):
    c, s = np.cos(theta), np.sin(theta); e = np.exp(1j * phi)
    U = np.zeros((4, 4), dtype=complex)
    U[0, 0] = 1.0; U[1, 1] = c; U[2, 2] = c * e; U[1, 2] = s; U[2, 1] = -s * e; U[3, 3] = e
    return UnitaryGate(U, label="Givens")

def as_matrix(instr):
    qc = QuantumCircuit(2); qc.append(instr, [0, 1]); return Operator(qc).data

# ================================================================== 1. workflow diagram
fig, ax = plt.subplots(figsize=(11, 2.6)); ax.axis("off")
boxes = [
    (0.02, "軌域係數矩陣 Q\n(ν × N，Q Q† = I)", "#eef3f8"),
    (0.27, "OpenFermion\nslater_determinant_\npreparation_circuit(Q)", "#fff8e6"),
    (0.52, "旋轉序列\n{(j, k, θ, φ), …}\n(分成可平行的層)", "#eef3f8"),
    (0.77, "電路：X 閘放粒子\n+ 逐一套用 G(θ, φ)\n→ |Ψ⟩", "#eef8f0"),
]
for x, txt, col in boxes:
    ax.add_patch(FancyBboxPatch((x, 0.2), 0.2, 0.6, boxstyle="round,pad=0.01,rounding_size=0.02",
                                fc=col, ec="#334155", lw=1.2, transform=ax.transAxes))
    ax.text(x + 0.1, 0.5, txt, ha="center", va="center", fontsize=10.5, transform=ax.transAxes)
for x in (0.22, 0.47, 0.72):
    ax.add_patch(FancyArrowPatch((x + 0.005, 0.5), (x + 0.045, 0.5), transform=ax.transAxes,
                                 arrowstyle="-|>", mutation_scale=18, lw=1.5, color="#334155"))
ax.text(0.87, 0.06, "缺陷所在：這一步把 (θ, φ) 變成閘的方式", ha="center", fontsize=9.5, color=BEFORE, transform=ax.transAxes)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
save(fig, "fig_workflow")

# ================================================================== 2. circuit drawings
theta, phi = 0.7, 0.9
qc_leg = QuantumCircuit(2); qc_leg.rz(phi, 1); qc_leg.rz(-phi, 0); qc_leg.cx(0, 1); qc_leg.cry(-2.0 * theta, 1, 0); qc_leg.cx(0, 1)
fig = qc_leg.draw("mpl", style={"name": "bw"}, fold=-1, scale=1.1)
fig.suptitle("修正前：legacy_givens(θ, φ) 的閘分解（cry 角度為 −2θ）", fontsize=11, x=0.02, ha="left", color=BEFORE)
save(fig, "fig_circuit_legacy")

qc_fix = QuantumCircuit(2); qc_fix.append(fixed_givens(theta, phi), [0, 1])
fig = qc_fix.draw("mpl", style={"name": "bw"}, fold=-1, scale=1.1)
fig.suptitle("修正後：直接以 4×4 UnitaryGate 實作 G(θ, φ)", fontsize=11, x=0.02, ha="left", color=AFTER)
save(fig, "fig_circuit_fixed")

# ================================================================== 3. matrix comparison
G = as_matrix(fixed_givens(theta, phi)); Ulegacy = as_matrix(legacy_givens(theta, phi))
labels = ["|00⟩", "|01⟩", "|10⟩", "|11⟩"]
def draw_matrix(ax, M, title, color):
    ax.imshow(np.abs(M), cmap="Blues", vmin=0, vmax=1.15, alpha=0.55)
    for i in range(4):
        for j in range(4):
            z = M[i, j]
            if abs(z) < 1e-9: txt, c = "0", "#b0b0b0"
            elif abs(z.imag) < 1e-9: txt, c = f"{z.real:+.3f}", INK
            else: txt, c = f"{z.real:+.3f}\n{z.imag:+.3f}i", INK
            ax.text(j, i, txt, ha="center", va="center", fontsize=8.5, color=c)
    ax.set_xticks(range(4)); ax.set_yticks(range(4)); ax.set_xticklabels(labels); ax.set_yticklabels(labels)
    ax.tick_params(length=0); ax.set_title(title, fontsize=10.5, loc="left", color=color)
    for s in ax.spines.values(): s.set_color("#cccccc")
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
draw_matrix(axes[0], G, "OpenFermion 定義的 G(θ, φ)（修正後）", AFTER)
draw_matrix(axes[1], Ulegacy, "legacy_givens(θ, φ) 實際的矩陣（修正前）", BEFORE)
draw_matrix(axes[2], Ulegacy - G, "差異 legacy − G", INK)
fig.suptitle(f"θ = {theta}, φ = {phi}：基底順序 |q1 q0⟩，列 = 輸出、欄 = 輸入", fontsize=10.5, x=0.02, ha="left")
fig.tight_layout(); save(fig, "fig_matrix_compare")

# phi = 0 special case: pure sign flip of theta
G0 = as_matrix(fixed_givens(theta, 0.0)); U0 = as_matrix(legacy_givens(theta, 0.0)); Gm = as_matrix(fixed_givens(-theta, 0.0))
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
draw_matrix(axes[0], G0, "G(θ, 0)（應得）", AFTER)
draw_matrix(axes[1], U0, "legacy(θ, 0)", BEFORE)
draw_matrix(axes[2], Gm, "G(−θ, 0)：與 legacy 完全相同", BEFORE)
fig.suptitle(f"實數躍遷 φ = 0 時，legacy 恰等於把 θ 反號的 G（θ = {theta}）", fontsize=10.5, x=0.02, ha="left")
fig.tight_layout(); save(fig, "fig_matrix_phi0")
print("max|legacy(θ,0) − G(−θ,0)| =", np.max(np.abs(U0 - Gm)))

# ================================================================== 4. hydrogen-molecule-like example
# two atoms (A, B), two electrons, tight-binding H = -t (a†b + b†a) per spin; modes: 0=A↑ 1=B↑ 2=A↓ 3=B↓
t = 1.0
Qup = np.array([[1.0, 1.0]]) / np.sqrt(2)     # bonding orbital per spin
def prepare_h2(givens_fn):
    qc = QuantumCircuit(4); qc.x(0); qc.x(2)
    for ops in of.circuits.slater_determinant_preparation_circuit(Qup):
        for j, k, th, ph in ops:
            qc.append(givens_fn(th, ph), [j, k]); qc.append(givens_fn(th, ph), [j + 2, k + 2])
    return qc, Statevector(qc).data
qc_h2, psi_leg = prepare_h2(legacy_givens); _, psi_fix = prepare_h2(fixed_givens)
# reference: bonding⊗bonding = (|A↑⟩+|B↑⟩)(|A↓⟩+|B↓⟩)/2 -> occupation kets (little-endian index bits q0..q3)
ref = np.zeros(16, dtype=complex)
for up in (0, 1):
    for dn in (2, 3):
        ref[(1 << up) | (1 << dn)] = 0.5
# energy: H = -t sum_spin (a†b + h.c.) as 16x16 via openfermion, little-endian
op = of.FermionOperator()
for (i, j) in ((0, 1), (2, 3)):
    op += of.FermionOperator(f"{i}^ {j}", -t) + of.FermionOperator(f"{j}^ {i}", -t)
Hbig = of.get_sparse_operator(op, n_qubits=4).toarray()
perm = np.zeros(16, dtype=int)
for idx in range(16):
    for p in range(4): perm[idx] |= ((idx >> p) & 1) << (3 - p)
H = Hbig[np.ix_(perm, perm)]
E = lambda v: float(np.vdot(v, H @ v).real)
E_leg, E_fix, E_ref = E(psi_leg), E(psi_fix), E(ref)
fid = lambda a, b: float(abs(np.vdot(a, b)) ** 2)
print(f"H2-like: E_legacy={E_leg:+.4f} E_fixed={E_fix:+.4f} E_ref={E_ref:+.4f}; fid(legacy,ref)={fid(psi_leg, ref):.4f} fid(fixed,ref)={fid(psi_fix, ref):.4f}")

fig = qc_h2.draw("mpl", style={"name": "bw"}, fold=-1, scale=1.0)
fig.suptitle("氫分子模型的態製備電路：X 閘放入兩個電子，再對每個自旋套用一個 Givens 旋轉", fontsize=10.5, x=0.02, ha="left")
save(fig, "fig_circuit_h2")

kets = [(1 << up) | (1 << dn) for up in (0, 1) for dn in (2, 3)]
names = ["|A↑ A↓⟩", "|A↑ B↓⟩", "|B↑ A↓⟩", "|B↑ B↓⟩"]
# align global phase to reference
def align(v): ph = np.vdot(ref, v); return v * np.exp(-1j * np.angle(ph)) if abs(ph) > 1e-12 else v
vals = {"參考（鍵結 ⊗ 鍵結）": (ref, REF), "修正前 legacy": (align(psi_leg), BEFORE), "修正後 G(θ, φ)": (align(psi_fix), AFTER)}
fig, ax = plt.subplots(figsize=(9, 3.8)); x = np.arange(4); w = 0.26
for i, (lab, (v, col)) in enumerate(vals.items()):
    amps = [v[k].real for k in kets]
    bars = ax.bar(x + (i - 1) * (w + 0.01), amps, w, color=col, label=lab)
    for b, a in zip(bars, amps):
        ax.annotate(f"{a:+.2f}", (b.get_x() + b.get_width() / 2, a), xytext=(0, 3 if a >= 0 else -11), textcoords="offset points", ha="center", fontsize=8)
ax.axhline(0, color="#999", lw=0.8); ax.set_xticks(x); ax.set_xticklabels(names)
ax.set_ylabel("振幅（實部）"); ax.set_title(f"氫分子模型的製備態振幅：⟨H⟩ 修正前 {E_leg:+.2f}t，修正後 {E_fix:+.2f}t（基態 {E_ref:+.2f}t）", fontsize=10.5, loc="left")
ax.legend(frameon=False, fontsize=9, loc="lower left"); ax.margins(y=0.25); style(ax)
save(fig, "fig_h2_amplitudes")

# density is identical -> figure showing why occupation numbers cannot tell
dens = lambda v: [float(sum(abs(v[i]) ** 2 for i in range(16) if (i >> q) & 1)) for q in range(4)]
fig, ax = plt.subplots(figsize=(6.5, 3.2)); x = np.arange(4)
ax.bar(x - 0.2, dens(psi_leg), 0.38, color=BEFORE, label="修正前"); ax.bar(x + 0.2, dens(psi_fix), 0.38, color=AFTER, label="修正後")
ax.set_xticks(x); ax.set_xticklabels(["A↑", "B↑", "A↓", "B↓"]); ax.set_ylabel("⟨n_q⟩"); ax.set_ylim(0, 0.8)
ax.set_title("兩個態的佔據數完全相同（各 0.5）：只看密度看不出錯", fontsize=10.5, loc="left"); ax.legend(frameon=False, fontsize=9); style(ax)
save(fig, "fig_h2_density")

# ================================================================== 5. SSH spectrum: lowest vs highest filling
import sys; sys.path.insert(0, str(REPO / "src"))
import aqs.core as core
from aqs.core import SSHHModel, build_annealing_circuit
from aqs.backends import get_backend
bk = get_backend("qiskit")
@contextlib.contextmanager
def patched(mod, name, repl):
    old = getattr(mod, name); setattr(mod, name, repl)
    try: yield
    finally: setattr(mod, name, old)
m = SSHHModel.from_total_electrons(3, 0.5, 1.5, 8, PBC=False)
Qu, Qd, nu, nd = m.slater_Q_matrices(); L = m.L
ev = np.sort(np.linalg.eigvalsh(m.H[:L, :L]))
def prep_energy(fn):
    with patched(core, "_givens_instruction", fn):
        sv = bk.statevector(build_annealing_circuit(m, Qu, Qd, 0, 0, T_A=0, steps=0, ramp_U=False))
    opf = of.FermionOperator()
    for i in range(2 * L):
        for j in range(2 * L):
            if abs(m.H[i, j]) > 1e-12: opf += of.FermionOperator(f"{i}^ {j}", complex(m.H[i, j]))
    Hs = of.get_sparse_operator(opf, n_qubits=2 * L).tocsr()
    n = 2 * L; idx = np.arange(1 << n); pm = np.zeros_like(idx)
    for p in range(n): pm |= ((idx >> p) & 1) << (n - 1 - p)
    Hs = Hs[pm][:, pm]
    return float(np.vdot(sv, Hs @ sv).real)
E_before, E_after = prep_energy(legacy_givens), prep_energy(core._givens_instruction)
fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
for ax, title, col, fill_lowest, Ep in ((axes[1], "修正後：填最低 4 個軌域", AFTER, True, E_after), (axes[0], "修正前：實際填的是最高 4 個軌域", BEFORE, False, E_before)):
    for k, e in enumerate(ev):
        ax.hlines(e, 0.2, 0.8, color="#888", lw=1.4)
    occ = ev[:nu] if fill_lowest else ev[::-1][:nu]
    for e in occ:
        ax.plot([0.35, 0.65], [e, e], "o", color=col, ms=9, mec="white")
    ax.text(0.5, ev[-1] + 0.5, f"每個自旋填 {nu} 個電子\n⟨H⟩ = {Ep:+.3f}", ha="center", fontsize=9.5, color=col)
    ax.set_xlim(0, 1); ax.set_xticks([]); ax.set_title(title, fontsize=10.5, loc="left", color=col)
    ax.axhline(0, color="#ccc", lw=0.8, ls="--"); ax.spines[["top", "right", "bottom"]].set_visible(False)
axes[0].set_ylabel("單粒子能階 ε（SSH，N=3，v=0.5，w=1.5，OBC）"); axes[0].set_ylim(ev[0] - 0.8, ev[-1] + 1.6)
save(fig, "fig_spectrum")
print(f"SSH N=3 OBC: E_before={E_before:+.6f} E_after={E_after:+.6f} lowest-fill={ev[:nu].sum()*2:+.6f}")

# ================================================================== 6. sublattice gauge b -> -b sketch
fig, axes = plt.subplots(2, 1, figsize=(10, 3.6))
for ax, title, sgn, col in ((axes[0], "原模型 H(v, w, U)：舊碼製備的是它的最高本徵態", "+", INK),
                            (axes[1], "作 b → −b 之後：−H(v, w, −U)，最高本徵態變成「U 反號模型」的基態", "−", BEFORE)):
    ax.axis("off"); ax.set_xlim(-0.5, 6.5); ax.set_ylim(-1.35, 1.2)
    for j in range(3):
        xa, xb = 2 * j, 2 * j + 1
        ax.add_patch(Circle((xa, 0), 0.22, fc="#2a78d6", ec="none")); ax.text(xa, 0, "A", ha="center", va="center", color="white", fontsize=10, weight="bold")
        ax.add_patch(Circle((xb, 0), 0.22, fc="#eb6834" if sgn == "−" else "#666", ec="none")); ax.text(xb, 0, "B", ha="center", va="center", color="white", fontsize=10, weight="bold")
        ax.plot([xa + 0.22, xb - 0.22], [0, 0], color="#333", lw=2); ax.text((xa + xb) / 2, 0.32, f"{sgn if sgn=='−' else ''}v", ha="center", fontsize=10)
        if j < 2:
            ax.plot([xb + 0.22, xb + 1 - 0.22], [0, 0], color="#333", lw=2); ax.text(xb + 0.5, 0.32, f"{sgn if sgn=='−' else ''}w", ha="center", fontsize=10)
        ax.text(xa, -0.55, f"{'−' if sgn=='−' else ''}U_A n↑n↓" if sgn == "−" else "U_A n↑n↓", ha="center", fontsize=8.5, color="#444")
        ax.text(xb, -0.55, f"{'−' if sgn=='−' else ''}U_B n↑n↓" if sgn == "−" else "U_B n↑n↓", ha="center", fontsize=8.5, color="#444")
    ax.text(-0.4, 0.95, title, fontsize=10.5, color=col)
axes[1].text(3.0, -1.15, "躍遷全部反號、Hubbard 項不變 ⇒ 整體乘 −1 後等於 U → −U 的模型", ha="center", fontsize=9.5, color=BEFORE)
save(fig, "fig_gauge")

# ================================================================== 7. copy paper-impact figures
for name in ("A_polarization", "paper_fig3_bottom", "paper_fig4_top", "paper_fig4_bottom", "paper_fig3_top"):
    src = REPO / "docs" / "defect_report" / "figures" / f"{name}.svg"
    shutil.copy(src, FIG / f"{name}.svg"); print("copied", name)
for name in ("_page_9_Figure_0.jpeg", "_page_9_Figure_10.jpeg"):
    shutil.copy(REPO / ".claude" / "specs" / "paper_theory" / name, FIG / name); print("copied", name)
