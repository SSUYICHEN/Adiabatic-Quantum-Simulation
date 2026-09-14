// =====================================================================
//  AQS - Adiabatic Quantum Simulation
//  程式碼架構與物理演算法技術文件
//  以 typst compile docs/codebase_architecture.typ 編譯
// =====================================================================

#set document(
  title: "AQS 絕熱量子模擬框架：程式碼架構與物理演算法技術文件",
  author: "Adiabatic-Quantum-Simulation",
)

#set page(
  paper: "a4",
  margin: (x: 2.2cm, y: 2.4cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 8pt, fill: rgb("#666666"))
      AQS 技術文件 · 絕熱量子模擬框架
      #h(1fr)
      v2.0.0
      #line(length: 100%, stroke: 0.4pt + rgb("#cccccc"))
    ]
  },
)

#set text(
  font: ("New Computer Modern", "Noto Serif CJK TC"),
  size: 10.5pt,
  lang: "zh",
  region: "tw",
)

#set par(justify: true, leading: 0.78em, first-line-indent: 0em)
#set heading(numbering: "1.1")

#show heading: it => {
  set text(font: ("Noto Sans", "Noto Sans CJK TC"))
  block(above: 1.4em, below: 0.8em, it)
}
#show heading.where(level: 1): it => {
  set text(size: 17pt, fill: rgb("#0b3d61"))
  block(above: 1.8em, below: 1.0em)[
    #it
    #v(-0.5em)
    #line(length: 100%, stroke: 1.2pt + rgb("#0b3d61"))
  ]
}
#show heading.where(level: 2): set text(size: 13pt, fill: rgb("#14507d"))
#show heading.where(level: 3): set text(size: 11.5pt, fill: rgb("#1f6392"))

#show raw: set text(font: ("DejaVu Sans Mono", "Noto Sans Mono CJK TC"), size: 8.8pt)
#show raw.where(block: true): it => block(
  width: 100%,
  fill: rgb("#f6f8fa"),
  stroke: (left: 2.5pt + rgb("#0b3d61"), rest: 0.4pt + rgb("#dfe3e8")),
  radius: 2pt,
  inset: (x: 9pt, y: 8pt),
  it,
)
#show link: set text(fill: rgb("#14507d"))

// ------------------------------------------------- 共用小工具
#let note(title, body) = block(
  width: 100%, fill: rgb("#fff8e6"),
  stroke: (left: 2.5pt + rgb("#d9a300"), rest: 0.4pt + rgb("#e8dcb8")),
  radius: 2pt, inset: (x: 10pt, y: 9pt), above: 1em, below: 1em,
)[*#title* \ #body]

#let danger(title, body) = block(
  width: 100%, fill: rgb("#fdeff0"),
  stroke: (left: 2.5pt + rgb("#c0392b"), rest: 0.4pt + rgb("#eccfcb")),
  radius: 2pt, inset: (x: 10pt, y: 9pt), above: 1em, below: 1em,
)[*#title* \ #body]

#let ok(title, body) = block(
  width: 100%, fill: rgb("#eef8f0"),
  stroke: (left: 2.5pt + rgb("#1e8449"), rest: 0.4pt + rgb("#cfe6d5")),
  radius: 2pt, inset: (x: 10pt, y: 9pt), above: 1em, below: 1em,
)[*#title* \ #body]

#let tbl(cols, aligns, ..cells) = table(
  columns: cols,
  align: aligns,
  stroke: 0.4pt + rgb("#c9d1d9"),
  fill: (_, y) => if y == 0 { rgb("#eef2f6") } else { white },
  inset: (x: 7pt, y: 5.5pt),
  ..cells,
)

// =====================================================================
#align(center)[
  #v(1.5cm)
  #text(size: 24pt, font: ("Noto Sans", "Noto Sans CJK TC"), fill: rgb("#0b3d61"))[
    *AQS 絕熱量子模擬框架*
  ]
  #v(0.3cm)
  #text(size: 14pt, fill: rgb("#14507d"))[
    程式碼架構、物理演算法與缺陷修正技術文件
  ]
  #v(0.6cm)
  #text(size: 10.5pt)[
    自旋 SSH--Hubbard（$4N$ 量子位元）至\
    延伸無自旋 SSH--Ising（$2N$ 量子位元）之演進
  ]
  #v(0.5cm)
  #text(size: 9.5pt, fill: rgb("#666666"))[
    套件版本 `aqs 2.0.0` · Python $>=$ 3.10 · qiskit / CUDA-Q 雙後端
  ]
  #v(1.2cm)
]

#outline(depth: 3, indent: auto)
#pagebreak()

// =====================================================================
= 執行摘要

本專案 `aqs` 是一套以*絕熱量子演化*（Adiabatic Quantum Simulation, AQS）量測一維拓樸系統
不變量的模擬工具組，同時提供命令列介面與 Python API，並可在 CPU（qiskit）或
NVIDIA GPU（CUDA-Q）後端執行。整體程式碼約 2,100 行（含進入點共十個檔案），分為下列九個主要模組。

== 兩套模型的並存

本框架目前同時維護兩套*物理上不同*的模型，兩者共用閘極原語與可觀測量登錄表，
但電路拓樸完全不同：

#tbl(
  (auto, 1fr, 1fr),
  (left, left, left),
  [*項目*], [*內建自旋 SSH--Hubbard*], [*延伸無自旋 SSH--Ising*],
  [模組], [`core.py`], [`spinless.py`],
  [量子位元數], [$4N$（自旋區塊映射）], [$2N$（每格點一個位元）],
  [映射], [$[0,L)$ 為上自旋、$[L,2L)$ 為下自旋], [$q_(2j) = A_j$，$q_(2j+1) = B_j$],
  [交互作用], [在位 Hubbard $U n_arrow.t n_arrow.b$], [最近鄰 $V_v n_A n_B + V_w n_B n_A$],
  [CP 閘作用], [跨兩個自旋區塊], [沿鏈方向的相鄰格點],
  [斜坡參數], [$s = (l - 0.5) \/ L$], [$lambda_l = (2l-1) \/ (2L)$],
)

其中「延伸無自旋 SSH--Ising」指的是無自旋 SSH 鏈加上密度--密度（Ising 型）
最近鄰交互作用；由於 $n_i n_j$ 在計算基底下為對角算符，其量子閘實作即為受控相位閘，
故稱 Ising 型。

== 本次工作的三項核心成果

+ *新增 $2N$ 無自旋模型*：完整實作 `SpinlessSSHHSim`，涵蓋精確 Slater 行列式態製備、
  Trotter 化絕熱演化、PBC/OBC、複數躍遷振幅、以及精確／取樣兩種量測協定。

+ *修正三個靜默的正確性缺陷*：CUDA-Q 位元序反轉、OpenFermion 大端序與
  qiskit 小端序不一致、以及 Givens 旋轉符號錯誤導致製備出*能量最高*的行列式態。
  三者皆不會拋出例外，只會安靜地產生看似合理但錯誤的物理結果。

+ *建立可重現的驗證體系*：態製備保真度、Trotter 收斂、絕熱收斂、跨後端一致性、
  取樣收斂，全部以精確對角化為基準交叉驗證。

#danger("重要影響")[
  Givens 旋轉修正會*改變既有的極化（polarization）數值結果*，最大變化達 $+0.79$。
  任何先前產出的極化圖表都必須重新生成。詳見 @sec-impact。
]

// =====================================================================
= 核心架構與模組設計

== 模組總覽與相依關係

#tbl(
  (auto, auto, 1fr),
  (left, right, left),
  [*模組*], [*行數*], [*職責*],
  [`core.py`], [211], [自旋 SSH--Hubbard 模型、閘極原語、Givens 旋轉、退火電路],
  [`spinless.py`], [309], [無自旋 SSH + 最近鄰交互作用、`SpinlessSSHHSim`],
  [`hamiltonian.py`], [164], [任意哈密頓量檔案載入、JW 變換、位元序校正],
  [`observables.py`], [180], [扭轉不變量、密度分布、性質登錄表、取樣],
  [`backends.py`], [176], [qiskit（CPU）／cudaq（GPU）後端與自我測試],
  [`experiments.py`], [248], [五個高階實驗執行器],
  [`plotting.py`], [315], [六種出版級圖表],
  [`cli.py`], [447], [`aqs` 命令列介面],
  [`__init__.py`], [63], [公開 API 匯出],
)

相依方向為單向無環：

```
cli.py ──> experiments.py ──> spinless.py ──┐
   │            │                 │         ├──> core.py ──> (qiskit, openfermion)
   ├──> plotting.py               │         │
   └──> hamiltonian.py ───────────┴─────────┴──> observables.py
                    └──> backends.py
```

== `core.py`：閘極原語與精確態製備

=== 三個雙位元閘原語

規格書定義的三個參數化閘為：

$ R_(i,j)(theta) = exp(-i theta/2 (X_i X_j + Y_i Y_j)) $

$ G_(i,j)(theta) = exp(-i theta/2 (X_i Y_j - Y_i X_j)) $

$ "CP"_(i,j)(phi) = exp(-i phi/4 (I_i - Z_i)(I_j - Z_j)) = exp(-i phi thin n_i n_j) $

$R$ 對應*實數*躍遷項、$G$ 對應*虛數*躍遷項、$"CP"$ 對應密度--密度交互作用。

#danger("關鍵陷阱：程式中的輔助函式與規格定義相差常數因子")[
  `core.py` 的輔助函式*並非*直接實作上述定義，而是相差固定的因子與符號。
  以數值方式（`Operator` 對比 `expm`，容許全域相位）驗證得到：

  #v(0.4em)
  #tbl(
    (auto, auto, auto),
    (left, left, left),
    [*輔助函式*], [*實際么正算符*], [*與規格之差異*],
    [`create_R_gate(θ)`], [$R_"spec" (theta \/ 2)$], [因子 2],
    [`create_G_gate(θ)`], [$G_"spec" (-theta \/ 2)$], [因子 2 且反號],
    [`create_CP_gate(θ)`], [$"CP"_"spec" (-theta)$], [反號],
  )
  #v(0.4em)
  因此*絕不可*直接以規格中的角度呼叫這些函式。正確的呼叫慣例為：

  ```python
  create_R_gate(-2.0 * dt * t.real)          # 實數躍遷
  create_G_gate(-2.0 * dt * t.imag)          # 虛數躍遷（_G_SIGN = -1）
  create_CP_gate(-dt * V * lambda_l)         # 交互作用
  ```
]

=== 這些常數因子是設計選擇還是意外？

*結論：是意外，並非刻意的設計權衡。* 三個因子都沒有任何實作上的優勢，
理由如下。

*(一) 程式與自身的 docstring 互相矛盾。* `create_R_gate` 的 docstring 寫的是
`exp(-i theta/2 (XX + YY))`，這正是 $R_"spec" (theta)$；
但函式主體卻寫成 `qc.rxx(theta/2.0); qc.ryy(theta/2.0)`，實際只有一半。
若這是刻意選定的另一套慣例，docstring 應該記載該慣例才對。

*(二) 因子純粹來自包裝 qiskit 原語時未做正規化。* 由於
$"rxx"(phi) = exp(-i phi\/2 thin X X)$ 而 $X X$ 與 $Y Y$ 對易，故

$ R_"spec" (theta) = exp(-i theta/2 (X X + Y Y)) = "rxx"(theta) dot "ryy"(theta) $

也就是說，*只要把 `/2.0` 拿掉*，函式就與其 docstring 完全一致。
同理，`create_CP_gate` 只是直接轉呼叫 qiskit 的 `cp`，
而 $"cp"(lambda) = "diag"(1,1,1,e^(i lambda))$ 與
$"CP"_"spec" (phi) = e^(-i phi n_i n_j)$ 差一個負號，寫成 `cp(-theta)` 即可對齊。
三者的正規化版本皆已數值驗證：

#v(0.4em)
#tbl(
  (auto, auto, auto),
  (left, left, center),
  [*閘*], [*正規化後的主體*], [$= $ 規格式?],
  [$R$], [`rxx(t); ryy(t)`], [是],
  [$G$], [`sdg(0); rxx(-t); s(0); sdg(1); rxx(t); s(1)`], [是],
  [$"CP"$], [`cp(-t)`], [是],
)
#v(0.4em)

*(三) 正規化只會讓呼叫端更簡潔，不會更複雜。* 若三個原語都對齊規格，
Trotter 層的角度會退化成最直觀的形式，$2.0$ 這個因子直接消失：

```python
# 現況
create_R_gate(-2.0 * dt * t.real)
create_G_gate(-2.0 * dt * t.imag)
# 正規化後
create_R_gate(-dt * t.real)
create_G_gate(+dt * t.imag)
```

因此連「因子讓最常見的呼叫端變簡單」這個唯一可能的辯護理由也不成立。

#note("那為什麼文件仍保留這個警告？")[
  保留*警告*與認可*慣例*是兩回事。警告存在的理由是：

  + `create_R_gate` 的 docstring 目前具有誤導性，讀者若照字面呼叫會得到一半的角度，
    而且#underline[不會出錯、不會拋例外]，只會安靜地產生錯誤物理。
    在尚未正規化之前，這段警告是必要的補救。

  + 正規化是一次「靜默風險」重構：需同時改動 `core.py` 與 `spinless.py`
    共六處呼叫端，任何一處改錯都不會有任何測試失敗。
    本專案已經因為完全相同的失效模式產生過三個缺陷。

  建議的處理方式是*正規化，但必須先補上守門測試*：
  以 $H_"eff" = i log(U_l)\/dif t$ 擷取單步 Trotter 的有效哈密頓量，
  與小端序真值逐元比對。這個技術在本次工作中已被證實能同時抓出
  因子錯誤與虛部符號錯誤（見 @sec-lesson），是目前唯一可靠的驗證手段。
]

=== Givens 旋轉與精確 Slater 行列式態製備

初始態為 $H_"SSH"$ 的精確基態，即由最低 $N_f$ 個單粒子軌域構成的 Slater 行列式。
`openfermion.circuits.slater_determinant_preparation_circuit(Q)` 會回傳一串
$(i, j, theta, phi)$ 描述，其單粒子旋轉矩陣依 OpenFermion 慣例為：

$ G = mat(cos theta, -e^(i phi) sin theta; sin theta, e^(i phi) cos theta),
  quad det G = e^(i phi) $

提升至雙位元佔據數基底後，$lr(|00⟩)$ 不變、單佔據子空間套用 $G$、
$lr(|11⟩)$ 乘上 $det G$：

$ U_"Givens" = mat(
  1, 0, 0, 0;
  0, cos theta, sin theta, 0;
  0, -e^(i phi) sin theta, e^(i phi) cos theta, 0;
  0, 0, 0, e^(i phi)
) $

`_givens_instruction(theta, phi)` 即直接以此 $4 times 4$ 么正矩陣建構
`UnitaryGate`。此實作經 `spinless.py` 與 `core.py` 共用，避免兩份分歧的副本。

#ok("驗證")[
  對照以 OpenFermion 直接建構的 Slater 行列式
  $ product_(a=1)^(N_f) (sum_i Q_(a i) a_i^dagger) lr(|"vac"⟩) $
  在實數與複數躍遷、OBC 與 PBC、$N = 2 dots.c 4$ 的所有組合下保真度皆為 $1.0$。
]

== `hamiltonian.py`：JW 變換、靜態哈密頓量與位元序

=== 自描述的 JSON 哈密頓量格式

支援四種格式：`fermion_operator`、`qubit_operator`、`dense_matrix`、`sparse_matrix`。
前兩者經 Jordan--Wigner 變換由 `openfermion.get_sparse_operator` 轉為稀疏矩陣。

```json
{
  "type": "hamiltonian",
  "format": "fermion_operator",
  "metadata": { "n_qubits": 6, "n_cells": 3, "layout": "spinless" },
  "terms": [ {"coeff": [-0.5, 0.0], "ops": "0^ 1"} ]
}
```

=== 位元序：置換相似變換

#danger("大端序 vs 小端序")[
  `openfermion.get_sparse_operator` 將量子位元 $q$ 放在位元位置 $n-1-q$（*大端序*），
  而 qiskit 的 `Statevector`（以及 `observables.PROPERTIES` 中所有函式）
  將 $q$ 放在位元位置 $q$（*小端序*）。
]

`load_hamiltonian` 因此對算符格式套用位元反轉置換 $P$ 的相似變換：

$ H_"little" = P H_"big" P^T, quad P_(i j) = delta_(i, "bitrev"(j)) $

由於 $P$ 為置換矩陣，此變換為么正相似變換，*嚴格保持全部本徵值不變*，
僅改變基底標記。矩陣格式（`dense_matrix` / `sparse_matrix`）不做轉換，
因為其慣例由使用者自行決定。

```python
def _bit_reversal_permutation(n_qubits):
    idx = np.arange(1 << n_qubits, dtype=np.int64)
    perm = np.zeros_like(idx)
    for p in range(n_qubits):
        perm |= ((idx >> p) & 1) << (n_qubits - 1 - p)
    return perm

def _to_little_endian(H, n_qubits):
    perm = _bit_reversal_permutation(n_qubits)
    return H[perm][:, perm]
```

== `spinless.py`：延伸無自旋 SSH--Ising 模擬器

=== 量子位元映射

$2N$ 個量子位元，每個空間格點一個：

$ q_(2j) arrow.r A_j ("單位晶胞" j "的 A 次晶格"), quad
  q_(2j+1) arrow.r B_j ("單位晶胞" j "的 B 次晶格") $

鍵結分為兩類：*偶數鍵* $(2j, 2j+1)$ 為晶胞內躍遷 $v$；
*奇數鍵* $(2j+1, 2j+2)$ 為晶胞間躍遷 $w$，其中 $i = 2N-1$ 為 PBC 環繞鍵。

=== 三個主要類別與函式

#tbl(
  (auto, 1fr),
  (left, left),
  [*名稱*], [*職責*],
  [`SpinlessSSHModel`], [資料類別：單粒子矩陣 `H1`、`slater_Q()`、`fermion_operator()`、`sparse_hamiltonian()`],
  [`build_spinless_annealing_circuit`], [建構 Trotter 化絕熱電路（`steps=0` 時僅回傳態製備區塊）],
  [`SpinlessSSHHSim`], [門面類別：`build_circuit()`、`statevector()`、`measure()`、`exact_ground_state()`],
)

=== Trotter 層 $U_l$ 的建構

對每一步 $l = 1, dots, L$，令 $dif t = T \/ L$ 與斜坡權重
$lambda_l = (2l-1) \/ (2L)$，依序施加：

```python
for l in range(1, steps + 1):
    lam = (2.0 * l - 1.0) / (2.0 * steps)

    # (1) H_SSH 躍遷層：偶數鍵用 v，奇數鍵用 w
    for start, t_val in ((0, model.v), (1, model.w)):
        for i in range(start, L, 2):
            is_wrap = (i == L - 1)
            ...
            th_R = -2.0 * dt * t_re
            th_G = _G_SIGN * 2.0 * dt * t_im
            if is_wrap:                     # PBC 環繞鍵帶 JW 弦
                th_R = -parity * th_R
                th_G = -parity * th_G
            qc.append(create_R_gate(th_R), [i, j])
            qc.append(create_G_gate(th_G), [i, j])

    # (2) H_NN 交互作用層：由 lambda_l 線性斜坡
    qc.append(create_CP_gate(-dt * model.V_v * lam), [2*j, 2*j+1])
    qc.append(create_CP_gate(-dt * model.V_w * lam), [2*j+1, 2*((j+1) % N)])
```

#note("PBC 邊界的不對稱性（極易寫錯）")[
  在週期性邊界下，環繞鍵 $i = 2N-1 arrow.r 0$ 的*躍遷*項在 Jordan--Wigner 變換下
  會帶一整條穿越鏈的 $Z$ 弦，因粒子數守恆而化為固定的宇稱因子
  $-(-1)^(N_f)$，必須乘進躍遷角度。

  然而環繞鍵的*交互作用*項 $n_i n_j$ 在計算基底下是對角的，
  #underline[完全不帶 JW 弦]，因此 *`CP` 閘不可施加宇稱修正*。
  這個不對稱性是此類實作最常見的錯誤來源之一。
]

=== 精確基態參考：固定粒子數子空間

`exact_ground_state()` 預設在*固定粒子數子空間*內對角化，而非取全域基態：

```python
popcount = sum((idx >> p) & 1 for p in range(L))
sel = np.nonzero(popcount == n_particles)[0]
block = H[sel][:, sel]
```

#note("為何必須固定粒子數")[
  絕熱電路的兩類項（躍遷與密度--密度）皆守恆粒子數，故有意義的參考態是
  *同一粒子數扇區內*的基態。排斥型模型的全域基態通常落在不同填充：
  對本專案內附的 $N=3$ 交互作用範例，全域基態位於 2 個粒子而非 3 個，
  若直接使用會使比較完全失去意義。
]

== `observables.py`：量測協定與性質登錄表

=== 精確態向量評估

所有可觀測量直接以整數位元遮罩在原始振幅上向量化計算，
*刻意避開* `Statevector.probabilities_dict()`——後者會為 $2^n$ 個基態各自
產生一個字串標籤，在 $n = 24$ 時約需 5 GiB 記憶體。

=== 次晶格極化（OBC 邊緣診斷）

$ P_j^e = ⟨ n_(A,j) - n_(B,j) ⟩
  = (⟨ Z_(2j+1) ⟩ - ⟨ Z_(2j) ⟩) / 2 $

=== 多體 Berry 相位（PBC 體不變量）

以 Resta 型扭轉算符定義：

$ hat(X) = sum_j (j+1) thin n_("cell" j), quad
  z_N = ⟨ exp(i (2 pi)/N hat(X)) ⟩
  = sum_b p_b exp(i (2 pi)/N X_b) $

$ gamma = "Im" ln z_N quad (mod 2 pi) $

由於 $hat(X)$ 在計算基底下為對角算符，$z_N$ 僅取決於機率 $p_b = |c_b|^2$，
與相位無關。這個性質是取樣協定得以零成本重用同一套函式的關鍵。

=== 取樣協定：`empirical_probs`

```python
def empirical_probs(sv_np, shots, seed=None):
    probs = np.abs(np.asarray(sv_np)) ** 2
    counts = np.random.default_rng(seed).multinomial(int(shots), probs / probs.sum())
    return np.sqrt(counts / float(shots))
```

#ok("設計要點")[
  `PROPERTIES` 中的每個函式都只依賴 $|"amplitude"|^2$。因此以經驗機率
  $hat(p)$ 建構的*合成實數態向量* $sqrt(hat(p))$ 可精確重現取樣分布，
  使得整個 shot-based 協定僅需這一個輔助函式，
  #underline[無須改動任何性質函式的簽章]。
  丟棄相位正是計算基底量測的物理內涵。
]

=== 性質登錄表

```python
PROPERTIES = {
    "berry": _prop_berry,
    "polarization": _prop_polarization,
}
```

新增性質只需註冊一個 $f(#raw("sv"), #raw("layout")) arrow.r #raw("dict")$ 函式，
即可透過 `--property <name>` 在所有指令中使用。佈局（layout）由
`spin_block_layout()` 或 `spinless_layout()` 產生。

== `experiments.py`、`plotting.py` 與 `cli.py`

=== 實驗執行器

#tbl(
  (auto, 1fr),
  (left, left),
  [*執行器*], [*說明*],
  [`fidelity_scan`], [自旋模型 Trotter 收斂（參考態為製備態，自我參照）],
  [`berry_sweep`], [Berry 相位對晶胞間躍遷 $w$（PBC）],
  [`polarization_sweep`], [極化對 $Delta U$（OBC）],
  [`spinless_fidelity_scan`], [收斂至*精確交互作用基態*（真正的絕熱性檢驗）],
  [`spinless_interaction_sweep`], [$(V_v, V_w)$ 二維掃描，偵測拓樸崩潰],
)

#note("兩種 fidelity 掃描的本質差異")[
  `fidelity_scan` 以*製備態本身*為參考，屬自我參照，因此對態製備錯誤幾乎完全不敏感
  ——這正是 Givens 缺陷長期未被發現的原因。
  `spinless_fidelity_scan` 則以*精確交互作用基態*為參考，才是真正的絕熱性檢驗。
]

=== 視覺化

`plot_spinless_phase_map` 以 `pcolormesh` 繪製 $(V_v, V_w)$ 相圖，
可選量為 `Berry_Phase_pi_wrapped`、`Twist_Amplitude` 或 `edge_polarization`，
帶號量使用發散色階（`coolwarm`），非帶號量使用 `viridis`。

=== CLI 介面

```bash
# 單點量測，並與精確對角化比對
aqs spinless-measure --N 6 --v 0.5 --w 1.5 --Vv 1.0 --Vw 0.5 \
                     --T 1.0 --steps 40 --compare-exact

# 二維交互作用相圖（PBC，半填充）
aqs spinless-sweep --N 6 --v 1.0 --w 1.5 --Vv 0,1,2,3,4 --Vw 0,2 \
                   --T 1.0 --steps 40

# 邊緣響應：OBC，N+1 個粒子
aqs spinless-sweep --N 6 --obc --particles 7 --Vv 0,1,4 --Vw 0,1 \
                   --property polarization --plot-quantity edge_polarization

# 收斂掃描：T 與 L 要多大？
aqs spinless-fidelity --N 6 --Vv 1.0 --Vw 0.5 --T 1,5,15,80 --steps 10,40,100
```

躍遷振幅可為複數（`--v 0.5+0.3j`）；任何量測指令都可加 `--shots M` 與 `--seed`
切換為取樣估計；`--backend cudaq` 可切換至 GPU。

// =====================================================================
= 物理與演算法推導

== 目標靜態哈密頓量

系統為*時間無關*的延伸無自旋 SSH 模型：

$ H = H_"SSH" + H_"NN" $

=== SSH 躍遷項

$ H_"SSH" = - sum_(j=0)^(N-1) (v thin b_j^dagger a_j + v^* thin a_j^dagger b_j)
            - sum_(j) (w thin b_j^dagger a_(j+1) + w^* thin a_(j+1)^dagger b_j) $

其中 OBC 時第二個求和取 $j = 0 dots N-2$，PBC 時額外包含環繞項。

#note("符號慣例")[
  專案規格書的費米子形式寫成 $+v thin b^dagger a$，但其自身的 Pauli 算符形式卻寫成
  $-("Re" v)\/2 (X X + Y Y)$，兩者*互相矛盾*。本實作採用標準 SSH 的負號慣例
  $H = -v(dots.c) - w(dots.c)$，與 `core.SSHHModel._build_hamiltonian`、
  `examples/make_example_hamiltonian.py` 及 SSH 文獻一致。
]

=== Jordan--Wigner 變換

對相鄰模式，JW 弦相消，得到純局域的雙位元 Pauli 形式：

*以下一律採用 qiskit 的小端序慣例*（量子位元 $q$ 對應位元位置 $q$），
因為程式中所有電路與可觀測量都在此慣例下運作：

$ H_"intra" = -("Re" v)/2 sum_(j=0)^(N-1) (X_(2j) X_(2j+1) + Y_(2j) Y_(2j+1))
              +("Im" v)/2 sum_(j=0)^(N-1) (X_(2j) Y_(2j+1) - Y_(2j) X_(2j+1)) $

$ H_"inter" = -("Re" w)/2 sum_(j) (X_(2j+1) X_(2j+2) + Y_(2j+1) Y_(2j+2))
              +("Im" w)/2 sum_(j) (X_(2j+1) Y_(2j+2) - Y_(2j+1) X_(2j+2)) $

此 Pauli 分解已用數值方式驗證：對 $t = 0.7 + 0.45i$ 的兩格點鍵，
以 $"tr"(P H)\/4$ 取出的係數為 $X X: -0.35$、$Y Y: -0.35$、
$X Y: +0.225$、$Y X: -0.225$，與 $-("Re" t)\/2$ 及 $plus.minus ("Im" t)\/2$ 吻合。

#danger("虛部符號完全取決於位元序慣例")[
  若直接對 `openfermion.get_sparse_operator` 的輸出（*大端序*）做同樣的
  Pauli 分解，對兩模式系統而言模式 0 與 1 被交換，會得到
  $X Y: -0.225$、$Y X: +0.225$——#underline[虛部符號恰好相反]。

  實部項 $X X + Y Y$ 在交換下不變，因此*不會*暴露這個錯誤；
  只有虛部（複數躍遷）才會。這正是本專案 `_G_SIGN` 一度被設成錯誤符號的直接原因，
  見 @sec-lesson。任何要重新推導此符號的人，都必須先呼叫
  `hamiltonian._to_little_endian` 再做分解。
]

=== 最近鄰交互作用項

$ H_"NN" = V_v sum_(j=0)^(N-1) n_(A,j) n_(B,j)
         + V_w sum_(j) n_(B,j) n_(A,j+1) $

其中 $n_(A,j) = a_j^dagger a_j = (I - Z_(2j))\/2$，
$n_(B,j) = b_j^dagger b_j = (I - Z_(2j+1))\/2$。代入得 Ising 型：

$ H_"NN,intra" = V_v/4 sum_(j=0)^(N-1) (I_(2j) - Z_(2j))(I_(2j+1) - Z_(2j+1)) $

由於此項在計算基底下對角，其時間演化算符即為受控相位閘，
且*不帶任何 Jordan--Wigner 弦*——包含 PBC 環繞鍵在內。

== Trotter 分解

單步演化以一階 Trotter 分解近似：

$ U_l approx underbrace(product_j R_j G_j, H_"SSH" "層")
         dot underbrace(product_j "CP"_j (lambda_l), H_"NN" "層") $

需注意 $R$ 與 $G$ *不對易*（兩者分別對應同一費米子雙線性形式的
$sigma_x$ 與 $sigma_y$ 分量），故同一鍵上先後施加會引入
$cal(O)(dif t^2)$ 的 Trotter 誤差。此誤差已驗證隨 $dif t$ 二次收斂：
$dif t$ 減半時誤差約降為 $1\/4$。

== 絕熱演化協定

模擬自 $H_"SSH"$ 的精確基態出發，將交互作用強度沿 $L$ 個 Trotter 步線性升高：

$ lambda_l = (2l - 1)/(2 L), quad l = 1, dots, L $

此為*中點取樣*的線性斜坡（$lambda_(1\/2) = 0$ 至 $lambda_(L+1\/2) = 1$ 的中點），
與 `core.build_annealing_circuit` 使用的 $s = (l - 0.5)\/L$ 為同一公式。

誤差有兩個獨立來源，必須分開診斷：

#tbl(
  (auto, auto, 1fr),
  (left, left, left),
  [*誤差來源*], [*控制參數*], [*收斂行為*],
  [絕熱性（非絕熱激發）], [總演化時間 $T$], [固定 $dif t$，增大 $T$],
  [Trotter 離散化], [步長 $dif t = T \/ L$], [固定 $T$，增大 $L$],
)

#danger("常見診斷陷阱")[
  若同時增大 $T$ 與 $L$ 但*保持 $dif t$ 固定*，保真度會停滯在某個平台值，
  容易被誤判為絕熱性失效。實際上此時 Trotter 誤差為常數。
  正確做法是固定其中一項、單獨掃描另一項。
]

== 量測協定：體不變量 vs. 邊緣診斷

#tbl(
  (auto, auto, auto, 1fr),
  (left, center, center, left),
  [*可觀測量*], [*邊界*], [*填充*], [*物理意義*],
  [Berry 相位 $gamma$], [PBC], [$N$（半填充）], [體拓樸不變量],
  [次晶格極化 $P_j^e$], [OBC], [$N+1$], [邊緣模態響應],
)

#note("Berry 相位的填充／宇稱位移（極重要）")[
  以 $hat(X) = sum_j (j+1) n_j$ 定義的*絕對* Berry 相位帶有一個常數位移。
  在完全二聚化極限下，每個晶胞恰有一個電子，故
  $ X_b = sum_(j=0)^(N-1) (j+1) = N(N+1)/2 quad arrow.r.double quad
    z_N = exp(i pi (N+1)) $
  因此對*偶數* $N$，平庸相位落在 $gamma = pi$、拓樸相位落在 $gamma = 0$，
  與直覺標記恰好相反。

  *具物理意義的是跨越相變時 $pi$ 的跳躍*，而非絕對值。
  已對 $N = 4, 5, 6, 7$ 驗證：跳躍恆為 $1.0 pi$，
  且平庸側位置隨 $N$ 的奇偶交替，與上式預測完全一致。

  另注意：PBC 下極化恆等於零（無邊緣），邊緣響應必須使用 `--obc`。
]

// =====================================================================
= 關鍵缺陷修正與影響分析 <sec-impact>

本節記錄三個*靜默*缺陷。三者的共同特徵是：不會拋出例外、不會產生
NaN、輸出看起來完全合理，只是物理上錯誤。

== 缺陷一：CUDA-Q 後端位元序反轉

=== 症狀與根因

`aqs selftest --backend cudaq` 回報保真度 $0.18124364$（應為 $1.0$）。

原始程式將 qiskit 量子位元 $i$ 映射至 CUDA-Q 量子位元 $n-1-i$，
其註解宣稱「CUDA-Q 以量子位元 0 為最高有效位元」。此敘述對 CUDA-Q 的
*量測位元字串*成立，但對 `get_state()` 回傳的*態向量振幅*不成立。

直接測試三個量子位元的系統：

$ X "on" q[0] arrow.r "index" 1, quad
  X "on" q[1] arrow.r "index" 2, quad
  X "on" q[2] arrow.r "index" 4 $

即 CUDA-Q 的態向量索引為*小端序*，與 qiskit 完全一致
（已於 `nvidia`、`nvidia-fp64`、`qpp-cpu` 三個 target 確認）。

=== 修正

```python
# 修正前：q[n - 1 - qiskit_index]
# 修正後：
def cq(qiskit_index):
    return q[qiskit_index]
```

=== 影響

所有*非位元反轉對稱*電路的 GPU 結果皆為錯誤。修正後
`selftest` 保真度為 $1.00000000$，並額外以 25 組隨機非對稱電路
（2--5 量子位元，混合 `rx`/`ry`/`rz`/`cx`）驗證，最差保真度 $1.0000000000$。

== 缺陷二：OpenFermion 大端序與可觀測量小端序不一致

=== 根因

如 @sec-endian 所述，`get_sparse_operator` 為大端序而 `PROPERTIES` 為小端序，
而 `measure_hamiltonian` 未做任何轉換。

=== 診斷方法

由於專案內附的兩個範例其密度分布皆為*迴文對稱*，位元反轉對它們是恆等操作，
因此無法用來偵測此缺陷。改以刻意不對稱的哈密頓量診斷：在四模式鏈的模式 0
施加強吸引位能，基態應有 $n_0 approx 1$：

```python
op = of.FermionOperator('0^ 0', -5.0)   # 只在模式 0 設陷阱
...
# 修正前：density = [0.2555, 0.5034, 0.2423, 0.9988]  -> argmax = 3  (錯)
# 修正後：density = [0.9988, 0.2423, 0.5034, 0.2555]  -> argmax = 0  (對)
```

=== 修正與影響 <sec-endian>

`load_hamiltonian` 現以置換相似變換回傳小端序算符。已驗證本徵譜完全不變。
影響範圍為 `aqs measure` 的 `fermion_operator` 與 `qubit_operator` 格式；
矩陣格式不受影響。

== 缺陷三：Givens 旋轉符號錯誤——製備出能量最高的行列式態

=== 根因

原始 `_givens_instruction` 以下列電路分解實作：

```python
gv.rz(phi, 1); gv.rz(-phi, 0); gv.cx(0, 1)
gv.cry(-2.0 * theta, 1, 0)          # <-- 應為 +2.0 * theta
gv.cx(0, 1)
```

存在兩個獨立缺陷：

+ *`cry` 角度反號*：導致態製備建構出由*最高*單粒子軌域組成的行列式態。
  在具手徵對稱性的 SSH 鏈上，最高與最低軌域的能量互為精確負值，
  因此製備態的能量為 $+E$ 而非基態能量 $-E$。

+ *`rz` 結構對複數躍遷完全錯誤*：當 $phi eq.not 0$ 時無法重現 OpenFermion
  的相位慣例。此路徑目前無 CLI 可觸及（`--v`/`--w` 僅接受實數）。

=== 定量證據

以製備態的單體能量期望值 $⟨H⟩$ 對照最低／最高填充能量：

#tbl(
  (auto, auto, auto, auto, auto),
  (left, right, right, right, right),
  [*組態*], [*修正前 $E_"prep"$*], [*修正後 $E_"prep"$*], [*最低填充*], [*最高填充*],
  [$N=3$, $v=1.0$, $w=1.5$, PBC], [$+10.291503$], [$-10.291503$], [$-10.291503$], [$+10.291503$],
  [$N=3$, $v=0.5$, $w=1.5$, OBC], [$+6.416289$], [$-6.416289$], [$-6.416289$], [$+6.416289$],
  [$N=4$, $v=1.0$, $w=0.5$, PBC], [$+8.472136$], [$-8.472136$], [$-8.472136$], [$+8.472136$],
)

修正後 $E_"prep"$ 在所有組態下皆*精確等於*最低填充能量。

=== 為何長期未被發現

`fidelity_scan` 以*製備態本身*為參考比較演化後的態，屬自我參照的比較，
因此對「製備了哪一個行列式態」幾乎完全不敏感——修正前後其輸出
#underline[逐位元完全相同]。

=== 量化影響

#tbl(
  (auto, 1fr),
  (left, left),
  [*結果*], [*是否受影響*],
  [`aqs fidelity`], [*否*。自我參照，數值逐位元相同。],
  [`aqs berry`（$Delta U = 0$）], [*否*。粒子--電洞對稱性使上下能帶行列式態的扭轉不變量相同。],
  [`aqs berry`（$Delta U eq.not 0$）], [*是，但幅度小*，典型約 $0.03 pi$。],
  [`aqs polarization`], [*是，且影響顯著*。邊緣極化最大變化 $+0.79$，且分布形狀改變。],
)

極化的完整分布變化（$N=3$，$v=0.5$，$w=1.5$，$U_A=1.0$，OBC）：

#tbl(
  (auto, auto, auto),
  (right, left, left),
  [$Delta U$], [*修正前*], [*修正後*],
  [0], [`[0.9259, 0.0, -0.9259]`], [`[0.8351, 0.0, -0.8351]`],
  [0.1], [`[0.908, -0.0308, -0.9481]`], [`[0.8503, 0.025, -0.8152]`],
  [1], [`[0.7429, -0.3126, -1.144]`], [`[0.9771, 0.2325, -0.6449]`],
  [3], [`[0.3867, -0.905, -1.5111]`], [`[1.1807, 0.5623, -0.3432]`],
)

#danger("需重新生成的資料")[
  在 $Delta U = 3$ 時，中間晶胞的極化由 $-0.905$ 變為 $+0.5623$——*正負號翻轉*。
  任何先前發表的極化圖表都必須重新生成。
  $Delta U = 0$ 的 Berry 相位圖不受影響。
]

#note("關於最大的 Berry 相位變化")[
  最大單點變化出現在 $N=3$、$Delta U = 0.1$、$w = 2.0$，達 $+0.69 pi$。
  然而該點的扭轉振幅 $|z_N| = 0.0088$，
  即此相位是一個趨近於零的複數之輻角，
  #underline[修正前後皆不具物理意義]，不應據此解讀。
]

== 一則方法論教訓：驗證本身也可能有位元序錯誤 <sec-lesson>

在本次工作中，作者最初曾「證實」`core.py` 的 $G$ 閘符號有誤並據此設定
`_G_SIGN = +1`。該結論*本身是錯的*：驗證時直接使用了 OpenFermion 的
大端序矩陣而未轉換，對兩模式系統而言這會交換模式 0 與 1，
並悄悄地將待檢查的非對角元共軛。

真正的破綻是：複數躍遷電路的保真度停滯在 $0.80$ 且*不隨 $L$ 增大而改善*
——這代表系統性誤差而非離散化誤差。改以擷取單步 Trotter 的有效哈密頓量

$ H_"eff" = i log(R dot G) \/ dif t $

與真實鍵結哈密頓量在小端序基底下逐元比對，才確認虛部被共軛：
真值 $-0.7 - 0.45i$，電路給出 $-0.7 + 0.45i$。最終 `_G_SIGN = -1`，
與 `core.build_annealing_circuit` 原本的慣例一致。

#note("防呆註解")[
  此陷阱已寫入 `spinless.py` 的 `_G_SIGN` 註解，明確警告未來若要重新推導此符號，
  必須以小端序哈密頓量（`hamiltonian._to_little_endian`）為基準，否則會再次
  「確認」錯誤的符號。
]

// =====================================================================
= 驗證與基準測試

== 驗證矩陣

#tbl(
  (auto, 1fr, auto),
  (left, left, center),
  [*編號*], [*驗證項目*], [*結果*],
  [V1], [非交互作用極限：$gamma$ 跨 $w = v$ 跳躍 $pi$（$N = 4 dots 7$）], [PASS],
  [V2], [態製備 $=$ 精確非交互作用基態（實／複數、OBC／PBC）], [$1.0000000000$],
  [V3], [絕熱收斂至精確交互作用基態], [$>= 0.9997$],
  [V4], [跨後端一致性（qiskit vs. CUDA-Q）], [$< 10^(-4)$],
  [V5], [取樣收斂 $hat(z) arrow.r z$（$M arrow.r infinity$）], [$prop 1\/sqrt(M)$],
  [V6], [既有指令無回歸], [PASS],
  [V7], [目標運行 $N=6$、$T=1$、$L=40$], [PASS],
)

== 態製備保真度（V2）

#tbl(
  (auto, auto, auto, auto),
  (left, left, center, right),
  [$N$], [$v$], [*邊界*], [*保真度*],
  [6], [$0.5$], [PBC], [$1.0000000000$],
  [3], [$0.5$], [OBC], [$1.0000000000$],
  [4], [$0.4 + 0.9i$], [OBC], [$1.0000000000$],
  [3], [$0.6 + 0.2i$], [PBC], [$1.0000000000$],
)

== 電路正確性：與精確演化比對（`ramp=False`，$T = 0.4$）

固定 $T$ 增大 $L$，檢驗 Trotter 誤差是否隨 $dif t$ 收斂至零：

#tbl(
  (auto, auto, auto, auto, auto),
  (left, center, right, right, right),
  [*躍遷*], [*邊界*], [$L = 50$], [$L = 200$], [$L = 800$],
  [實數], [PBC], [$0.999977$], [$0.99999858$], [$0.9999999111$],
  [實數], [OBC], [$0.999986$], [$0.99999911$], [$0.9999999441$],
  [複數], [PBC], [$0.999966$], [$0.99999786$], [$0.9999998662$],
  [複數], [OBC], [$0.999978$], [$0.99999861$], [$0.9999999129$],
)

== 兩類誤差的分離診斷（$N=3$、$V_v=1.0$、$V_w=0.5$、OBC）

#grid(
  columns: (1fr, 1fr),
  gutter: 12pt,
  [
    *(a) 固定 $dif t = 0.025$，增大 $T$*（絕熱性）
    #v(0.3em)
    #tbl(
      (auto, auto, auto),
      (right, right, right),
      [$T$], [$L$], [*保真度*],
      [1], [40], [$0.97936981$],
      [2], [80], [$0.99205010$],
      [5], [200], [$0.99959422$],
      [10], [400], [$0.99975783$],
      [20], [800], [$0.99988331$],
    )
  ],
  [
    *(b) 固定 $T = 20$，增大 $L$*（Trotter）
    #v(0.3em)
    #tbl(
      (auto, auto, auto),
      (right, right, right),
      [$L$], [$dif t$], [*保真度*],
      [80], [$0.2500$], [$0.98689181$],
      [160], [$0.1250$], [$0.99700806$],
      [400], [$0.0500$], [$0.99956161$],
      [800], [$0.0250$], [$0.99988331$],
      [1600], [$0.0125$], [$0.99994753$],
    )
  ],
)

(b) 欄相鄰列的誤差比值依序為 $4.38, 6.82, 3.76, 2.22$，
對照 $cal(O)(dif t^2)$ 預期的 $(L_(i+1)\/L_i)^2 = 4, 6.25, 4, 4$
（$L$ 的遞增倍率並非固定，第二步為 $2.5$ 倍）。
前三步吻合良好；最後一步偏低是因為誤差已降至 $5 times 10^(-5)$，
逼近數值精度下限，並非收斂階數改變。

== 規格指定的基準點：$T = 1.0$、$L = 40$

在 $N = 6$（12 量子位元）、$V_v = 1.0$、$V_w = 0.5$ 下，
與精確交互作用基態的保真度為 $0.99934891$。

`aqs spinless-fidelity` 的收斂表（$N=6$）：

#tbl(
  (auto, auto, auto, auto),
  (right, right, right, right),
  [$T$], [$L = 10$], [$L = 40$], [$L = 100$],
  [1], [$0.993$], [$0.999$], [$1.000$],
  [5], [$0.944$], [$0.996$], [$0.999$],
  [15], [$0.129$], [$0.919$], [$0.987$],
  [80], [$0.005$], [$0.065$], [$0.704$],
)

#note("如何解讀此表")[
  $T$ 增大而 $L$ 固定時 $dif t$ 隨之增大，故右下角的低保真度反映的是
  *Trotter 誤差*而非絕熱性不足。規格指定的 $(T, L) = (1.0, 40)$
  落在收斂良好的區域。
]

== 跨後端驗證（V4）

#tbl(
  (auto, auto, auto, auto),
  (left, right, right, right),
  [*後端*], [$gamma$ ($pi$)], [$|z_N|$], [*邊緣極化*],
  [qiskit（CPU）], [$-0.000000$], [$0.308828$], [$+0.000000$],
  [cudaq（GPU）], [$+0.000000$], [$0.308828$], [$+0.000000$],
)

#note("精度說明")[
  CUDA-Q 的 `nvidia` target 預設為 *fp32 單精度*，
  故兩後端的差異約在 $10^(-4)$ 量級（如 $0.8244$ vs. $0.8243$）。
  這是精度而非缺陷。如需雙精度可使用
  `cudaq.set_target('nvidia', option='fp64')`。
]

== 取樣收斂（V5）

精確值 $|z_N| = 0.308828$，每組取 5 個獨立種子：

#tbl(
  (auto, auto, auto),
  (right, right, right),
  [*shots $M$*], [*平均 $|z_N|$*], [*標準差*],
  [$10^2$], [$0.334403$], [$0.030897$],
  [$10^3$], [$0.317566$], [$0.008263$],
  [$10^4$], [$0.312001$], [$0.008533$],
  [$10^5$], [$0.309156$], [$0.000888$],
  [$10^6$], [$0.308459$], [$0.000616$],
)

誤差隨 $1\/sqrt(M)$ 下降，符合預期的散粒雜訊標度。

== 拓樸相變崩潰邊界（V7）

在 $N = 6$（12 量子位元）、$v = 1.0$、$w = 1.5$（拓樸側）、
$T = 1.0$、$L = 40$、PBC 半填充下的 $(V_v, V_w)$ 掃描：

#tbl(
  (auto, auto, auto, auto),
  (right, right, right, right),
  [$V_v$], [$V_w$], [$gamma$ ($pi$)], [$|z_N|$],
  [0], [0], [$-0.0000$], [$0.3148$],
  [0], [4], [$-0.0000$], [$0.4349$],
  [1], [0], [$+0.0000$], [$0.1988$],
  [2], [0], [$-0.0000$], [$0.0752$],
  [2], [4], [$-0.0000$], [$0.1943$],
  [4], [0], [$+1.0000$], [$0.1645$],
  [4], [2], [$+1.0000$], [$0.1374$],
  [4], [4], [$+1.0000$], [$0.0954$],
)

#ok("物理解讀")[
  Berry 相位在 $V_v <= 2$ 時維持量子化值 $0$，於 $V_v = 4$ 翻轉至 $pi$，
  且*扭轉振幅 $|z_N|$ 沿同一邊界塌縮*（$V_w = 0$ 時由 $0.3148$ 降至 $0.0752$）。

  這與物理直覺一致：強晶胞內排斥 $V_v$ 懲罰同一晶胞內的雙佔據，
  使系統偏好晶胞內二聚化的平庸組態，因而驅動拓樸相崩潰。
  相對地，$V_w$ 對相位邊界的影響微弱得多。

  兩個獨立訊號（相位翻轉與振幅塌縮）在同一邊界一致，
  是此結果可信度的重要佐證。
]

// =====================================================================
= 附錄

== 快速上手

```bash
# 環境建置（uv）
uv sync --extra cu12          # CUDA 12.x 驅動
uv sync --extra cu13          # CUDA 13.x 驅動
uv sync                       # 純 CPU

# 後端驗證
uv run aqs selftest --backend cudaq
```

== Python API 範例

```python
from aqs import SpinlessSSHHSim

sim = SpinlessSSHHSim(N_cells=6, v=0.5, w=1.5, V_v=1.0, V_w=0.5, PBC=True)

res = sim.measure(["berry", "polarization"], T=1.0, steps=40)
print(res["berry"]["Berry_Phase_pi_wrapped"])

# 與精確對角化交叉驗證（固定粒子數子空間）
psi_exact = sim.exact_ground_state()

# 取樣協定（規格書 §5）
res_shots = sim.measure(["berry"], T=1.0, steps=40, shots=100_000, seed=42)
```

== 缺陷修正一覽

#tbl(
  (auto, 1fr, auto),
  (left, left, left),
  [*Commit*], [*內容*], [*狀態*],
  [`1e05232`], [CUDA-Q 後端位元序反轉], [已修正],
  [`d692254`], [OpenFermion 大端序 / Givens 旋轉符號], [已修正],
  [`2be4328`], [新增無自旋模組與可觀測量擴充], [新功能],
  [`91daae1`], [實驗執行器、相圖與 CLI 子指令], [新功能],
  [`7f9dce3`], [文件與公開 API 匯出], [新功能],
)

== 已知限制與後續工作

- CUDA-Q `nvidia` target 預設 fp32；`--fp64` 旗標尚未接入 CLI。
- 專案目前*沒有自動化測試套件*，回歸驗證仰賴 `aqs selftest` 與本文件所述的
  手動交叉驗證流程。建議將 V1--V7 固化為 `pytest` 測試。
- `dist/` 內的 wheel 為 2.0.0 版本，與 `src/` 可能不同步。
- 絕對 Berry 相位帶有 $pi(N+1)$ 的填充／宇稱位移，目前以文件說明而非程式校正。
