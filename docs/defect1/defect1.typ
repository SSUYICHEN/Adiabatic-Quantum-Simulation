// ╔══════════════════════════════════════════════════════════════╗
// ║   EDA Lab · Quantum Group   Defect Report #1                  ║
// ║   Givens 旋轉反號：Slater 行列式態製備出最高能帶                ║
// ║   Template: weekly-report.typ · Theme: themes/JellyFish        ║
// ║   Figures: python docs/defect1/make_figures.py                 ║
// ╚══════════════════════════════════════════════════════════════╝

#import "@preview/physica:0.9.7": *

// ── Cover Colour Palette (used ONLY on cover) ─────────────────
#let indigo = rgb("#4F46E5")
#let indigo-lt = rgb("#EEF2FF")
#let indigo-mid = rgb("#C7D2FE")
#let green-acc = rgb("#059669")
#let amber-acc = rgb("#D97706")

// ── Inner Content Palette (grayscale) ─────────────────────────
#let black = rgb("#1A1A1A")
#let dark-gray = rgb("#333333")
#let mid-gray = rgb("#666666")
#let border-gray = rgb("#CCCCCC")
#let light-gray = rgb("#F5F5F5")
#let white = rgb("#FFFFFF")

// ── Code Block Colours ────────────────────────────────────────
#let red-acc = rgb("#DC2626")
#let red-button = rgb("#FF5F58")
#let yellow-button = rgb("#FFBD2D")
#let green-button = rgb("#27C840")
#let code-bg = rgb("#00002c")
#let code-fg = rgb("#EEFFFF")
#let code-line = rgb("#0E0021")

// ── Figure colours（與圖片一致：橘 = 修正前、藍 = 修正後）──────
#let c-before = rgb("#c2410c")
#let c-after = rgb("#1d5fb3")
#let before = text(fill: c-before, weight: "bold")[修正前]
#let after = text(fill: c-after, weight: "bold")[修正後]

// ── Report Metadata ───────────────────────────────────────────
#let author = "鄧恩陞"
#let group = "Quantum Group"
#let lab = "ALCom Lab"
#let report-title = "Defect Report"
#let week-no = "Defect #1"
#let date-str = "2026-09-17"
#let advisor = "Prof. 江介宏"
#let subtitle = "Givens 旋轉反號：Slater 行列式態製備出最高能帶"

// ══════════════════════════════════════════════════════════════
//  Global document settings
// ══════════════════════════════════════════════════════════════
#set raw(theme: "themes/JellyFish.tmTheme")
#set document(author: author, title: week-no + " " + report-title)

#set page(
  paper: "a4",
  margin: (top: 2.4cm, bottom: 2.4cm, left: 2.6cm, right: 2.6cm),
  header: context {
    if counter(page).get().first() > 1 {
      set text(size: 8pt, fill: border-gray)
      grid(
        columns: (1fr, auto),
        [#lab · #group · #week-no], [#author],
      )
      line(length: 100%, stroke: 0.5pt + border-gray)
    }
  },
  footer: context {
    set text(size: 8pt, fill: border-gray)
    line(length: 100%, stroke: 0.5pt + border-gray)
    v(2pt)
    grid(
      columns: (1fr, auto),
      [#date-str], [#counter(page).display("1 / 1", both: true)],
    )
  },
)

#set text(
  font: ("Linux Libertine O", "Noto Serif CJK TC", "Noto Serif"),
  size: 11pt,
  fill: black,
  lang: "zh",
  region: "TW",
)
#set par(justify: true, leading: 0.85em, first-line-indent: 0em)

#set table(
  stroke: none,
  inset: (x: 10pt, y: 5pt),
  align: center,
)
#show table: set text(size: 10pt)
#show figure.caption: set text(size: 9.5pt, fill: dark-gray)
#show figure: set block(above: 1.2em, below: 1.4em)

#set heading(numbering: "1.1.1")
#show heading.where(level: 1): it => {
  v(18pt)
  block(below: 12pt)[
    #text(weight: "bold", size: 14.4pt, fill: black)[
      #counter(heading).display("1.")
      #h(0.4em)
      #it.body
    ]
  ]
}
#show heading.where(level: 2): it => {
  v(16pt)
  block(below: 8pt)[
    #text(weight: "bold", size: 12pt, fill: black)[
      #counter(heading).display("1.1")
      #h(0.4em)
      #it.body
    ]
  ]
}
#show heading.where(level: 3): it => {
  v(16pt)
  block(below: 8pt)[
    #text(weight: "bold", size: 11pt, fill: black)[
      #counter(heading).display("1.1.1")
      #h(0.4em)
      #it.body
    ]
  ]
}

// ── Modern Code Block ─────────────────────────────────────────
#show raw.where(block: true): it => {
  block(
    width: 100%,
    radius: 8pt,
    clip: true,
    stroke: 1pt + code-line,
    breakable: false,
  )[
    #set block(spacing: 0pt)
    #block(
      width: 100%,
      fill: code-line,
      inset: (x: 12pt, y: 6pt),
    )[
      #set text(size: 8pt, fill: code-fg.lighten(20%))
      #stack(
        dir: ltr,
        spacing: 6pt,
        circle(radius: 4pt, fill: red-button),
        circle(radius: 4pt, fill: yellow-button),
        circle(radius: 4pt, fill: green-button),
        h(1fr),
        if it.lang != none { text(fill: indigo-mid)[#it.lang] },
      )
    ]
    #block(
      width: 100%,
      fill: code-bg,
      inset: (x: 14pt, y: 12pt),
    )[
      #set text(
        font: ("MesloLGS Nerd Font Mono", "DejaVu Sans Mono", "Noto Sans Mono CJK TC"),
        size: 9.5pt,
        fill: code-fg,
      )
      #set par(leading: 1.2em)
      #it
    ]
  ]
}

#show raw.where(block: false): it => {
  box(
    fill: light-gray,
    stroke: 0.5pt + border-gray,
    radius: 3pt,
    inset: (x: 4pt, y: 2pt),
    text(font: ("MesloLGS Nerd Font Mono", "DejaVu Sans Mono", "Noto Sans Mono CJK TC"), size: 9pt, fill: dark-gray)[#it],
  )
}

// ── Callout Boxes ─────────────────────────────────────────────
#let callout(title: "", color: mid-gray, icon: "ℹ", body) = {
  block(
    width: 100%,
    radius: 4pt,
    stroke: (left: 3pt + color, rest: 0.5pt + color.lighten(60%)),
    fill: color.lighten(90%),
    inset: (x: 12pt, y: 10pt),
    breakable: false,
  )[
    #text(fill: color.darken(20%), weight: "bold", size: 9.5pt)[#icon #h(3pt) #title]
    #v(4pt)
    #body
  ]
}
#let info(title: "Note", body) = callout(title: title, color: rgb("#6B7280"), icon: "ℹ", body)
#let warn(title: "Problem", body) = callout(title: title, color: rgb("#D97706"), icon: "⚠", body)
#let done(title: "Done", body) = callout(title: title, color: rgb("#059669"), icon: "✓", body)

#let badge-cover(label, color: indigo) = box(
  fill: color.lighten(80%),
  stroke: 0.5pt + color,
  radius: 10pt,
  inset: (x: 8pt, y: 3pt),
  text(size: 8.5pt, fill: color, weight: "bold")[#label],
)

// booktabs-style table helper
#let btable(cols, header, ..rows) = block(width: 100%)[
  #line(length: 100%, stroke: 1pt + black)
  #table(
    columns: cols,
    table.header(..header),
    table.hline(stroke: 0.5pt + black),
    ..rows,
  )
  #v(-8pt)
  #line(length: 100%, stroke: 1pt + black)
]

// ══════════════════════════════════════════════════════════════
//  COVER PAGE
// ══════════════════════════════════════════════════════════════
#let slate = rgb("#334155")

#page(
  margin: (top: 0pt, bottom: 0pt, left: 0pt, right: 0pt),
  header: none,
  footer: none,
)[
  #block(
    width: 100%,
    height: 38%,
    fill: gradient.linear(indigo, rgb("#7C3AED"), angle: 135deg),
  )[
    #place(bottom + left, dx: 52pt, dy: -30pt)[
      #text(fill: rgb("#FFFFFF").transparentize(20%), size: 72pt, weight: "bold")[Q]
      #h(-12pt)
      #text(fill: rgb("#FFFFFF"), size: 28pt, weight: "bold")[uantum]
      #linebreak()
      #text(fill: rgb("#FFFFFF").transparentize(30%), size: 13pt)[#lab · #group]
    ]
  ]

  #place(top + left, dx: 46pt, dy: 32%)[
    #block(
      width: 520pt,
      fill: rgb("#FFFFFF"),
      radius: 12pt,
      stroke: 1.5pt + indigo-mid,
      inset: (x: 36pt, y: 30pt),
    )[
      #text(size: 26pt, weight: "bold", fill: indigo)[#report-title]
      #v(2pt)
      #text(size: 13pt, fill: slate)[#subtitle]
      #v(4pt)
      #line(length: 100%, stroke: 1.5pt + indigo-mid)
      #v(14pt)

      #grid(
        columns: (auto, 1fr),
        column-gutter: 12pt,
        row-gutter: 10pt,
        text(fill: slate.lighten(30%), size: 9.5pt)[ID], text(weight: "bold", size: 15pt, fill: slate)[#week-no],
        text(fill: slate.lighten(30%), size: 9.5pt)[DATE], text(size: 10pt)[#date-str],
        text(fill: slate.lighten(30%), size: 9.5pt)[AUTHOR], text(size: 10pt)[#author],
        text(fill: slate.lighten(30%), size: 9.5pt)[ADVISOR], text(size: 10pt)[#advisor],
        text(fill: slate.lighten(30%), size: 9.5pt)[GROUP], text(size: 10pt)[#lab · #group],
        text(fill: slate.lighten(30%), size: 9.5pt)[REPO], text(size: 10pt)[`Adiabatic-Quantum-Simulation` · 分支 `fix/core-silent-convention-bugs`],
      )

      #v(18pt)
      #badge-cover("Slater 態製備", color: indigo)
      #h(4pt)
      #badge-cover("Givens 旋轉", color: green-acc)
      #h(4pt)
      #badge-cover("論文影響評估", color: amber-acc)
    ]
  ]

  #place(bottom + center, dy: -16pt)[
    #text(size: 8pt, fill: slate.lighten(50%))[
      Confidential · #lab · #datetime.today().display("[year]")
    ]
  ]
]

// ══════════════════════════════════════════════════════════════
//  MAIN CONTENT
// ══════════════════════════════════════════════════════════════

= 📋 摘要

#warn(title: "缺陷一句話")[
  `main` 分支的 `_givens_instruction` 用 `rz / cx / cry / cx` 手工分解 Givens 旋轉，
  其中 `cry` 的角度多了一個負號。實數躍遷下它做出的是 $G(-theta)$ 而不是 $G(theta)$，
  結果 Slater 行列式態製備出的是單粒子能譜中*最高*的一組軌域，而非基態。
]

#done(title: "本報告確認的事實")[
  - 以最簡單的氫分子模型（兩原子、兩電子）測試，修正前製備出反鍵結態，能量 $+2t$；修正後為基態 $-2t$。
  - 以矩陣推導證明：舊分解等於 $G(-theta, dot)$，且 $phi eq.not 0$ 時漏掉 $ket(11)$ 上的 $e^(i phi)$ 相位。
  - 舊碼的整條電路輸出，與「修正後程式、但 $U_A, U_B$ 全部反號」逐位元相同。
    因此論文 Fig. 3、Fig. 4 中 $Delta U eq.not 0$ 的曲線畫的是吸引型不平衡的物理，卻標成排斥型。
  - Fig. 2 與 $Delta U = 0$ 的曲線不受影響；「弱不平衡下拓樸特徵穩健」的定性結論成立。
]

= 🧩 背景：OpenFermion 的 Slater 行列式態製備

== 演算法的輸入與輸出

要在量子電路上製備一個由 $nu$ 個已知單粒子軌域填滿的 Slater 行列式態

$
  ket(Psi) = b_0^dagger b_1^dagger dots.c b_(nu-1)^dagger ket("vac"), quad
  b_j^dagger = sum_(k=0)^(N-1) Q_(j k) a_k^dagger ,
$

OpenFermion 提供 `slater_determinant_preparation_circuit(Q)`：

- *輸入*：$nu times N$ 的軌域係數矩陣 $Q$，每一列是一個被佔據的軌域，滿足 $Q Q^dagger = I_nu$（$N$ 為軌域總數、$nu$ 為電子數）。
- *輸出*：一串兩模式 Givens 旋轉 ${(j_1, k_1, theta_1, phi_1), (j_2, k_2, theta_2, phi_2), dots}$，分成可以平行執行的層。

演算法本質是對 $Q$ 做 QR 型的消去：每一個 Givens 旋轉把 $Q$ 的一個元素消成零，最後只剩前 $nu$ 個模式被佔據。
把消去順序反過來執行，就從 $ket(1 dots 1 0 dots 0)$ 製備出目標態。
只在被佔據軌域之間旋轉不改變行列式（只差整體相位），所以這個程序是精確的。

#figure(
  image("figures/fig_workflow.svg", width: 100%),
  caption: [態製備流程。OpenFermion 只回傳角度序列，如何把每一組 $(theta, phi)$ 變成兩量子位元閘是使用者的責任，也正是缺陷所在。],
)

== 每一個 Givens 旋轉應該是什麼

在兩個 qubit 的佔據基底 ${ket(00), ket(01), ket(10), ket(11)}$（寫法 $ket(q_1 q_0)$，右邊是低位 qubit）上，
OpenFermion 定義的旋轉是

$
  G(theta, phi) = mat(
    1, 0, 0, 0;
    0, cos theta, sin theta, 0;
    0, -e^(i phi) sin theta, e^(i phi) cos theta, 0;
    0, 0, 0, e^(i phi);
  ).
$

三個特徵缺一不可：$ket(00)$ 不動；單佔據子空間 ${ket(01), ket(10)}$ 做一個帶相位的旋轉；
$ket(11)$ 乘上 $det G = e^(i phi)$。最後一項容易被忽略，但複數軌域時它決定了行列式的相位一致性。

#figure(
  image("figures/fig_circuit_h2.svg", width: 62%),
  caption: [氫分子模型（4 個模式、2 個電子）的製備電路：先以 X 閘放入粒子，再對每個自旋區塊套用 OpenFermion 給出的 Givens 旋轉。],
)

= 🔬 兩種實作

== 修正前：`main` 分支的閘分解

```python
# ---- 修正前：逐字抄自 main 0fe9d1d 的 core._givens_instruction ----
def legacy_givens(theta, phi):
    gv = QuantumCircuit(2, name="Givens")
    gv.rz(phi, 1)
    gv.rz(-phi, 0)
    gv.cx(0, 1)
    gv.cry(-2.0 * theta, 1, 0)      # <-- 角度反號
    gv.cx(0, 1)
    return gv.to_instruction()
```

#figure(
  image("figures/fig_circuit_legacy.svg", width: 58%),
  caption: [修正前的分解（$theta = 0.7$、$phi = 0.9$）。`cry` 的旋轉角是 $-2theta$，這個負號就是整個缺陷的來源。],
)

== 修正後：直接寫下 $4 times 4$ 矩陣

```python
def _givens_instruction(theta, phi):
    c, s = np.cos(theta), np.sin(theta)
    e = np.exp(1j * phi)
    U = np.zeros((4, 4), dtype=complex)
    U[0, 0] = 1.0
    U[1, 1] = c
    U[2, 2] = c * e
    U[1, 2] = s
    U[2, 1] = -s * e
    U[3, 3] = e
    return UnitaryGate(U, label="Givens")
```

#figure(
  image("figures/fig_circuit_fixed.svg", width: 42%),
  caption: [修正後不再手工分解，直接把 $G(theta, phi)$ 當成一個 `UnitaryGate` 交給 qiskit 轉譜，避免任何符號或相位慣例錯誤。],
)

= 🧪 氫分子模型：兩者製備出的態不同

最簡單的檢驗：兩個原子 A、B，各有自旋上下兩個模式，共 4 個模式、2 個電子，
$H = -t sum_s (a_s^dagger b_s + b_s^dagger a_s)$。每個自旋的基態是鍵結軌域 $(ket(A) + ket(B)) \/ sqrt(2)$，
所以 $Q = (1, 1) \/ sqrt(2)$，目標態是「鍵結 $⊗$ 鍵結」，能量 $-2t$。

#figure(
  image("figures/fig_h2_amplitudes.svg", width: 92%),
  caption: [四個雙電子基底態上的振幅。修正前在 $ket(A arrow.t B arrow.b)$ 與 $ket(B arrow.t A arrow.b)$ 上的符號相反，那是「反鍵結 $⊗$ 反鍵結」，能量 $+2t$；修正後與參考態完全一致。],
)

#figure(
  image("figures/fig_h2_density.svg", width: 62%),
  caption: [兩個態的佔據數 $expval(n_q)$ 完全相同，每個模式都是 0.5。只看密度或粒子數，永遠看不出這個錯誤。],
)

#btable(
  (auto, auto, auto, auto),
  ([], [*$expval(H) \/ t$*], [*與參考態的 fidelity*], [*判定*]),
  [#before], [$+2.000$], [$0.000$], [反鍵結態],
  [#after], [$-2.000$], [$1.000$], [基態],
)

= 📐 數學檢查：舊分解到底做出了什麼

把修正前的電路逐閘寫成矩陣（基底順序 $ket(q_1 q_0)$，列為輸出、欄為輸入）。
兩個 $R_z$ 合起來是對角相位

$
  U_(R_z) = mat(1, 0, 0, 0; 0, e^(-i phi), 0, 0; 0, 0, e^(i phi), 0; 0, 0, 0, 1),
$

中間的 CNOT–$C R_y(-2theta)$–CNOT 三明治，在單佔據子空間上等於一個實旋轉。
把它乘開：

$
  U_"legacy" =
  mat(1, 0, 0, 0; 0, cos theta, -sin theta, 0; 0, sin theta, cos theta, 0; 0, 0, 0, 1)
  mat(1, 0, 0, 0; 0, e^(-i phi), 0, 0; 0, 0, e^(i phi), 0; 0, 0, 0, 1)
  =
  mat(1, 0, 0, 0;
      0, e^(-i phi) cos theta, -e^(i phi) sin theta, 0;
      0, e^(-i phi) sin theta, e^(i phi) cos theta, 0;
      0, 0, 0, 1).
$

與 $G(theta, phi)$ 對照有兩處不同：

+ *單佔據區塊的 $sin theta$ 符號相反*：$G$ 的 $(01, 10)$ 元是 $+sin theta$、$(10, 01)$ 元是 $-e^(i phi) sin theta$，
  舊分解剛好反過來。這正是 $theta arrow.r -theta$。
+ *$ket(11)$ 上的相位*：$G$ 給 $e^(i phi)$，舊分解給 $1$；單佔據區塊的相位分布也不同（$e^(-i phi)$ 出現在不該出現的地方）。
  沒有任何 $(theta, phi)$ 的重新參數化能同時修好這兩處。

#figure(
  image("figures/fig_matrix_compare.svg", width: 100%),
  caption: [數值驗證（$theta = 0.7$、$phi = 0.9$）：左為 OpenFermion 定義的 $G$，中為舊分解實際的矩陣，右為兩者之差。非對角元與 $ket(11)$ 相位都不同。],
)

#figure(
  image("figures/fig_matrix_phi0.svg", width: 100%),
  caption: [實數躍遷（$phi = 0$）時舊分解恰好等於 $G(-theta, 0)$，逐元差為 0。論文的 SSH 模型躍遷全為實數，所以缺陷以「角度反號」這個最單純的形式出現。],
)

#info(title: "為什麼反號會變成「最高能帶」")[
  SSH 鏈有手徵對稱：單粒子能階成對出現 $plus.minus epsilon$，而且把 $theta arrow.r -theta$ 相當於在 A、B 子格點上翻轉相對符號，
  把每個 $-epsilon$ 的軌域換成對應的 $+epsilon$ 軌域。因此舊碼不是「差一點」，而是把該填的最低 $nu$ 個軌域整組換成最高 $nu$ 個。
]

#figure(
  image("figures/fig_spectrum.svg", width: 92%),
  caption: [SSH 鏈（$N = 3$、$v = 0.5$、$w = 1.5$、OBC、8 電子）的單粒子能階。修正後每個自旋填最低 4 個軌域，$expval(H) = -6.416$；修正前實際填的是最高 4 個，$expval(H) = +6.416$，正好是負值。],
)

= 📄 對學姊論文的影響

== 一個精確的等價關係

子格點變換 $C: b_(j s) arrow.r -b_(j s)$ 在佔據基底上是對角相位 $(-1)^(sum n_B)$。
它把所有躍遷項反號、卻不動任何 Hubbard 項 $n_arrow.t n_arrow.b$，所以

$
  C , H(v, w, U_A, U_B) , C^dagger = -H(v, w, -U_A, -U_B).
$

舊碼製備的最高本徵態正是 $C ket("GS")$；接下來每一步 Trotter 閘在 $C$ 下逐一對應到「$U$ 反號、時間反向」的閘。
因此*整條電路*（含 Trotter 誤差）與修正後程式在 $(-U_A, -U_B)$ 下的電路只差一個對角相位 $C$，
而 Berry 相位與極化都只由位元字串機率決定，$C$ 對它們完全不可見。

#figure(
  image("figures/fig_gauge.svg", width: 100%),
  caption: [子格點變換的圖示。上：原模型；下：把每個 B 格點的產生算符反號後，躍遷全部變號、Hubbard 項不變，整體乘 $-1$ 就是 $U arrow.r -U$ 的模型。],
)

#done(title: "數值驗證")[
  - $N = 3$（qiskit 精確態向量）：$max |p_"legacy"(U) - p_"fixed"(-U)| <= 3 times 10^(-16)$，四組 PBC／OBC、$Delta U = 0.3 dash 3$ 的參數皆然。
  - $N = 6$ 論文參數（CUDA-Q fp64）：Fig. 4 全部七個 $Delta U$ 的極化剖面，兩者差 $<= 7 times 10^(-10)$。
]

== 逐圖對照（$N = 6$、$T = 1$、$L = 40$，論文參數重算）

#btable(
  (auto, auto, 1fr),
  ([*圖*], [*是否受影響*], [*原因*]),
  [Fig. 2 Trotter 收斂], [否], [自我參照的 fidelity；最高能帶行列式態同樣是 $H_0$ 本徵態，兩版本差 $4 times 10^(-16)$],
  [Fig. 3、4 的 $Delta U = 0$ 曲線], [否], [$U_A = 0.01 arrow.r -0.01$ 對量子化的 Berry 相位與邊緣極化無可見影響],
  [Fig. 3 的 $Delta U eq.not 0$ 曲線], [*是*], [偏離量子化值的方向相反],
  [Fig. 4 的 $Delta U eq.not 0$ 曲線], [*是*], [體內偏移符號相反、兩端邊緣角色互換],
)

#figure(
  image("figures/_page_9_Figure_0.jpeg", width: 62%),
  caption: [論文原圖 Fig. 3：上圖 $gamma$ 對 $w$、下圖 $w$ 固定時 $gamma$ 對 $Delta U$。$Delta U eq.not 0$ 時曲線一律*往下*偏離量子化值。],
)

#figure(
  image("figures/paper_fig3_top.svg", width: 100%),
  caption: [以論文參數重算 Fig. 3 上圖。左（修正前）重現論文的往下偏離；右（修正後）偏離方向相反。$Delta U = 0$ 的階梯兩者相同。],
)

#figure(
  image("figures/paper_fig3_bottom.svg", width: 100%),
  caption: [以論文參數重算 Fig. 3 下圖。虛線橘為修正前、實線藍為修正後，叉號是「修正後但 $U arrow.r -U$」，與修正前逐點重合。論文曲線往下偏到 $-0.49pi$，正確物理往上偏到 $+0.35pi$。],
)

#btable(
  (auto, auto, auto),
  ([*$Delta U$（$w = 2$）*], [*#before（論文）$gamma \/ pi$*], [*#after $gamma \/ pi$*]),
  [0.1], [$-0.015$], [$+0.015$],
  [0.3], [$-0.046$], [$+0.044$],
  [1], [$-0.158$], [$+0.139$],
  [3], [$-0.490$], [$+0.349$],
)

#figure(
  image("figures/_page_9_Figure_10.jpeg", width: 62%),
  caption: [論文原圖 Fig. 4：上圖極化剖面（體內偏移為*負*）、下圖第一晶胞極化隨 $Delta U$ *下降*到 0.73。],
)

#figure(
  image("figures/paper_fig4_top.svg", width: 100%),
  caption: [以論文參數重算 Fig. 4 上圖。修正前（橘）重現論文：體內偏移為負、cell 1 邊緣極化減弱；修正後（藍）體內偏移為正、增強的是 cell 1、減弱的是 cell 6。],
)

#figure(
  image("figures/paper_fig4_bottom.svg", width: 72%),
  caption: [以論文參數重算 Fig. 4 下圖。論文曲線（橘）下降到 0.72，正確物理（藍）上升到 1.23；叉號再次確認修正前 $equiv$ 修正後的 $U arrow.r -U$。],
)

#btable(
  (auto, 1fr, 1fr),
  ([*$Delta U$*], [*#before（論文）cell 1 … cell 6*], [*#after*]),
  [0], [0.99, 0.01, 0.00, 0.00, −0.01, −0.99], [相同],
  [0.1], [0.976, −0.016, −0.026, −0.027, −0.037, −1.004], [1.003, 0.036, 0.026, 0.026, 0.016, −0.976],
  [1], [0.854, −0.258, −0.269, −0.269, −0.279, −1.130], [1.117, 0.261, 0.251, 0.251, 0.242, −0.858],
  [2], [0.721, −0.521, −0.532, −0.532, −0.543, −1.266], [1.226, 0.477, 0.468, 0.468, 0.459, −0.744],
)

#warn(title: "物理上一眼可辨的破綻")[
  $Delta U = U_B - U_A > 0$ 表示 B 格點排斥較強；超過半填充（14 電子）時電子應往 A 移動，體內 $expval(n_A - n_B)$ 應為*正*。
  論文圖中體內偏移為負，正是 $U$ 反號的徵兆，審稿人用平均場直覺即可看出。
]

== 需不需要修正

需要。不是因為結論被推翻，而是圖與文字描述的是與標示不同的哈密頓量：

+ 圖例標 $U_A = 0.01$、$Delta U > 0$（排斥），實際對應 $U_A = -0.01$、$Delta U < 0$（吸引）。
+ Fig. 4 體內偏移的符號、Fig. 3 偏離量子化值的方向，都與排斥型物理相反。
+ Sec. III-A 與 Eq. (15) 宣稱初態是基態並演化到目標基態；實際模擬追蹤的是扇區內的最高本徵態。方法敘述是對的，是實作偏離了敘述。

不需要修正：Fig. 1（電路）、Fig. 2、$Delta U = 0$ 的曲線、Sec. II–IV 全部解析推導，以及「弱不平衡下拓樸特徵穩健、強不平衡下破壞」的主要結論。

= 💬 初步討論方向

+ *修正範圍*：用修正後程式以相同參數重跑 `aqs berry` 與 `aqs polarization` 換掉 Fig. 3、Fig. 4，並改寫 Sec. V-B、V-C 描述偏離方向的句子（$gamma$ 往上偏；A 邊緣極化增強、B 邊緣極化減弱；體內偏移為正）。
+ *形式取決於發表狀態*：審稿中或 camera-ready 前直接換圖改文；已刊出則需與共同作者討論 erratum 並更新 arXiv 版本。建議修訂稿註明原圖對應 $U arrow.r -U$ 的模型。
+ *不建議的捷徑*：只把圖例改成 $-Delta U$、$U_A = -0.01$。文字中的排斥型敘述與式 (25) 的自旋翻轉安全條件仍會與圖不符。
+ *附帶議題*：程式預設 `mode="up_spin"` 只用自旋上區塊算 $z_N$，論文 Eq. (46) 則兩個自旋相加；重算數字顯示論文用的是前者（total 模式的偏離量恰為兩倍）。重跑前應先確認定義。
+ *`QCE26/` 目錄*：`origin/main` 上的會議程式碼仍是有缺陷的版本；用它新產生的 $Delta U eq.not 0$ 結果會重現同樣的符號問題。

= 📚 重現方式與相關檔案

#info(title: "重跑本報告的圖與數字")[
  ```bash
  # 本報告的圖（電路圖、矩陣、氫分子、能階、規範變換示意）
  uv run --no-sync --with pylatexenc python docs/defect1/make_figures.py
  # 論文參數 N=6 重算（GPU，約 18 分鐘）與繪圖
  uv run --no-sync python docs/defect_report/paper_impact_recompute.py
  uv run --no-sync python docs/defect_report/paper_impact_plot.py
  # 編譯本文件
  typst compile docs/defect1/defect1.typ
  ```
]

- `docs/defect_report/defects_before_after.ipynb`：三個缺陷的完整可執行對照，第 2 節為本缺陷。
- `docs/defect_report/slater_determinant_circuit_example.ipynb`：OpenFermion 演算法的獨立範例（Cirq），含以行列式公式驗證振幅。
- `docs/defect_report/paper_impact.md`、`paper_impact.json`：論文影響分析與全部重算數據。
- `tests/test_state_preparation.py`：守門測試；把 `cry(-2θ)` 放回去會有 4/6 個測試變紅。
