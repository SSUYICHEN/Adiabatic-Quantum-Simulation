// =====================================================================
//  AQS - 三個會改變物理結果的靜默缺陷：修正前後對照（PDF 摘要）
//  數字全部來自 results.json，由 defects_before_after.ipynb 執行時寫出。
//  以 typst compile docs/defect_report/defects_before_after.typ 編譯
// =====================================================================

#let r = json("results.json")

#set document(
  title: "AQS 靜默缺陷修正前後對照",
  author: "Adiabatic-Quantum-Simulation",
)

#set page(
  paper: "a4",
  margin: (x: 2.2cm, y: 2.4cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 8pt, fill: rgb("#666666"))
      AQS · 靜默缺陷修正前後對照
      #h(1fr)
      修正後 #raw(r.meta.after_branch) \@ #raw(r.meta.after_sha)
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

#show table: set block(breakable: false)
#show grid: set block(breakable: false)

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
  stroke: (x, y) => if y == 0 { (bottom: 0.8pt + rgb("#0b3d61")) } else { (bottom: 0.3pt + rgb("#dddddd")) },
  fill: (x, y) => if y == 0 { rgb("#eef3f8") } else { none },
  inset: (x: 6pt, y: 5pt),
  ..cells,
)

// 固定小數位數的數字格式（帶正負號）
#let f(x, d: 4) = {
  let a = calc.abs(x)
  let scale = calc.pow(10, d)
  let n = calc.round(a * scale)
  let ip = calc.floor(n / scale)
  let frac = int(calc.round(n - ip * scale))
  let fs = str(frac)
  while fs.len() < d { fs = "0" + fs }
  let sign = if x < 0 and n > 0 { "−" } else { "+" }
  sign + str(int(ip)) + "." + fs
}
#let fl(xs, d: 4) = "[" + xs.map(x => f(x, d: d)).join(", ") + "]"
#let before = text(fill: rgb("#c2410c"))[修正前]
#let after = text(fill: rgb("#1d5fb3"))[修正後]

// =====================================================================
#align(center)[
  #text(size: 20pt, weight: "bold", fill: rgb("#0b3d61"), font: ("Noto Sans", "Noto Sans CJK TC"))[
    三個會改變物理結果的靜默缺陷
  ]
  #v(0.3em)
  #text(size: 13pt, fill: rgb("#14507d"))[修正前後對照 · PDF 摘要]
  #v(0.6em)
  #text(size: 9.5pt, fill: rgb("#555555"))[
    修正前 = `main` \@ #raw(r.meta.before_sha)（與 `QCE26/src/aqs` 逐位元相同）
    #h(1.5em)
    修正後 = #raw(r.meta.after_branch) \@ #raw(r.meta.after_sha)
  ]
]
#v(1em)

#note("這份摘要怎麼用")[
  完整推導、可執行的舊碼與每一格的輸出都在同目錄的 `defects_before_after.ipynb`；
  本 PDF 只摘錄結論與關鍵數字。所有數字由 notebook 執行時寫入 `results.json`，
  本文件直接讀取，不另行手抄。三個缺陷的共同特徵：*不拋例外、不產生 NaN、輸出看起來合理，只是物理錯了。*
]

= 執行摘要

#tbl(
  (auto, 1fr, 1fr, 1fr),
  (left, left, left, left),
  [*代號*], [*缺陷*], [*修正前的錯法*], [*受影響的功能*],
  [*A*], [Slater 態製備的 Givens 旋轉角度反號], [製備出*能量最高*的行列式態（$+E$ 而非 $-E$）],
    [`aqs polarization`（極化翻號）、$Delta U eq.not 0$ 的 `aqs berry`],
  [*B*], [CUDA-Q 後端 qubit 映射 $i arrow.r n-1-i$], [態向量位元反轉，密度剖面左右顛倒],
    [*所有* `--backend cudaq` 的 GPU 結果],
  [*C*], [`load_hamiltonian` 未把 OpenFermion 大端序轉成小端序], [非對稱哈密頓量的觀測量讀在錯的 qubit；複數躍遷虛部被共軛],
    [`aqs measure` 的 `fermion_operator` / `qubit_operator` 格式],
)

#danger("合併後必須重新產生的結果")[
  + 所有極化圖（缺陷 A；$Delta U = 3$ 時中間晶胞由 #f(r.A.polarization.at(3).before.at(1), d: 3) 變為 #f(r.A.polarization.at(3).after.at(1), d: 3)，*正負號翻轉*）。
  + 所有以 `--backend cudaq` 得到的數字（缺陷 B；`selftest` 修正前 fidelity #f(r.B.selftest_before, d: 3)）。
  + 任何以 `aqs measure` 讀入非均勻密度哈密頓量的結果（缺陷 C）。內附兩個均勻範例不受影響。
]

// =====================================================================
= 缺陷 A：Givens 旋轉反號，製備出能量最高的行列式態

== 根因

態製備用 OpenFermion 的 `slater_determinant_preparation_circuit(Q)` 給出一串 Givens 旋轉 $(j, k, theta, phi)$。
舊碼把每個旋轉手工分解成 `rz / cx / cry / cx`，`cry` 的角度多了負號，實際做出的是 $G(-theta)$。
在手徵對稱的 SSH 鏈上，$theta arrow.r -theta$ 正好把最低軌域換成最高軌域。

#grid(columns: (1fr, 1fr), gutter: 10pt,
  [#before
  ```python
  def _givens_instruction(theta, phi):
      gv = QuantumCircuit(2)
      gv.rz(phi, 1); gv.rz(-phi, 0)
      gv.cx(0, 1)
      gv.cry(-2.0 * theta, 1, 0)   # 反號
      gv.cx(0, 1)
      return gv.to_instruction()
  ```],
  [#after
  ```python
  def _givens_instruction(theta, phi):
      c, s = np.cos(theta), np.sin(theta)
      e = np.exp(1j * phi)
      U = np.zeros((4, 4), dtype=complex)
      U[0,0] = 1;  U[1,1] = c;  U[1,2] = s
      U[2,1] = -s*e;  U[2,2] = c*e;  U[3,3] = e
      return UnitaryGate(U, label="Givens")
  ```],
)

== 玩具例：兩個模式、一個粒子，$H = -t(a^dagger b + b^dagger a)$

基態是對稱組合 $(|01⟩ + |10⟩) \/ sqrt(2)$，能量 $-t$；反對稱組合的能量是 $+t$。

#tbl(
  (auto, auto, auto, auto, auto, auto),
  (left, right, right, right, right, right),
  [], [$|00⟩$], [$|01⟩$], [$|10⟩$], [$|11⟩$], [$⟨H⟩ \/ t$],
  [#before], ..r.A.toy.before.map(x => f(x)), [#f(r.A.toy.E_before)],
  [#after],  ..r.A.toy.after.map(x => f(x)),  [#f(r.A.toy.E_after)],
)

兩個態的*密度完全相同*（每個模式各 0.5），只看佔據數看不出差別；差在振幅的相對符號，也就是能量。

== 真實例：製備態能量對照最低／最高填充

判準：正確的 Slater 行列式，其 $⟨H⟩$ 必須*精確等於*被佔據的最低單粒子能階之和 $sum epsilon$。
手徵對稱讓最高填充能量恰為其負值，所以錯的版本不是差一點，而是整個符號相反。
另以 OpenFermion 直接建的參考行列式（轉成小端序後）算 fidelity。

#tbl(
  (1fr, auto, auto, auto, auto, auto),
  (left, right, right, right, right, right),
  [*組態*], [*$E_"prep"$ 前*], [*$E_"prep"$ 後*], [*最低填充*], [*fid. 前*], [*fid. 後*],
  ..r.A.prep.map(p => (
    [#p.label], [#f(p.E_before)], [#f(p.E_after)], [#f(p.E_lowest)], [#f(p.fid_before, d: 3)], [#f(p.fid_after, d: 3)],
  )).flatten(),
)

== 對 `aqs polarization` 的影響

$N=3$、$v=0.5$、$w=1.5$、$U_A=1.0$、OBC、8 電子、$T=1$、$L=40$。初態錯了，絕熱演化就從錯的態出發，
終態的邊緣極化 $P_j = ⟨n_(A,j) - n_(B,j)⟩$ 整個改變。

#figure(image("figures/A_polarization.svg", width: 100%))

#tbl(
  (auto, 1fr, 1fr, auto),
  (right, left, left, right),
  [*$Delta U$*], [*#before $P_j$（cell 0, 1, 2）*], [*#after $P_j$*], [*最大變化*],
  ..r.A.polarization.map(p => (
    [#p.delta_U], [#fl(p.before)], [#fl(p.after)],
    [#f(calc.max(..p.before.zip(p.after).map(xy => calc.abs(xy.at(0) - xy.at(1)))))],
  )).flatten(),
)

#note("為何 `aqs fidelity` 沒有發現它")[
  `fidelity_scan` 拿演化後的態跟*製備態本身*比，屬自我參照；兩版本輸出最大差
  約 #str(r.A.fidelity_scan_max_diff)，只是浮點捨入。參考值必須來自程式之外。
]

// =====================================================================
= 缺陷 B：CUDA-Q 後端把 qubit 順序反過來

== 根因

舊碼註解寫「CUDA-Q 以 qubit 0 為最高有效位元」，於是把 qiskit qubit $i$ 接到 CUDA-Q qubit $n-1-i$。
這句話*只對 `sample()` 回傳的位元字串成立*；`get_state()` 回傳的態向量索引與 qiskit 一樣是小端序。
從字串慣例推論索引慣例，就多做了一次反轉。

#grid(columns: (1fr, 1fr), gutter: 10pt,
  [#before
  ```python
  def cq(qiskit_index):
      return q[n - 1 - qiskit_index]
  ```],
  [#after
  ```python
  def cq(qiskit_index):
      return q[qiskit_index]
  ```],
)

== 直接量測：X 閘作用在 3 個 qubit 上

#tbl(
  (auto, auto, auto, auto),
  (left, right, right, right),
  [*電路*], [*態向量非零索引（qiskit 與 CUDA-Q 相同）*], [*qiskit counts 字串*], [*CUDA-Q sample 字串*],
  ..r.B.x_table.map(x => (
    [X on qubit #x.q], [#x.index], [#raw(x.qiskit_string)], [#raw(x.cudaq_string)],
  )).flatten(),
)

索引相同、字串反向。對齊態向量必須用恆等映射。

== 舊映射對態向量做了什麼

把 qubit $i$ 接到 $n-1-i$ 再原樣回傳，等價於對正確態向量做位元反轉置換 $P$；不需要 GPU 就能精確重現。
以 `aqs selftest` 的 4-qubit 電路為例：

#tbl(
  (1fr, auto),
  (left, right),
  [], [*selftest fidelity*],
  [#before（$q[n-1-i]$，等價於位元反轉）], [#f(r.B.selftest_before, d: 8)],
  [#after（恆等映射，qiskit 模擬）], [#f(r.B.selftest_after_sim, d: 8)],
  [#after（真實 CUDA-Q 後端）], [#if r.B.selftest_after_gpu == none [（執行時無 CUDA-Q）] else [#f(r.B.selftest_after_gpu, d: 8)]],
)

== 對物理結果的影響：密度剖面整個顛倒

同一個 SSH–Hubbard 終態（$N=3$、OBC、8 電子、$Delta U=1$，12 qubit），舊映射回傳的態向量算出的
每個 qubit 佔據數左右顛倒：自旋區塊互換、格點順序反轉。

#figure(image("figures/B_density.svg", width: 100%))

#tbl(
  (auto, 1fr),
  (left, left),
  [], [*極化 $P_j$（cell 0, 1, 2）*],
  [#before], [#fl(r.B.density.pol_before)],
  [#after],  [#fl(r.B.density.pol_after)],
)

#danger("影響範圍")[
  凡是修正前用 `--backend cudaq` 跑出的數字都來自顛倒的態向量，*全部 GPU 結果都需要重跑*。
  對稱的測試態看不出這個問題，所以 `aqs selftest` 現在刻意用不對稱電路。
]

// =====================================================================
= 缺陷 C：OpenFermion 大端序沒有轉成 qiskit 小端序

== 根因

`aqs measure` 讀 `fermion_operator` / `qubit_operator` 格式時經過 OpenFermion 的 `get_sparse_operator`，
它把模式 $q$ 放在第 $n-1-q$ 個 bit（大端序）；`observables.PROPERTIES` 全部按 qiskit 慣例讀第 $q$ 個 bit。
舊 `load_hamiltonian` 直接把矩陣交出去，基態的每個 bit 都被反著讀。修正後多一步置換相似變換
$H arrow.r P H P^top$（`_to_little_endian`），只換基底標籤，本徵譜不變。

== 診斷例：把粒子關在模式 0

四模式鏈，在模式 0 放深度 $-8$ 的位能井，鍵結 $-0.2$。基態粒子必在模式 0，密度應在 *qubit 0* 最大。

#figure(image("figures/C_trap_density.svg", width: 78%))

#tbl(
  (auto, auto, auto, auto, auto, auto, 1fr),
  (left, right, right, right, right, left, left),
  [], [$⟨n_0⟩$], [$⟨n_1⟩$], [$⟨n_2⟩$], [$⟨n_3⟩$], [*密度最大*], [*極化 $P_j$（cell 0, 1）*],
  [#before], ..r.C.trap.before_density.map(x => f(x)), [qubit 3 ✗], [#fl(r.C.trap.before_pol)],
  [#after],  ..r.C.trap.after_density.map(x => f(x)),  [qubit 0 ✓], [#fl(r.C.trap.after_pol)],
)

兩版本的最低三個本徵值完全相同（#fl(r.C.trap.eig_after, d: 3)），只看能量察覺不到；
密度剖面卻左右顛倒，極化的兩個晶胞順序反轉且符號相反。

== 為何內附範例看不出來

`examples/ssh_spinless_N3_*.json` 是半填滿、密度處處 0.5 的均勻鏈。位元反轉對這種對稱分布是恆等操作。

#tbl(
  (1fr, auto, auto, auto, auto),
  (left, right, right, right, right),
  [*範例*], [*$|z_N|$ 前*], [*$|z_N|$ 後*], [*max$|P_j|$ 前*], [*max$|P_j|$ 後*],
  ..r.C.examples.map(e => (
    [#raw(e.name)], [#f(e.z_before, d: 6)], [#f(e.z_after, d: 6)], [#f(e.max_pol_before, d: 6)], [#f(e.max_pol_after, d: 6)],
  )).flatten(),
)

== 複數躍遷：實部看不到、虛部被共軛

兩模式的大端序 $arrow.l.r$ 小端序轉換等於交換兩個模式的標籤，會把非對角元 $-t$ 變成 $-t^*$。
實數躍遷毫無差別；$t$ 一有虛部，舊碼交出的哈密頓量虛部符號相反，本徵值卻仍相同。

#let cx(z) = f(z.at(0)) + " " + f(z.at(1)) + "i"
#tbl(
  (1fr, auto, auto),
  (left, right, right),
  [$t = #cx(r.C.complex.t)$], [*$H[|01⟩, |10⟩]$*], [*本徵值*],
  [#before（OpenFermion 原樣）], [#cx(r.C.complex.before_offdiag)], [#fl(r.C.complex.eig)],
  [#after（`_to_little_endian`）], [#cx(r.C.complex.after_offdiag)], [#fl(r.C.complex.eig)],
  [應得 $-t$], [#cx((-r.C.complex.t.at(0), -r.C.complex.t.at(1)))], [],
)

// =====================================================================
= 總結與守門測試

#tbl(
  (auto, 1fr, 1fr),
  (left, left, left),
  [], [*怎麼確認自己沒踩到*], [*本分支的守門測試（反向驗證過）*],
  [*A*], [製備態 $⟨H⟩$ 必須等於最低填充 $sum epsilon$], [`tests/test_state_preparation.py`：放回 `cry(-2θ)` 有 4/6 個測試變紅],
  [*B*], [`aqs selftest --backend cudaq` fidelity 必須為 1.0（不對稱電路）], [`tests/test_backend_endianness.py`：放回 `q[n-1-i]` selftest 掉到 0.181],
  [*C*], [陷阱模式測試：粒子關在模式 $m$ 要從 qubit $m$ 讀出], [`tests/test_hamiltonian_endianness.py`：拿掉 `_to_little_endian` 四個模式全部讀錯],
)

#danger("給 QCE26/ 的提醒")[
  `origin/main` 上 `QCE26/src/aqs/` 內的 `core.py`、`hamiltonian.py`、`backends.py` 與修正前版本逐位元相同，
  `cudaq_core.py` 亦使用 $n-1-i$ 映射。以該目錄產生的數字同時帶有 A、B、C 三個缺陷。
]

#ok("重跑本文件")[
  ```bash
  uv run --with jupyter jupyter nbconvert --to notebook --execute --inplace \
      docs/defect_report/defects_before_after.ipynb
  typst compile docs/defect_report/defects_before_after.typ
  ```
]
