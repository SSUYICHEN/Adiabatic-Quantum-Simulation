# 交接筆記：AQS 專案的除錯經驗

> 寫給接手的人（含未來的我）。這份文件記錄的不是「程式怎麼跑」——那在 `README.md`
> 與 `docs/codebase_architecture.typ` 裡——而是**這個 repo 會怎麼騙你**，以及哪些驗證
> 手法真的抓得到問題。內含我實際犯過的錯，請當作地雷圖看。

---

## 0. 一句話總結

**這個 repo 的缺陷全是「靜默的慣例不符」：不會拋例外、不會有測試失敗、輸出看起來
完全合理，只是物理是錯的。** 已找到並修正四個，全部屬於同一類：

> 某個函式實際產生的東西 $\neq$ 它宣稱要產生的東西，而呼叫端可能剛好抵銷、也可能沒有。

因此**不要用「跑起來沒錯」當作正確性證據**。每一項都要對照獨立計算的參考值。

---

## 1. 位元序總表（全部實測，勿憑記憶）

這是本專案最容易出錯的單一議題，已咬過三次。

### 索引型（振幅陣列的整數索引）

| 介面 | 慣例 | X 作用在 q0（3 qubits） |
|---|---|---|
| qiskit `Statevector` | 🟢 Little（q → bit q） | index `1` |
| qiskit `Operator` / kron | 🟢 Little（q0 是**最右**因子） | `kron(I, X)` |
| CUDA-Q `get_state()` | 🟢 Little | index `1` |
| OpenFermion `get_sparse_operator` | 🔴 **Big**（q → bit n−1−q） | mode 0 → bit `2` |
| `aqs.observables._bit(idx,pos)` | 🟢 Little | `(idx >> pos) & 1` |

`FermionOperator` 與 `QubitOperator` 皆為大端序。

### 字串型（量測位元字串）

| 介面 | qubit 0 的位置 | X@q0 |
|---|---|---|
| qiskit `counts` / `to_dict()` | **最右** | `'001'` |
| CUDA-Q `sample()` | **最左** | `'100'` |

### 三個必須記住的推論

1. **CUDA-Q 自己內部就不一致**：`get_state()` 小端序、`sample()` 字串卻是 q0 在最左。
   原始碼曾有註解「CUDA-Q 以 qubit 0 為 MSB」——對字串成立、對態向量不成立。
   這句半對的話造成了 `selftest` fidelity 0.181。

2. **qiskit 與 CUDA-Q 字串方向相反，但索引方向相同。**
   若從「字串長得不一樣」推論「索引也要反轉」，就會寫出多餘的 `q[n-1-i]`。
   **對齊態向量請用恆等映射。**

3. **OpenFermion 是唯一的大端序，而且實部項不會暴露它。**
   兩模式系統下它交換 mode 0/1 並**共軛**帶 $\mathrm{Im}(t)$ 的非對角元；
   $XX+YY$ 在交換下不變，只有複數躍遷才會現形。

專案內的轉換工具：`hamiltonian._to_little_endian(H, n)`（置換相似變換，保譜）。

---

## 2. 真正有效的驗證技術

### 2.1 有效哈密頓量往返（抓系統性誤差）

$$H_\text{eff} = i\log(U)/\mathrm{d}t$$

與真值逐元比對。**這是唯一能抓到「虛部被共軛」的方法** ——
那種錯誤是系統性的，**不隨 $\mathrm{d}t \to 0$ 縮小**，普通收斂測試看不到。

### 2.2 用「是否隨 L 收斂」區分兩類誤差

| 現象 | 判讀 |
|---|---|
| 保真度隨 $L$ 增大而改善 | Trotter 離散化誤差，正常 |
| 保真度**停在某個平台不動** | **系統性錯誤**，去找慣例問題 |

我就是靠「複數躍遷停在 0.80 且不隨 $L$ 改善」才發現 `_G_SIGN` 設錯。

### 2.3 分離絕熱誤差與 Trotter 誤差

固定 $\mathrm{d}t = T/L$ 同時增大 $T$ 和 $L$，保真度會停在平台——**這不是絕熱性失效**，
是 Trotter 誤差為常數。必須固定一項、單獨掃描另一項。

### 2.4 參考值絕不能自我參照

`fidelity_scan` 以**製備態本身**為參考比較演化後的態，因此對「製備了哪一個行列式態」
**結構上就不敏感**——Givens 缺陷躲了很久就是因為它。

好的參考：
- 製備態能量 $=$ 最低填充能量（手徵對稱性讓錯誤答案恰好是 $-E$，斷言會**翻號**而非漂移）
- 對照 OpenFermion 直接建構的 Slater 行列式
- 固定粒子數子空間內的精確對角化

### 2.5 測試輸入必須「不對稱」

**內附的兩個範例密度均勻（每格 0.5），位元反轉在其上是恆等操作**——它們在結構上
就無法偵測反轉類錯誤。同理，反轉對稱的電路分辨不出兩種 qubit 慣例。

診斷手法：在某個 mode 設深位能井，看密度從哪個 qubit 讀出來。

### 2.6 反向驗證每一個守門測試

寫完測試後，**把舊的錯誤程式碼放回去，確認測試真的會失敗**。
本分支四組守門測試全部這樣驗過（舊碼 4/6、4/5、3/5、2/7 失敗）。
沒反向驗過的測試，等於沒測。

---

## 3. 我犯過的錯（請勿重蹈）

### 3.1 用大端序矩陣「證實」了錯誤的符號 ⚠️ 最嚴重

我曾信誓旦旦說 `core.py` 的 G 閘符號有錯，並據此設 `_G_SIGN = +1`。
**那個結論本身是錯的** —— 驗證時直接用了 OpenFermion 的大端序矩陣沒轉換，
兩模式下這會交換 mode 0/1 並悄悄共軛待檢查的非對角元。

**教訓：驗證程式本身也會有位元序錯誤。** 任何涉及 $\mathrm{Im}$ 的推導，
一律先 `_to_little_endian` 再比對。

### 3.2 「bit-identical」講得太滿

我在 commit message 與 MR 裡寫 `fidelity_scan` 修正前後「逐位元相同」。
實測最大差 $6.7\times10^{-16}$ —— 以 `UnitaryGate` 取代 `rz`/`cx`/`cry` 分解會走不同
的浮點捨入路徑。物理無差異，但**不是**逐位元相同。

**教訓：宣稱 bit-identical 前，用 `==` 比對原始浮點，不要用 `:.6f` 格式化後目視。**

### 3.3 `git add docs/` 掃進別的分支的檔案（犯了兩次）

`docs/codebase_architecture.typ` 是另一分支的未追蹤檔案，會跟著 checkout 移動。
兩次都被 `git add docs/` 掃進 commit，事後 amend 移除。

**教訓：一律 `git add <明確檔名>`，不要 `git add <目錄>`。**

### 3.4 `git stash -u` 把 4.2 GB 的 `.venv` 收走

`-u` 會 stash 未追蹤檔案，包含整個虛擬環境。雖然 pop 後完好，但極其危險且緩慢。

**教訓：只 stash 明確的追蹤檔；要跨 revision 比對請用 `git worktree`：**

```bash
git worktree add -q --detach /path/to/wt main
PYTHONPATH=/path/to/wt/src .venv/bin/python script.py out.json
git worktree remove --force /path/to/wt
```

### 3.5 stash 堆疊取錯層

之前有一個殘留的 `core.py` stash，我 `git stash pop` 時取到它，
結果把 gate 正規化的修改 pop 到了文件分支上。

**教訓：`git stash pop` 前先 `git stash list` 確認。**

---

## 4. 物理慣例陷阱

### 4.1 Berry 相位帶有填充／宇稱位移

$\hat{X} = \sum_j (j+1) n_j$ 定義的**絕對** Berry 相位帶常數位移 $\pi(N+1)$。
**偶數 $N$ 時平庸相位在 $\gamma = \pi$、拓樸相位在 $0$**，與直覺標記相反。

**有物理意義的是跨相變時 $\pi$ 的跳躍，不是絕對值。**（已對 $N=4..7$ 驗證。）

### 4.2 PBC 環繞鍵的不對稱性

- **躍遷**項跨環繞鍵帶 Jordan–Wigner 弦 → 角度需乘宇稱因子 $-(-1)^{N_f}$
- **交互作用**項 $n_i n_j$ 是對角的 → **完全不帶弦，不可加宇稱修正**

這個不對稱極易寫錯。

### 4.3 PBC 下極化恆為零

沒有邊緣。邊緣響應必須用 OBC。

### 4.4 精確基態要在固定粒子數子空間取

電路守恆粒子數，所以參考態必須是**同一扇區**的基態。
排斥型模型的全域基態通常在不同填充（$N=3$ 交互作用範例落在 2 個粒子而非 3 個）。

### 4.5 閘原語的角度慣例

`create_R_gate` / `create_G_gate` / `create_CP_gate` 現在**與 docstring 一致**
（本分支修正）。若在 `main` 或另一分支上工作，它們差因子 2 或負號，呼叫端以
`-2.0*tau*w` 抵銷。**兩種慣例不可混用**。

---

## 5. 目前分支狀態

```
main                                   0fe9d1d   （未動）
fix/core-silent-convention-bugs        b8bc200   9 commits ← 目前所在
fix/cudaq-bit-order-and-uv-migration   c7057d1   7 commits
```

### `fix/core-silent-convention-bugs`（四項缺陷 + 測試）

| 修正 | 對既有結果 |
|---|---|
| 閘原語與 docstring 不符 | 無（逐位元相同） |
| Givens 製備出最高能帶 | **極化顯著改變，需重繪** |
| CUDA-Q 態向量反轉 | 所有 GPU 電路結果先前皆錯 |
| OpenFermion 大端序 | 內附範例不變；其他非對稱哈密頓量先前皆錯 |

測試 112 個、覆蓋率 99%。MR 內文在 `docs/MR_core_silent_convention_bugs.md`。

### `fix/cudaq-bit-order-and-uv-migration`（新功能）

uv 遷移、`spinless.py`（$2N$ 無自旋 SSH + 最近鄰交互作用）、實驗執行器、
相圖、CLI 子指令。也**獨立修過** CUDA-Q 與 OpenFermion 位元序、Givens。

### ⚠️ 兩分支合併會靜默衝突

`spinless.py` 用**舊**慣例呼叫閘原語：

```python
th_R = -2.0 * dt * t_re
th_G = _G_SIGN * 2.0 * dt * t_im     # _G_SIGN = -1
create_CP_gate(-dt * V * lam)
```

若兩分支各自落地，`spinless.py` 會靜默錯掉（實測絕熱收斂 0.99988 → **0.930398**）。
**git 不會標記衝突** —— 改的是不同檔案，只有語義衝突。後合併者必須改成：

```python
th_R = -dt * t_re
th_G = dt * t_im                     # _G_SIGN 不再需要
create_CP_gate(dt * V * lam)
```

建議：先合 `fix/core-silent-convention-bugs`，再 rebase 另一分支並重跑其驗證
（態製備保真度 1.0、絕熱收斂 $\geq 0.9997$）。

### 待辦

- [ ] 兩分支皆未推送（無 `gh` CLI，需手動開 PR）
- [ ] 合併後重新生成所有極化圖表
- [ ] `docs/codebase_architecture.typ`（1063 行繁中技術文件）在另一分支上仍**未提交**
- [ ] `dist/` 內是 2.0.0 的舊 wheel，與 `src/` 可能不同步
- [ ] CUDA-Q `nvidia` target 預設 fp32；`--fp64` 未接入 CLI

---

## 6. 環境須知

| 項目 | 狀態 |
|---|---|
| GPU | RTX 4090，driver 560.35.05 = **CUDA 12.6** |
| GPU 套件 | 必須用 **cu12** wheel。`requirements-gpu.txt` 原本硬釘 CUDA 13，會壞掉 |
| Python | venv 是 3.10；`cuda-quantum-cu*` $\geq$ 0.13 需要 $\geq$ 3.11 |
| 磁碟 | **很緊**。`/` 97% used、`/mnt/sata` 100% used。CUDA wheel 是 GB 級，動手前先 `df -h` |
| `uv` | 0.9.11。`uv sync` **不帶 extra 會移除 ~4 GB GPU 套件**，GPU 機器上一律 `--extra cu12` |
| `gh` | **未安裝**，無法自動開 PR |
| `typst` | 0.15.0（snap）。**讀不到 `/tmp`**，暫存檔要放在專案目錄內 |
| 字型 | 有 Noto CJK TC。**沒有** `New Computer Modern Sans`（只有 serif + math） |
| pytest / coverage | 我手動裝進 venv；`uv sync` 會移除，見 `requirements-dev.txt` |

### Typst 語法注意（0.15）

- `angle.l` / `angle.r` **不存在** → 直接用 Unicode `⟨` `⟩`，ket 用 `lr(|00⟩)`
- `bra` / `ket` 函式**不存在**
- `[tool.uv] default-extras` **不是有效欄位**（uv 0.9.11 會警告並忽略）

---

## 7. 指令速查

```bash
PY=.venv/bin/python                      # 直接用 venv，避免 uv 重新解析 pyproject

# 測試與覆蓋率
$PY -m pytest tests/ -q
$PY -m coverage run --source=src/aqs -m pytest tests/ && $PY -m coverage report -m

# 兩後端自我測試（cudaq 那項是位元序的守門員）
$PY -c "from aqs.backends import selftest; print(selftest('cudaq'))"

# 跨 revision 一致性對照
$PY tests/verify_consistency_vs_main.py out.json
$PY tests/verify_consistency_vs_main.py --compare before.json after.json

# 其他一次性驗證腳本
$PY tests/verify_gate_norm_equivalence.py out.json    # 閘正規化：逐位元相同
$PY tests/verify_givens_impact.py out.json            # Givens：量化影響

# GPU 環境（另一分支才有 uv 設定）
uv sync --extra cu12
uv run aqs selftest --backend cudaq
```

---

## 8. 給接手者的三條建議

1. **先跑 `$PY -m pytest tests/ -q`。** 112 個測試涵蓋了四個已知缺陷的所有慣例。
   若有紅燈，先弄清楚是不是又踩到位元序。

2. **改任何閘或態製備前，先讀 §1 與 §2。** 這個 repo 的錯誤不會自己浮現，
   你必須主動去對照獨立參考值。

3. **不要相信「看起來合理」的數字。** 極化 $[0.926, 0.0, -0.926]$ 看起來很漂亮，
   但它來自能量最高的行列式態，是錯的。正確答案是 $[0.835, 0.0, -0.835]$ ——
   一樣漂亮。**只有對照精確對角化才分得出來。**
