# AQS 目前進度

日期：2026-09-14。分支：`fix/core-silent-convention-bugs`。
此次將既有工作依需求納入不同本地 commit，沒有重寫舊 commit，也沒有 push。

## 已完成與提交分類

| 需求 | 內容 | Commit |
|---|---|---|
| 重建規格與物理規則 | spinful 論文、spinless 原始規格、圖片與七份細則 | `0417abe` |
| Spinless 模型與驗證 | 模型、態製備、Trotter、觀測量、sweep／benchmark、H_eff 與固定扇區等測試；套件納入 aqs.models | `32888fc` |
| 物理特徵測試 | 能帶、Berry 跳躍、邊緣極化、交互作用與收斂特徵 | `a4f2ae4` |
| 技術文件 | 1063 行 Typst 架構與物理演算法文件 | `42efa34` |
| Claude Code／Codex 共用與交接 | AGENTS.md、CLAUDE.md、core-rules.md、專案記憶備份、遷移報告及本進度 | 包含本文件的 commit |

先前分支已包含 gate primitive、Givens 態製備、CUDA-Q 與 OpenFermion
端序修正，以及 uv／dev dependency 設定；詳見 `docs/experience.md` 與 Git 歷史。

## 本次驗證

- `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest tests/ -q`：
  沙箱外 **156 passed、1 skipped**。第一次沙箱內為 152 passed、1 skipped、
  4 failed，四項失敗均為 CUDA 裝置不可見；沒有為此修改測試或降低斷言。
- Qiskit 與 CUDA-Q `aqs selftest` 均 PASS，fidelity 顯示 1.00000000。
  CUDA-Q selftest 在沙箱內可 fallback，不能單憑此項宣稱使用 NVIDIA GPU；
  GPU 存取相關測試由沙箱外完整測試涵蓋。
- 在暫存副本使用 setuptools 建置 wheel，確認包含 `aqs/models/spinless.py`，
  解包後由該 wheel 內容成功 import。依賴沿用既有 venv，未驗證全新機器依賴安裝。
- `typst compile docs/codebase_architecture.typ docs/codebase_architecture.pdf` 成功。
  PDF 依既有 .gitignore 忽略，只提交 Typst 原始檔；未逐頁複核物理敘述。
- 共用入口、七份细則、記憶原文備份已做靜態檢查。

## 尚未完成／下一步

1. **完整架構重構**：spinful 仍在 `core.py`；`mapping.py`、`circuits/`、
   `backends/` 目錄化尚未完成。目標架構表不能視為已實作清單。
2. **Spinless 執行整合**：目前 driver 直接用 Qiskit Statevector；尚未接入
   backend 選擇、共用 experiments／CLI。保留 CUDA-Q 是需求，不代表此 driver 已支援。
3. **完整參考圖重現**：測試涵蓋定性特徵、H_eff、L 收斂與固定扇區參考，
   但本次未完成全部 N=6、T=1、L=40 圖表產出、逐圖比對或重做 mutation 驗證。
   未重新量測覆蓋率，不沿用歷史 99% 作本次結論。
4. **搬移驗收**：在新 Claude Code／Codex 工作階段確認實際規則載入；
   沒有專案級 Skills，未來需要時才新增共用 SKILL.md。
5. **歷史文件與發行物**：細則與 experience 中的舊分支／環境資訊仍須按現況核對；
   `dist/` 的既有 2.0.0 檔案未更新，不能代表目前 src。
6. **遠端備份**：本次僅完成本地 commits，尚未 push 或合併其他分支。
