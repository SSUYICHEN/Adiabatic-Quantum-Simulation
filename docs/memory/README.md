# 專案記憶索引

- [重建決策原始備份](aqs-reconstruction-governance.md)：2026-08-11 的 Claude
  專案記憶，於 2026-09-14 原文備份。來源為
  `~/.claude/projects/-mnt-sata-b11202015-Adiabatic-Quantum-Simulation/memory/aqs-reconstruction-governance.md`。
  保留兩模型、兩套規格權威、雙後端與 gate convention 注意事項。
- [除錯經驗](../experience.md)：既有本地歷史紀錄，含端序、閘與態製備教訓。
- [遷移檢查](../agent-migration-audit.md)：本次觀察、修正與待處理事項。
- [缺陷前後對照](../defect_report/defects_before_after.ipynb)：2026-09-16 製作。三個會改變物理結果的缺陷
  （Givens 反號、CUDA-Q 位元序、OpenFermion 端序）以可執行實例並排修正前後數字；
  同目錄 `.typ` 讀取 `results.json` 產生 PDF 摘要。合併 PR 前給協作者閱讀。
- [論文影響分析](../defect_report/paper_impact.md)：2026-09-16。缺陷 A 精確等價於 U 反號；論文 Fig. 3、4 的 ΔU≠0 曲線需重製。
  全部圖形前後對照見 [paper_all_figures.md](../defect_report/paper_all_figures.md)，
  給協作者的完整報告見 [defect1.typ](../defect1/defect1.typ)。

原始備份描述的是當時狀態；「尚未重構、src 全為平面」已不完全適用。
2026-09-14 已提交 `src/aqs/models/spinless.py`；重構仍未全部完成。
分支、硬體、磁碟、測試數與覆蓋率必須重新確認，不能視歷史值為目前證據。

後續重要決策請記錄日期、原因、涉及檔案及驗證結果，並加入此索引。
共用工作規則維護於根目錄 `core-rules.md`；不要在記憶內維護另一套規則。

目前進度與本次提交分類見 [project-progress.md](../project-progress.md)。
