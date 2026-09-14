# Claude Code → Codex 遷移檢查

檢查日期：2026-09-14。結論：適合遷移；本次已補上共同入口與專案記憶備份，
但仍需將未追蹤檔案納入版本控制，並驗證新工作階段的規則讀取。
本報告原始檢查是設定／文件檢查，沒有驗證物理結果，也沒有修改 Python 程式。
後續分類提交與驗證結果見 [目前進度](project-progress.md)；下文原始觀察保留作歷史紀錄。

## 六項檢查

| 項目 | 原始狀態 | 本次結果 |
|---|---|---|
| AGENTS.md | 專案內沒有 | 新增根目錄入口 |
| 共用核心規則 | 只有 CLAUDE.md，沒有共同來源 | 原有內容抽至 core-rules.md；兩入口明確要求完整讀取 |
| Skills | 專案沒有 SKILL.md，也沒有 .claude/skills 或 .agents/skills | 沒有既有專案 skills 可搬；已記錄未來單一來源配置 |
| 記憶 | docs/experience.md 已追蹤；另有專案外 Claude memory | 原文備份治理決策至 docs/memory/，加入索引與時效說明 |
| 技術債 | 見下節 | 記錄優先級，修正入口中的過時結構描述 |
| core-rules.md | 不存在 | 已建立；保留原物理規則、細則索引與驗證流程 |

## 技術債與後續工作

1. **高：重要工作尚未進 Git。** 檢查開始時 `CLAUDE.md`、整個 `.claude/`、
   `src/aqs/models/`、兩個 spinless 測試與 `docs/codebase_architecture.typ` 都未追蹤；
   `src/aqs/__init__.py` 有既有修改。新的 clone 不會包含這些工作。
   本次新增文件同樣尚未提交。應逐一審查並以明確檔名加入 Git；本次未 stage 或 commit。
2. **高：套件清單可能漏包。** `pyproject.toml` 的
   `[tool.setuptools] packages = ["aqs"]` 沒有列出 `aqs.models`，但目前
   `src/aqs/__init__.py` 已 import `.models.spinless`。正式 wheel 可能無法匯入。
   後續應修正套件發現設定並在乾淨環境建置／安裝驗證；此處是靜態發現，未做 wheel 測試。
3. **高：跨分支 gate convention 語義衝突。** 歷史文件記錄不同分支使用不同角度／符號。
   換工具或搬分支時不能只看 merge 是否衝突；涉及物理變更仍須依細則做獨立參考驗證。
4. **中：細則與規格仍位於 .claude/。** 這些都是可直接讀取的本地檔案，並非 Claude
   專用格式，但不能假設 Codex 自動載入 `paths` 規則。入口現已要求明確讀取全部七份。
   保留路徑是為了維持現有程式與規格引用；若未來移到中立目錄，需一次更新所有引用。
5. **中：歷史資料混入即時狀態。** 舊入口宣稱 spinless 模組不存在，已修正。
   `docs/experience.md` 和細則中的分支 SHA、測試數、工具版本、磁碟與硬體資訊仍是
   歷史紀錄。此次 `df -h .` 顯示專案磁碟使用率 81%、可用 169G，與「近 100%」不同。
   入口已提醒重新量測，不把歷史值當現況；本次未重新確認 GPU 或套件版本。
6. **中：執行環境不可只複製。** 保留 `pyproject.toml`、`uv.lock`、`.python-version`；
   新機器依 driver 選 extra，現有 GPU 環境維持 cu12。
   核心規則中的執行指令已加 `uv run --no-sync`，避免執行檢查時隱式同步依賴。
   `dist/` 有 2.0.0 wheel／tarball；不能假設與未提交的 src 一致。
7. **中：平台工具不等於專案 Skills。** 家目錄中存在 Claude／Codex 的 plugin cache，
   它們不會隨專案 clone，也不能由 cache 的存在推論已啟用。未逐一審查帳號設定、
   全域 hooks、權限或 MCP；若工作流程依赖它們，需要另外建立依賴清單。

## Skills 的共用方式

目前沒有專案 skills，因此本次不建立空的或虛構的 SKILL.md。新增時建議：

```text
.agents/skills/<name>/SKILL.md                 # 唯一正文
.claude/skills/<name> -> ../../.agents/skills/<name>
```

Codex 官方文件列出 `.agents/skills` 的專案探索位置，並支援 symlink skill folders。
Claude 端的相對連結配置仍需在實際安裝版本確認探索結果，尤其是 Windows checkout。
不要複製兩份正文；共享規則不包含平台專屬工具名稱、權限與 hooks 的假設。
本次沒有建立或實測任何 skill 連結。

## 記憶範圍與限制

已讀取此專案對應的 Claude `memory/MEMORY.md` 與其治理記憶。
主要決策大多已有本地規則對應，但原始治理記憶沒有在 repo 中完整備份，現已補上。
Claude 的此類記憶本身位於本機家目錄；「本機存在」不代表「專案可攜」。
沒有掃描全部歷史對話、其他專案的記憶或不可見的雲端記憶，因此不能保證不存在其他遺漏。

## 驗證與接手

- 靜態檢查兩入口引用同一份 core-rules.md、七份細則存在、記憶備份與來源逐字相同。
- 本次只改文件，未執行 pytest、GPU selftest 或重新建置套件。
- 新開 Codex 與 Claude Code 工作階段，請各自列出讀取的規則来源、兩模型規格權威與
  gate convention 注意事項，確認實際載入；此次未啟動另一個客戶端做端到端驗證。
- 完成檔案審查與提交後，再以乾淨 checkout 驗證安裝與對應測試。

官方依據：[Codex AGENTS.md 探索規則](https://learn.chatgpt.com/docs/agent-configuration/agents-md)、
[Codex Skills 位置與連結支援](https://learn.chatgpt.com/docs/build-skills)。

## 後續處理：2026-09-14 分類提交

已依需求將規格／規則、spinless 模組與守門測試、物理特徵測試、架構文件分別提交；
共用入口、記憶與本報告在遷移文件 commit 納入。原先未追蹤資料的風險已透過本地提交處理，
尚未 push。`pyproject.toml` 已補列 `aqs.models`，wheel 建置與解包後匯入驗證通過。
完整測試在沙箱外為 156 passed、1 skipped；其餘驗證限制詳見進度文件。
