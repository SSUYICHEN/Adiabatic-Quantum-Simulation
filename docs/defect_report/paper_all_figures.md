# 論文全部圖形：修正前 | 修正後

日期：2026-09-17。以論文參數（N=6；Fig. 2 掃 T∈{1,5,15,80}、L∈{1,…,150}；Fig. 3、4 用 T=1、L=40）
在 CUDA-Q `nvidia-fp64` 上重算每一張圖的兩個版本：

- **修正前**：`main`（0fe9d1d）的態製備，即論文所用程式。重算結果與論文原圖逐點一致
  （Fig. 2 上圖 L=1、T=1 的 0.14 與 T=80 的 0.75 平台；Fig. 3 下圖 −0.49；Fig. 4 下圖 0.72），
  確認整條產圖流程都被重現。
- **修正後**：本分支，三個缺陷皆已修正。

| 檔案 | 內容 |
|---|---|
| `figures/paper_all_compare.svg` | 六列兩欄總覽：Fig. 2 上／下、Fig. 3 上／下、Fig. 4 上／下 |
| `figures/paper_fig2_compare.svg`、`paper_fig3_compare.svg`、`paper_fig4_compare.svg` | 每張論文圖各自的 2×2 對照 |
| `paper_fig2.json`、`paper_impact.json` | 全部數據 |
| `paper_fig2_recompute.py`、`paper_impact_recompute.py`、`paper_all_figures_plot.py` | 重算與繪圖腳本 |

## 逐圖結論

| 圖 | 修正後的變化 | 幅度 |
|---|---|---|
| Fig. 2 上（U=0）與下（U=1） | **完全不變**。fidelity_scan 以製備態自身為參考；最高能帶行列式態同樣是 H₀ 的本徵態 | 全部 88 點最大差 3e-13 |
| Fig. 3 上，ΔU=0 階梯 | 不變 | — |
| Fig. 3 上，ΔU≠0 曲線 | 拓樸相（w>1）的偏離由「往下」變「往上」；平庸相（w<1）由略低於 1 變略高於 1 | ΔU=0.3 時 ∓0.045 |
| Fig. 3 下（w=2） | 曲線翻到 γ>0，且大 ΔU 時幅度較小 | ΔU=3：−0.49 → +0.35 |
| Fig. 4 上，ΔU=0 | 不變 | — |
| Fig. 4 上，ΔU≠0 | 體內偏移由負變正；cell 1 邊緣極化改為增強、cell 6 改為減弱 | ΔU=2 體內 −0.53 → +0.47 |
| Fig. 4 下（cell 1） | 由下降改為上升 | ΔU=2：0.72 → 1.23 |

## 三個缺陷各自對這些圖的關係

- **缺陷 A（Givens 反號）**：唯一改變論文圖的缺陷。等價於 U_A、U_B 全部反號（見 `paper_impact.md`）。
- **缺陷 B（CUDA-Q 位元序）**：論文以 qiskit 後端產圖，不經過這條路徑。若當時用了舊的 GPU 後端，
  態向量會位元反轉：極化剖面變成 P′_j = −P_{N−1−j}、Berry 相位變成 γ′ = −γ，論文圖會再翻一次；
  重算的「修正前」與論文原圖逐點一致，證明論文沒有走這條路。
- **缺陷 C（OpenFermion 端序）**：只影響 `aqs measure` 讀 JSON 哈密頓量的路徑，論文三張圖都不使用。

重跑：

```bash
uv run --no-sync python docs/defect_report/paper_fig2_recompute.py     # ~18 min GPU
uv run --no-sync python docs/defect_report/paper_impact_recompute.py   # ~18 min GPU
uv run --no-sync python docs/defect_report/paper_all_figures_plot.py
```
