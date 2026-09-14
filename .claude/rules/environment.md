---
paths:
  - "pyproject.toml"
  - "uv.lock"
  - "requirements*.txt"
  - ".python-version"
---

# Environment & Workflow Constraints

## uv / GPU packages

- On this GPU machine (RTX 4090, driver = CUDA 12.6), **always** `uv sync --extra cu12`.
  Plain `uv sync` uninstalls ~4 GB of GPU wheels; disk is nearly full (`/` and `/mnt/sata`
  both near 100%), so re-fetching may fail — check `df -h` before anything wheel-sized.
- Python venv is 3.10; the `cu13` extra requires ≥3.11 and pins CUDA 13 (breaks on this
  driver). `cu12`/`cu13` are mutually exclusive by design (`[tool.uv] conflicts`) — do not
  "upgrade" the extra or the pin.
- `pytest`/`coverage` live in the `dev` dependency group (installed by `uv sync` by
  default); `requirements-dev.txt` is only a pip fallback.
- `[tool.uv] default-extras` is not a valid field in uv 0.9.x — don't add it.

## Git hygiene (each rule below cost real damage once)

- `git add` explicit filenames only — never `git add <directory>` (it swept another
  branch's untracked files into commits, twice).
- Never `git stash -u` (it once stashed the 4.2 GB `.venv`); `git stash list` before any
  `pop`. For cross-revision comparison use a worktree instead:
  ```bash
  git worktree add -q --detach /path/wt <rev>
  PYTHONPATH=/path/wt/src .venv/bin/python script.py out.json
  git worktree remove --force /path/wt
  ```
- No `gh` CLI on this machine — PRs are opened manually.
- ⚠️ Branch semantics: `fix/core-silent-convention-bugs` and
  `fix/cudaq-bit-order-and-uv-migration` use **incompatible gate-angle conventions** in
  different files; merging them is a silent semantic conflict git will not flag (measured
  damage: adiabatic convergence 0.99988 → 0.930398). Merge order and required call-site
  edits are documented in `docs/experience.md` §5.

## Misc tooling

- typst 0.15 (snap) cannot read `/tmp` — put its temp files inside the project; no
  `angle.l`/`angle.r`/`bra`/`ket` in this typst version; Noto CJK TC available, New
  Computer Modern **Sans** is not.
