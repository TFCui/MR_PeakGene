# Manuscript rules (Claude Code; CodeX for the final review only)

Project MR_PeakGene: MR-based peak–gene links from CIMA single-cell QTLs. Target: PLOS Genetics,
Research Article. Academic English, LaTeX; the author writes on Overleaf. The author decides the
structure; Claude's core job is polishing the English (fluent, idiomatic, precise, academic) while
keeping track of what each paragraph should and should not say given the whole paper.

## Read first
1. `notes/1_WORKFLOW.md`: roles, one-round workflow (W1 review → author answers → W2 polish →
   W3 close), polishing rules (incl. scMultiMap as a style model), red lines, prompt templates.
2. `notes/2_BOARD.md`: sections 1–2 (unified wording, whole-paper map), then only the current unit
   in section 4.
3. `notes/3_LEDGER.md` when needed: verified numbers (section 1), figure ledger, submission,
   Overleaf sync, update log.
Do not create other notes files. Continue from the actual files, not from chat history.

## Red lines
- Never change numbers, P values, thresholds, sample sizes, `\cite` keys, `\label`s, equations or
  figure/table content; report suspected errors instead.
- Never invent results, method details, citations or conclusions; use `\TODO{}` or `[CITE: claim]`.
- Interpretive sentences come only from the author's points.
- Write only inside `Manuscript/` (helper scripts: `../tmp/Manuscript_tables/`; round snapshots:
  `../tmp/ms_rounds/`). `../Scripts/`, `../Results/`, `../paper_reference/` and all data are
  read-only; no analyses or SLURM jobs unless asked. If a step needs a write elsewhere, stop and ask.
- No commit or push (GitHub or Overleaf) unless the author asks in this session.

## Build and sync
- Build both documents: `bash -l ../tmp/build_ms.sh` (`main.tex` and `supplementary.tex`).
  Supporting items live in `supplementary.tex` and are cited with `\nameref{S1_Fig}`.
- Checks: `python tools/ms_check.py invariants [--snapshot DIR | --base REF]` and `placeholders`.
- Overleaf: before editing, run `bash tmp/overleaf_sync.sh status` from the project root and
  `pull` if Overleaf is newer. `push` only when asked. Never print the token. Details:
  `notes/3_LEDGER.md` section 4.

## End of every task
W2: report the change log, expression notes, checklist and diff command. W3: close the work card,
update `2_BOARD` 2.2 / 2.3, and append one line to `notes/3_LEDGER.md` section 5
(date | tool and model | what | how checked). Figure, format and sync changes are logged there too.
