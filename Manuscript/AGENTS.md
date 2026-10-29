# Manuscript collaboration instructions

Rules for the AI assistants working in `Manuscript/` (Claude Code and CodeX).
This file holds shared rules only and should stay stable. Paper facts and
progress live in each tool's own notes file (see Roles).

## Roles

| Tool | Responsible for | Own notes file |
|---|---|---|
| Claude Code | Organizing LaTeX, formatting, checking cross-references, compiling the paper; figure and table provenance | `Notes_Claude_Code.md` |
| CodeX | Reading Methods, checking equations, discussing scientific logic, drafting paragraph structure | `CODEX_METHODS_NOTES.md` (organized by paper section) |

- Each tool reads both files but writes only its own. Do not edit the other
  tool's notes file.
- Editing `sections/*.tex`: CodeX edits Methods/Results/Discussion wording and
  equations when the user asks; Claude Code does everything else in LaTeX
  (layout, floats, labels, `\ref`/`\cite` checks, tables, figure inclusion,
  build). Never revert the other tool's text changes; if something looks
  wrong, say so to the user.
- `Notes_Claude_Code.md` section D describes what the `.tex` files currently
  contain; Claude Code refreshes it from the actual files, including text
  that CodeX changed. CodeX records its scientific reasoning, formula checks,
  and paragraph plans in `CODEX_METHODS_NOTES.md`.

## Project identity

- Project: **MR_Peak_Gene**.
- Field: Statistical Genetics and Computational Genomics.
- Target: **PLOS Genetics**, Research Article.
- Language: Academic English. Format: LaTeX.
- Primary manuscript editing platform: Overleaf.

## Start and finish of every task

1. Read your own notes file, then skim the other tool's notes file.
   `Notes_Claude_Code.md` records, for every figure and table, where it comes
   from, how it was changed before entering the paper, and what it means
   (section C), and what each manuscript section currently contains and still
   needs (section D).
2. Inspect the relevant LaTeX sections and `main.tex`.
3. Verify which version is current, including any newer Overleaf revisions;
   do not overwrite newer work with an older local copy or ZIP.
4. Continue from the actual files, not solely from conversation history.
5. When finished, update your own notes file:
   - Claude Code: a figure or table added, replaced, or changed -> its entry in
     section C of `Notes_Claude_Code.md` (source, modification, meaning, hash;
     refresh hashes with `sha256sum figures/*.pdf`); a section changed ->
     section D; always one line in section E (date | tool | what changed).
   - CodeX: add the discussion, checks, and decisions under the matching
     section heading in `CODEX_METHODS_NOTES.md`, with the date and what was
     changed in LaTeX, if anything.
6. Summarize important modifications and relevant validation to the user.

Platforms do not share session history. Hand work across tools through the two
notes files and the manuscript sources, not through chat history.

## Writing principles

- Follow PLOS Genetics author guidelines and preserve the official PLOS
  formatting conventions already adopted here.
- Preserve scientific accuracy. Never invent results, methods, references,
  authorship details, declarations, or scientific conclusions.
- Use actual project results as evidence. Documentation is background;
  verify important methodological details against the relevant source code
  when writing them. Distinguish verified, documented, exploratory, and
  missing information. Use TODO placeholders where information is missing.
- Preserve existing scientific content unless its revision is requested.
  Preserve LaTeX labels, citation keys, and cross-references.
- Keep sources modular and editable. Retain visible figures during drafting;
  prepare a separate submission version when requested. Preserve research
  figure assets and the official bibliography style.

## File safety

- Edit files only inside `Manuscript/`; work incrementally and touch only
  files relevant to the current task. Do not rewrite unrelated sections.
- Everything outside `Manuscript/` is read-only: `../Scripts/`, `../Results/`,
  `../paper_reference/`, `../Scripts_reference/`, and all data. Do not modify
  research scripts or analysis results as part of manuscript work, and do not
  run analyses or submit SLURM jobs without an explicit task. If a step seems
  to require writing outside `Manuscript/`, stop and ask the user.
- Temporary files and generator scripts for manuscript tables belong in
  `../tmp/` (for example `../tmp/Manuscript_tables/`); compilation products
  stay inside `Manuscript/`.
- Build: `bash -l ../tmp/build_ms.sh` (loads the TeX module and builds both
  root documents: `main.tex` with pdfLaTeX/BibTeX, and `supplementary.tex`
  with pdfLaTeX twice; the script documents its own settings). Supporting
  figures and tables live in `supplementary.tex`; `main.tex` only lists them
  in `sections/supporting_info.tex` and cites them with `\nameref{S1_Fig}`.
- Do not commit or push without explicit permission.

## Reference papers

`../paper_reference/` holds the core papers behind the Methods: MR-Egger
(Burgess and Thompson 2017 and appendix), GRIP, SMR (Zhu et al. 2016), and the
CIMA Science paper. Read-only. `Notes_Claude_Code.md` section A says where
each is used and which `references.bib` key it maps to. Prefer them over
memory when writing or checking Methods, Introduction, and Discussion.

## Network access

- Prefer local documentation and the existing `template_reference/` files.
  Avoid unnecessary network requests from the HPC server.
- Respect task-specific network restrictions. If access is necessary and
  permitted, use explicit, bounded timeouts; do not repeatedly retry failures.
  Report unavailable resources instead of leaving requests hanging.
