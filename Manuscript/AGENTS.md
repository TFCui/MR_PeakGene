# Manuscript collaboration instructions

## Project identity

- Project: **MR_Peak_Gene**.
- Field: Statistical Genetics and Computational Genomics.
- Target: **PLOS Genetics**, Research Article.
- Language: Academic English. Format: LaTeX.
- Primary manuscript editing platform: Overleaf.

## Recover context before substantial work

1. Read `MANUSCRIPT_CONTEXT.md` in this directory.
2. Inspect the relevant current LaTeX sections and `main.tex`.
3. Verify which version is current, including any newer Overleaf revisions;
   do not overwrite newer work with an older local copy or ZIP.
4. Consult `README.md`, `submission_requirements.md`, and project documentation
   only as needed. Research documentation is indexed in the context file.
5. Continue from the actual files, not solely from conversation history.

Other platforms do not automatically share these files or session history.
When handing work across platforms, provide the context file and the latest
relevant manuscript sources explicitly.

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

## Collaboration and maintenance

- Work incrementally and edit only files relevant to the current task.
  Do not rewrite unrelated sections.
- Do not modify research scripts or analysis results as part of manuscript
  work. Do not run analyses or submit SLURM jobs without an explicit task.
- Keep compilation products and temporary manuscript files within
  `Manuscript/`; follow the existing server build notes in `README.md`.
- Summarize important modifications and relevant validation after each task.
- Update `MANUSCRIPT_CONTEXT.md` after major writing milestones or structural
  revisions, not minor wording changes. Keep this `AGENTS.md` stable.
- Do not commit or push without explicit permission.

## Network access

- Prefer local documentation and the existing `template_reference/` files.
  Avoid unnecessary network requests from the HPC server.
- Respect task-specific network restrictions. If access is necessary and
  permitted, use explicit, bounded timeouts; do not repeatedly retry failures.
  Report unavailable resources instead of leaving requests hanging.
