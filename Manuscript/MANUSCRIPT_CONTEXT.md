# MR_Peak_Gene manuscript context

Last updated: **2026-10-15** (America/Chicago). Project root:
`/scratch/user/tfcui2025/CIMA/Claude_project`; paths below are relative to
`Manuscript/` unless prefixed by `../`.

This is a recovery record, not scientific approval. It reflects inspected
local files. Phase 3 did not re-audit analysis code, rerun analyses, access
websites, recompile LaTeX, or establish synchronization with Overleaf.

## A. Research overview

**Motivation/objective:** investigate whether genetic evidence can help link
accessible regulatory peaks to target genes, and compare MR-based links with
marginal correlations without presupposing universal superiority of MR.
Study cell types: `Bn_TCL1A`, `CD4_Tn_CCR7`, `CD8_Tn_CCR7`.

**Data:** CIMA caQTL/eQTL summary statistics (OMIX009786), candidate peak–gene
tables, donor-level pseudobulk RNA/ATAC, and external PBMC multiome WNN
metacells. MR uses an hg38 1000 Genomes EAS LD reference; S-LDSC uses separate
EUR resources. Primary data locations are documented in `../CLAUDE.md` §6;
the upstream provenance and sample-size description remain writing TODOs.

**Methods recorded in the current manuscript:** LD-aware GLS-IVW and GLS-MR
Egger, top-instrument SMR, equal-weight Cauchy combination of their p-values
(CCT), and WNN/pseudobulk Spearman correlations. Filtering precedes CCT and
BH adjustment within each gene and method; the main discovery threshold is
`q < 0.1`. CCT supplies combined significance, not a combined effect estimate.

**Validation/characterization:** regulatory-reference hypergeometric tests
and pooled enrichment; gene-wise ROC-AUC against EnhancerAtlas, FANTOM5,
4DGenome, GTEx, HiChIP, and PCHiC; Signac/JASPAR motif enrichment and
TF-expression checks; Fisher overlap with CIMA SCENIC+ networks; GO BP GSEA
and a separate gene-set resampling analysis; S-LDSC across nine GWAS;
positional GWAS overlap and locus illustrations; null-SNP z-score QC.

**Evidence map:** start with `README.md` (figure map and evidence for draft
numbers), `sections/results.tex`, and `sections/methods.tex`. Consult
`../Scripts/Pipeline.md` for current workflow descriptions,
`../Scripts/Visualization.md` for figure exports, and
`../Scripts_reference/00_README_Pipeline_Summary.md` for historical context.
Relevant implementation directories are `../Scripts/{Model,QC,Validation,Visualization}/`;
derived evidence is under `../Results/{MR,Validation,QC,Figures}/`.
Historical module numbering and completion checkmarks are not evidence of
current manuscript review. Verify selected code details when developing Methods.

Interpretation safeguards already present in the draft:

- Reference enrichment and ranking/AUC answer different questions; the
  existing ROC analysis selects direction automatically.
- CIMA GRN corroboration is within-study; positional GWAS overlap does not
  establish colocalization. The null-SNP diagnostic does not prove calibration.
- S-LDSC findings are exploratory: displayed stars use raw `p < 0.1`, there
  are no SE bars, and separate annotation fits do not directly test their
  difference. Conventional GO and empirical resampling nulls are distinct.
- Do not merge the two MR-only definitions: GO/S-LDSC use IVW/Egger/SMR
  `q < 0.1` and both marginal `q >= 0.1`; the GWAS screen includes CCT,
  uses MR `q < 0.001` and both marginal `q >= 0.05`
  (`sections/results.tex`, functional, heritability, and GWAS subsections).

## B. Target journal and requirements

**PLOS Genetics — Research Article.** Recorded official URLs (not fetched
in Phase 3):

- https://journals.plos.org/plosgenetics/s/submission-guidelines
- https://journals.plos.org/plosgenetics/s/latex

**Locally verified:** the supplied template identifies itself as version
3.8, April 2026. Its source requires `plos2025.bst`, a single submission
source without embedded graphics, an Abstract below 300 words, and an Author
Summary of 150–200 words. It directs funding, data availability, and author
contributions to submission-system fields. See
`template_reference/plos_latex_template.tex` and `submission_requirements.md`.

**Established Genetics/project targets:** full title ≤200 characters, short
title ≤70; Abstract ≤300 words; Author Summary 150–200 words; double spacing,
continuous line/page numbering, single column, at most three heading levels,
and numerical Vancouver references. The existing titles were measured as
100 and 48 characters. The generic template says 250 characters for the full
title; retain the agreed, stricter Genetics limit of 200.

**Agreed and implemented order:** Title Page → Abstract → Author Summary →
Introduction → Results → Discussion → Materials and Methods → Acknowledgments
→ References → Supporting Information Captions. This differs from the generic
template's illustrative order.

**STROBE-MR:** the prior journal review recorded a completed checklist with
manuscript page numbers and relevant text as a required supporting item.
`sections/supplementary.tex` reserves **S1 Checklist**, but no checklist exists.
This requirement and current Genetics-specific details were not independently
reconfirmed from the website in Phase 3. Reconfirm current journal instructions
and submission-system fields before submission; do not treat a TODO as compliance.

## C. Current local structure and build state

```text
Manuscript/
  AGENTS.md, MANUSCRIPT_CONTEXT.md       # instructions and living record
  main.tex                             # active entry point
  sections/
    title_page.tex, abstract.tex, author_summary.tex
    introduction.tex, results.tex, discussion.tex, methods.tex
    acknowledgments.tex, supplementary.tex
  references.bib, plos2025.bst          # database and official style
  figures/                             # eight PDFs, two XLSX workbooks
  figure_manifest.tsv                  # source paths, selected pages, hashes
  README.md, submission_requirements.md
  template_reference/                  # user-supplied ZIP and extracted template
  main.pdf, main.{aux,bbl,blg,log,out}   # published local draft/build outputs
  plos_draft.*                          # retained validation build
  overleaf_project.zip                  # source/asset upload bundle
  .plos-review/                        # prior backups and validation artifacts
```

The **last successful build** is 21 pages using pdfLaTeX/BibTeX. Inspection
of `main.log`, `main.aux`, and `.plos-review/validation_report.json` confirms
five resolved bibliography entries, nine figure labels, no LaTeX errors or
unresolved references/overflow, and one benign `epstopdf` shell-escape warning.
`main.pdf` matches the validated `plos_draft.pdf`; this was checked without
recompiling in Phase 3. The active official style matches the uploaded ZIP.

The current Overleaf ZIP has **25 files**, matching its local sources/assets
at inspection. It excludes build outputs and template-reference examples.
**Phase 3 does not regenerate the ZIP: these two new context files are not
inside it.** Transfer them separately or include them in a future authorized
package refresh. No newer Overleaf version has been checked.

For builds, follow `README.md` (including the server's indexed TeX search
workaround). Overleaf uses `main.tex`, pdfLaTeX, standard packages, relative
paths, and BibTeX; it needs no analysis environment or server-specific paths.

## D. Writing progress

Statuses describe the inspected content; no section is marked Reviewed.
Allowed vocabulary: Not started, Outline, Initial draft, Under revision, Reviewed.

| Section | Status | Existing content / important unfinished work |
|---|---|---|
| Title Page | Initial draft | Full and provisional short titles; authors, affiliations, correspondence, and applicable notes missing. |
| Abstract | Not started | Placeholder and guidance only; write after the Results narrative/figure choices are reviewed. |
| Author Summary | Not started | TODO only; needs a nontechnical 150–200-word summary. |
| Introduction | Outline | Five-point rationale/objective outline; needs prose and background/data-source citations. |
| Results | Initial draft | Seven subsections, numerical descriptions, two tables, five figures; biological interpretation, example selection, and scientific TODOs remain. |
| Discussion | Outline | Findings, complementary evidence, inference limits, and future directions outlined; full prose pending. |
| Materials and Methods | Outline | Thirteen subsections with short method summaries; mathematical details, instrument/input provenance, software/resources, and citations incomplete. |
| Acknowledgments | Not started | Placeholder; author-approved text, if applicable, missing. |
| References | Initial draft | Five entries: Bowden, Zhu, Liu/Xie, Benjamini/Hochberg, Finucane. Resource/software/GWAS references and metadata completion pending. |
| Supporting Information Captions | Initial draft | Four existing figure captions; provisional workbook descriptions and S1 Checklist placeholder. Final supporting-table captions/checklist missing. |

## E. Figures and tables

All items remain **draft**, including selection and placement; none is
recorded as finalized. All nine figure displays have existing captions.
Filenames retain historical numbers, which differ from manuscript numbering.
Source/page/hash details are in `figure_manifest.tsv`; preserve original assets.

| Item | Current file/page | Related section | Included in LaTeX? / pending work |
|---|---|---|---|
| Fig 1 | `figures/Fig1_MR_overview.pdf`, pp. 1, 4 | Results: shared/discordant links | Yes, two unchanged pages in one display; review presentation. |
| Fig 2 | `figures/Fig2_examples.pdf`, p. 1 | Results: examples | Yes; dense panels, example choice and caption need review. |
| Fig 3 | `figures/Fig4_AUC.pdf`, p. 1 (Bn) | Results: reference support | Yes; CD4/CD8 pages retained in source; agree final cell-type coverage. |
| Fig 4 | `figures/Fig5_SLDSC.pdf`, p. 2 (Combined) | Results: heritability | Yes; decide uncertainty presentation and exploratory wording. |
| Fig 5 | `figures/Fig8_Link_GWAS.pdf`, p. 7 | Results: GWAS/TMEM258 | Yes; confirm locus choice, interpretation, and literature. |
| S1 Fig | `figures/Fig7_GRN.pdf`, p. 1 | Results: GRN; supporting section | Yes; decide high-quality-network detail and placement. |
| S2 Fig | `figures/Fig9_GO_GSEA.pdf`, p. 2 | Results: GO; supporting section | Yes; conventional BH comparison displayed; resolve intended GO narrative. |
| S3 Fig | `figures/Fig5_SLDSC.pdf`, p. 1 | Results: heritability; supporting section | Yes; retain distinction from Combined and review uncertainty. |
| S4 Fig | `figures/Sup1_OneSampleMR.pdf`, p. 1 (Bn) | Methods/QC; supporting section | Yes; CD4/CD8 retained in source; descriptive diagnostic only. |
| Table 1 | `sections/results.tex`, `tab:analysis-set` | Results: common analysis set | Yes, native LaTeX; review counts and wording with Results. |
| Table 2 | `sections/results.tex`, `tab:significant-links` | Results: discovery yield | Yes, native LaTeX; retain within-gene BH interpretation. |
| S1 Table | `figures/Fig3_Hyper_table.xlsx` | Results: reference enrichment | Mentioned, not rendered; finalize caption and table selection/format. |
| S2 Table | `figures/Fig6_Motif_table.xlsx` | Results: motifs | Mentioned, not rendered; inherited annotations and misleading “MR-only” title need review. |
| S1 Checklist | No file yet; reservation in `sections/supplementary.tex` | STROBE-MR supporting item | Placeholder only; complete and map to final manuscript pages. |

## F. Established writing decisions

- Target PLOS Genetics Research Article; write in Academic English.
- Maintain `main.tex` plus modular section files for development; do not
  recreate the manuscript or flatten it during routine drafting.
- Use the official, unmodified `plos2025.bst` with `references.bib` and
  numerical citation-order references. The `\citep` compatibility alias
  preserves existing section sources. Some entries have only three named
  authors plus `others`; full metadata and journal abbreviations need review.
- Keep figures visible (`\draftfigurestrue`); use selected PDF pages without
  extraction/conversion or asset changes. Publication TIFF/EPS preparation is
  deferred. Final submission needs separate graphics and a single source.
- Keep TODOs visible for drafting (`\showtodostrue`). Turning them off does
  not remove literal placeholders or inline caption TODOs.
- Preserve relative paths and ordinary pdfLaTeX/BibTeX compatibility.
  Declarations are tracked in `submission_requirements.md`; do not fabricate
  them or place funding information in Acknowledgments.

## G. Cross-platform collaboration

| Environment | Intended role |
|---|---|
| Codex in VS Code | LaTeX editing, figure integration, compilation/formatting, project maintenance. |
| ChatGPT Work | Literature synthesis, manuscript drafting/revision, work with reference papers and project files. |
| ChatGPT Chat Mode | Scientific discussion, statistical reasoning, paragraph-level writing/editing. |
| Overleaf | Primary manuscript editing, PDF compilation, final review. |

These environments **do not automatically share all files or histories**.
Before substantial edits, compare the local files with the latest user-provided
Overleaf export/revisions and agree which version to edit. Never import this
ZIP over newer Overleaf work blindly. Provide these context files and current
relevant sections when handing off work; record major agreed changes here.

## H. Priority open tasks

1. Review the existing Results narrative, central message, and main/supporting
   figure and locus selections with the researcher before expanding prose.
2. Complete biological interpretation of motifs, GRN, GO, and disease examples;
   resolve GO comparisons and S-LDSC uncertainty reporting.
3. Develop Methods with targeted source-code checks: candidate/instrument
   provenance, upstream covariates/sample sizes, GLS/SMR details, resources,
   coordinate rules, backgrounds, and uncertainty conventions.
4. Complete missing study/resource/software/GWAS citations and bibliography
   metadata; do not infer reference identities from filenames.
5. Expand Introduction and Discussion; write Abstract and Author Summary
   after the Results/figure message is agreed.
6. Supply authorship/declarations, final supporting-table captions and
   STROBE-MR checklist; review final figure legends/readability and journal
   requirements. Use `submission_requirements.md` for the detailed checklist.
7. At a later authorized submission milestone, prepare single-source,
   caption-only submission files and separate publication figures, and refresh
   the Overleaf package including current context files as appropriate.

## I. Latest progress and resume point

- **2026-10-15, Phase 3:** created these two persistent-context files after
  inspecting actual sources, documentation, existing build evidence, assets,
  and ZIP contents. No manuscript content or other files changed.
- **Last manuscript milestone:** official local PLOS template adaptation;
  successful 21-page draft and matching Overleaf source package. Formatting
  is established; scientific writing remains an early draft, not submission-ready.
- **Recommended next writing task:** review/refine Results and figure selection
  with the researcher, then prioritize complete Materials and Methods.
- **Resume:** read `AGENTS.md` and this file, verify the latest local/Overleaf
  version, inspect only the sections needed for the requested task, and proceed
  within that scope. Update this record after major milestones or structural
  changes, not individual wording edits. Wait for the user's next task.
