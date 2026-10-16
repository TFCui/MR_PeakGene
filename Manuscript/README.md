# PLOS Genetics manuscript working draft

The existing Phase 2 manuscript has been adapted to the user-supplied official
PLOS LaTeX template (version 3.8, April 2026), preserved under
`template_reference/`. The active bibliography style is the unmodified official
`plos2025.bst`. No network access or template download was used for this update.
See `submission_requirements.md` for remaining manuscript, supporting-file,
and submission-system requirements and the final single-source/figure-free
submission step. Modular section sources and visible figures are retained for
development. Scientific text, existing caption text, bibliography entries,
and all research figure/table assets are preserved; only the Methods and
supporting-section headings were renamed.

The section order is Title Page, Abstract, Author Summary, Introduction,
Results, Discussion, Materials and Methods, Acknowledgments, References, and
Supporting Information Captions. The abstract and new Author Summary remain
explicit placeholders. The PLOS generic template's 250-character title comment
is superseded by the requested Genetics limit of 200 (short title: 70).

This is a concise, editable working draft, prepared on 15 October 2026. It is
not a publication-ready manuscript. The abstract is a placeholder; the
Introduction, Methods, and Discussion are outlines or short summaries. Results
contain supported descriptive findings and visible TODOs. No authors,
affiliations, new analyses, or regenerated analysis figures were added.

## Files and editing

- `main.tex`: official-template-based `article` entry point and shared formatting;
  double spacing, continuous line numbering, single column, page numbering.
- `plos2025.bst`: official numerical, citation-order Vancouver bibliography style.
- `sections/title_page.tex`, `author_summary.tex`, `acknowledgments.tex`:
  new front-matter and acknowledgment placeholders; no authors invented.
- `sections/abstract.tex`, `introduction.tex`, `results.tex`, `methods.tex`,
  `discussion.tex`: modular manuscript text.
- `sections/supplementary.tex`: provisional supplementary displays and native
  supporting-table descriptions.
- `references.bib`: five verified foundational method references; the remaining
  source/software/GWAS bibliography is explicitly TODO.
- `figures/`: byte-identical copies of eight existing multipage/single-page
  PDFs and two supporting XLSX workbooks. No raster review previews are used
  in the manuscript.
- `figure_manifest.tsv`: original paths, page counts, selected pages, placement,
  and SHA-256 hashes of copied assets.
- `main.pdf`: compiled review copy; see build status below.
- `overleaf_project.zip`: source-only upload bundle; excludes build products
  and temporary review files.

Edit section files directly. Set `\showtodosfalse` in `main.tex` to hide
editorial TODO paragraphs; caption TODOs remain explicit text until edited.
Keep `\draftfigurestrue` while drafting. The `\draftfiguresfalse` switch hides
images while retaining captions and labels for later submission preparation.
Figure/table numbering is assigned by LaTeX and is independent of the
historical figure numbers embedded in filenames. `q` means the existing
within-gene, method-specific BH-adjusted link p-value, not a genome-wide FDR.

## Overleaf and local compilation

Upload `overleaf_project.zip` to a new Overleaf project, set
`main.tex` as the main document, and choose **pdfLaTeX**. The project uses only
standard packages and relative paths; no server data, R, Python, shell escape,
or analysis pipeline is required to compile. XLSX files are native supporting
tables, not LaTeX graphics.

For a local build, from `Manuscript/`:

```bash
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

On this server, the already-installed compiler can be exposed with
`module load GCC/13.2.0 texlive/20230313`. No software installation is needed.
Compilation writes should remain inside `Manuscript/`.

The server module sets `TEXMFHOME` to the full shared TeX distribution, which
caused repeated recursive searches and a four-minute timeout before processing
the document body. Local validation used an empty manuscript-local `TEXMFHOME`
and the installed filename database instead. To avoid that delay here, after
loading the module and before the commands above, use:

```bash
mkdir -p .plos-review/texmf-home .plos-review/texmf-var .plos-review/texmf-config
export TEXMFHOME="$PWD/.plos-review/texmf-home"
export TEXMFVAR="$PWD/.plos-review/texmf-var"
export TEXMFCONFIG="$PWD/.plos-review/texmf-config"
export TEXINPUTS=".:!!/sw/eb/sw/texlive/20230313-GCC-13.2.0/texmf-dist/tex//"
```

These server-specific search overrides are unnecessary on Overleaf. All local
LaTeX passes were bounded by a six-minute timeout (BibTeX: one minute).

The PLOS adaptation compiles successfully with pdfLaTeX and BibTeX to a
21-page PDF. All five bibliography entries and nine figure labels resolve;
there are no LaTeX errors, undefined references/citations, or text/float
overflow warnings. The retained template `epstopdf` package issues an expected
warning because local validation disables shell escape; the draft uses PDF
graphics and requires no conversion. Source and asset hashes matched the
saved pre-formatting baseline; scientific prose/captions, all bibliography
entries, and all ten PDF/XLSX assets were checked for preservation.

The Overleaf ZIP includes the active LaTeX sources, official bibliography
style, bibliography database, figure assets, manifest, and drafting notes.
It excludes compiler products, review caches/backups, and the local template
reference archive/examples. Set `main.tex` as the main document.

## Proposed scientific narrative

Preliminary title: **Mendelian randomization and multi-omic evidence for
regulatory peak–gene links in human immune cells**.

1. Establish the common tested pair sets and six-method comparison.
2. Describe shared and method-specific discoveries and effect directions.
3. Evaluate reference support using enrichment and gene-wise ROC-AUC, keeping
   the metrics distinct.
4. Use CIMA regulatory networks as within-cohort corroboration.
5. Explore motif and gene-set relevance without interpreting motif counts as
   proof of specificity or treating the two GO null comparisons as identical.
6. Explore disease heritability of MR-only regions.
7. Present GWAS-overlapping regulatory examples as candidates for follow-up.

The observed tradeoff is more marginal-multiome discoveries versus stronger
pooled regulatory-reference enrichment for MR-selected sets; AUC differences
are smaller and benchmark-dependent. This is descriptive and does not imply a
universal superiority claim. Detailed mathematical development and a complete
literature review are deferred.

## Figure-to-Results map and status

Original paths below are relative to the project root (one directory above
this README). `new/` abbreviates `Results/Figures/Figures_20261019/`, and `old/`
abbreviates `Results/Figures/Figures_20261012/`.

All placements, captions, and final panel choices are **draft**. A documented
completion checkmark establishes that an export exists, not publication
approval. Representative pages of all selected PDF assets were visually
reviewed, all source PDF page counts were checked, and underlying selected
values were compared with result tables. Not every page of each source PDF was
visually reviewed.

| Source artifact | Underlying analysis and scientific purpose | Results subsection / current placement | Status |
|---|---|---|---|
| `new/Fig1_MR_overview.pdf` (4 pp.) | Six-method discovery overlap and effect-direction concordance | Common analysis set; shared/discordant links. Main Fig. 1 uses pages 1 and 4. Pages 2–3 retained in the copy. | Current export; manuscript arrangement provisional |
| `new/Fig2_examples.pdf` (6 pp.) | QTL/MR/marginal measurements for selected pairs | Shared/discordant links. Main Fig. 2 uses page 1; other example pages retained. | Current export; dense panels need later editorial review |
| `new/Fig3_Hyper_global.tsv` (108 rows) and `new/Fig3_Hyper_table.xlsx` (3 sheets) | Pooled reference overlap, global hypergeometric p-values, and fold enrichment | Regulatory reference support; workbook copied as supporting Table S1. TSV remains in its original location. | Current table output; no current Fig3 PDF exists |
| `old/Fig4_AUC.pdf` (3 pp.) | Gene-wise ROC-AUC against six references | Regulatory reference support. Main Fig. 3 uses Bn page 1; CD4/CD8 retained as pages 2–3. | Latest available AUC export, despite older directory; automatic ROC direction noted |
| `new/Fig5_SLDSC.pdf` (2 pp.) | MROnly versus all-candidate peak heritability enrichment | Disease heritability. Main Fig. 4 uses Combined page 2; Supplementary Fig. S3 uses cell-type page 1. | Current exploratory export; raw p<0.1 stars, no SE bars |
| `new/Motif_all.tsv` (13,428 rows) and `new/Fig6_Motif_table.xlsx` (1 sheet) | Method-specific motif enrichment and selected motif annotations | Functional relevance; workbook copied as supporting Table S2. TSV remains in its original location. | Current numeric cells; inherited manual annotations and misleading “MR-only” title require review |
| `old/Fig7_GRN.pdf` (1 p.) | Counts and Fisher odds ratios against CIMA SCENIC+ networks | CIMA network corroboration; Supplementary Fig. S1 | Latest available GRN export; within-CIMA evidence |
| `new/Fig8_Link_GWAS.pdf` (27 pp.) | Disease-associated SNPs inside peaks of strong MR-supported links | GWAS examples. Main Fig. 5 uses page 7 (CD4 TMEM258/asthma); other loci retained. | Current export; illustrative locus selection provisional |
| `new/Fig9_GO_GSEA.pdf` (3 pp.) | GO BP GSEA and a separate empirical gene-set comparison | Functional relevance; Supplementary Fig. S2 uses conventional-BH page 2. Pages 1 and 3 retained. | Current exploratory output; conventional and empirical comparisons distinguished |
| `new/Sup1_OneSampleMR.pdf` (3 pp.) | Null-SNP z-score correlation diagnostic | Methods and supplementary QC; Supplementary Fig. S4 uses Bn page 1; CD4/CD8 retained. | Current diagnostic export; not a calibration proof |

The copied complete PDFs retain all their original pages. The manuscript uses
the `page=` option of `\includegraphics` to display selected pages without
extracting, cropping, changing, or regenerating the underlying assets.

### Other existing artifacts reviewed as historical alternatives

| Artifact in `old/` | Contents and analysis | Decision / status |
|---|---|---|
| `Fig1_MR_overview.pdf` (5 pp.) | Earlier discovery and concordance presentation | Superseded by current 4-page overview; not copied |
| `Fig2_examples_Bn_TCL1A.pdf` and `_v2.pdf` (22 pp. each) | Earlier Bn example layouts | Superseded by the combined current example PDF; not copied |
| `Fig2_examples_CD4_Tn_CCR7.pdf` and `_v2.pdf` (23 pp. each) | Earlier CD4 example layouts | Superseded; not copied |
| `Fig2_examples_CD8_Tn_CCR7.pdf` and `_v2.pdf` (24 pp. each) | Earlier CD8 example layouts | Superseded; not copied |
| `Fig3_Hyper.pdf` (3 pp.) | Earlier gene-level hypergeometric p-value histograms | Replaced by the current table presentation; not copied |
| `Fig5_SLDSC.pdf` (1 p.) | Earlier six-method, FDR 0.05 heritability presentation | Superseded by current MROnly/pool analysis; not copied |
| `Fig6_Motif_table.xlsx` and `Motif_all.tsv` | Older motif presentation and selected annotations | Current directory chosen; older versions not copied |
| `Fig8_Link_GWAS.pdf` (32 pp.) | Earlier GWAS-locus presentation | Superseded by current 27-page PDF; not copied |

## Evidence supporting draft numbers

These were lightweight summaries of existing derived results, not reruns of
statistical analysis pipelines. Source code was consulted to interpret the
figures; no additional statistical audit was undertaken in this phase.

| Draft statement | Existing evidence, relative to project root |
|---|---|
| Tested pairs, genes, peaks, donor/metacell counts; six-method significance counts | `Results/MR/<CT>_AllPairs_results.tsv.gz`; interpretation in `Scripts/Model/MR_AllPairs.R:292–411` |
| All-six intersections and sign concordance | The same MR tables; sign definitions in `Scripts/Visualization/Plot1_MRresults.R:68–79` |
| Example gene/peak identities | `Scripts/Visualization/Example_pairs.tsv` (`Fig2_selected`, page 1); existing PDF and MR tables |
| Pooled folds >1; SMR largest in 17/18 cell-type/reference combinations | `Results/Figures/Figures_20261019/Fig3_Hyper_global.tsv`; `Plot3_Hyper.R:20–30` |
| Mean ROC-AUC examples | `Results/Validation/Valid4_AUC/<CT>_AUC_long.tsv.gz`; `Valid4_AUC.R:83–96` |
| GRN membership, overlap counts, and odds ratios | `Results/Validation/Valid5_GRN/GRN_summary.tsv` and `GRN_sig_counts_by_method.tsv` |
| Motif counts | `Results/Validation/Valid3_Motif/<CT>_motif_enrichment.tsv.gz`, using p.adjust<0.05 and fold.enrichment>1 |
| GO gene/term counts and separate significance measures | `Results/Validation/Valid7_GO/GO_GSEA_MROnly.tsv.gz`, using `p.adjust` versus `emp_fdr`; `Valid7_GO_1_gsea.R`, `_2_perm.R`, `_3_merge.R` |
| S-LDSC peak counts and enrichment/p-values | `Results/Validation/Valid2_SLDSC/set_stats.tsv` and `SLDSC_all.tsv`; prep and merging scripts |
| 27 GWAS-overlapping pairs, 21 genes, 14 peaks, 37 rows; TMEM258 example | `Results/Validation/Valid6_GWAS/GWAS_MRonly_overlap.tsv` and `GWAS_pair_counts.tsv`; `Valid6_GWAS.R` |
| Null-SNP diagnostic definition | `Scripts/QC/QC_OneSampleMR.R`; `Results/QC/<CT>_global_rho.tsv`; existing supplementary PDF |

## References and writing TODOs

Verified foundational sources are linked in each BibTeX entry: MR-Egger
([Bowden et al.](https://academic.oup.com/ije/article/44/2/512/754653)), SMR
([Zhu et al.](https://www.nature.com/articles/ng.3538)), CCT
([Liu and Xie](https://doi.org/10.1080/01621459.2018.1554485)), BH
([Benjamini and Hochberg](https://onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1995.tb02031.x)),
and S-LDSC ([Finucane et al.](https://www.nature.com/articles/ng.3404)). These
foundational references do not imply that all implementation details are
identical to the original papers. Multi-author entries currently use an
explicit `and others` abbreviation.

Priority TODOs for the next writing round:

1. Review title, Results order, selected examples, and main/supplement placement.
2. Complete source citations for CIMA, external multiome, regulatory resources,
   software, motifs, GO, GRN, and each GWAS; do not infer them from filenames.
3. Develop TF and disease-locus interpretation with supporting literature.
4. Complete instrument/input provenance and concise mathematical Methods.
5. Decide the intended GO comparison and wording of exploratory enrichment.
6. Replace the abstract placeholder and develop the introductory/discussion
   outlines after the Results message is agreed.
7. Later revise dense example panels, blank space in inherited PDFs, supporting
   table formats, and final figure legends. No such figure revisions were made
   in Phase 2.

Original figures, analysis results, scripts, documentation, and Git state must
remain unchanged during manuscript-only editing. Compilation products and
temporary visual checks are confined to this directory.
