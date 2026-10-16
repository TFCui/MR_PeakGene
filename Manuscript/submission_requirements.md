# PLOS Genetics submission preparation

This is a formatting adaptation of an early draft, not a submission-ready
manuscript. No author details or declarations have been inferred. No external
sites were accessed during this adaptation.

## Template and development conventions

The user supplied `template_reference/PLOS_latex_template.zip`. Its LaTeX
source identifies itself as **version 3.8, April 2026**. `main.tex` retains
the template's required packages, typography/layout conventions, right-side
line numbering, caption styling, and page/total-page footer. Double spacing
is enabled. The root `plos2025.bst` is an unmodified copy from that ZIP;
the existing `references.bib` entries and citation keys are unchanged.

The generic template says a title may have 250 characters and illustrates
a different section order. For this Genetics draft, the requested **200-character
full title**, **70-character short title**, and journal-specific section order
take precedence. The full title is unchanged; the short title is provisional.

Development deliberately retains `\input` section files, visible PDF figures,
and editorial placeholders. The template calls for a **single `.tex` source**
and **no embedded graphics** at submission. Prepare that version separately
after scientific review: inline the section sources, use
`\draftfiguresfalse` to suppress graphics while retaining captions/labels,
and remove editorial TODOs and literal placeholders. Hiding TODO commands
alone does not remove every placeholder. Keep continuous line numbering.

The draft's final Supporting Information Captions section still displays the
four supporting figures with their original captions. Final submission needs
caption-only supporting entries and separately supplied supporting files.
Figure captions and cell-based tables should follow their first textual
citation; existing page floats remain a drafting convenience. Final figure
citations should use “Fig” and supporting item citations should use PLOS
item names. Existing prose has deliberately not been rewritten for these
editorial changes. Three heading levels at most are permitted; the draft
uses only section and subsection headings.

## Information for the manuscript and supporting files

- **Title page:** confirm title/short title; supply author names, order,
  affiliations, corresponding author email, and any applicable author notes.
- **Abstract:** replace the existing placeholder; unstructured, at most
  300 words, without references.
- **Author Summary:** replace the placeholder with 150–200 nontechnical
  words, in the first person. A compliant summary has not been authored here.
- **Scientific text:** complete the existing Introduction, Materials and
  Methods, Discussion, and scientific TODOs after review. Existing Results
  and all figure-caption text have been preserved.
- **Acknowledgments:** supply actual agreed acknowledgments, if applicable;
  this section must not be used for funding or submission declarations.
- **References:** complete missing resource, software, GWAS, and other source
  citations. Check complete author metadata (several existing entries end
  after three named authors plus `others`) and journal abbreviations before
  submission. Citation order and numerical Vancouver formatting now use the
  official style, but incomplete metadata has not been invented.
- **STROBE-MR:** prepare the completed reporting checklist as supporting
  information, including manuscript page numbers and relevant text. The
  draft reserves “S1 Checklist”; this is not a completed checklist.
- **Supporting information captions:** finalize S1–S4 Fig captions and
  S1–S2 Table captions, agree the supporting-file arrangement, and verify the
  inherited annotations and foreground description in the motif workbook.
- **Figures:** review panel readability and prepare separate publication
  TIFF/EPS files later. Existing PDF/XLSX assets must remain unchanged.

## Required submission-system information

The local template explicitly directs funding, data availability, and author
contributions to dedicated PLOS Editorial Manager fields, which are incorporated
into the article in production. The following are information placeholders,
not statements for submission:

| Submission field | Information still required |
|---|---|
| Data Availability Statement | **TODO:** author-approved availability statement, accession IDs/DOIs or repository links, and any verified access restrictions and access procedure. Reconcile with the manuscript's data-source description. |
| Author Contributions | **TODO:** author-approved CRediT roles for every author, agreed with the final author list. |
| Funding / Financial Disclosure | **TODO:** verified funders, grant identifiers, recipients, and funder roles. Do not assume no funding or no funder involvement; keep this out of Acknowledgments. |
| Competing Interests | **TODO:** declaration agreed by all authors. Do not assume there are no competing interests. |
| Author metadata | **TODO:** author names, affiliations, correspondence details, ORCID information as required, and consistency with the title page. |

Before upload, review these fields and supporting-file requirements against
the journal's current instructions. Previously established guideline URLs
(recorded here for later review; not fetched during this task):

- https://journals.plos.org/plosgenetics/s/submission-guidelines
- https://journals.plos.org/plosgenetics/s/latex
