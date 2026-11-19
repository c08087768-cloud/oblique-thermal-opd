# Applied Optics submission format

The English main manuscript uses the publisher's current `optica-article` universal submission template, with `opticajournal` and `Research Article` settings. The publisher's template page, updated 3 November 2026, lists this for Applied Optics. It is a single-column submission format; the final Applied Optics journal is typeset in two columns during production. This file does not present the draft as a published article.

The supplement uses the separate `osa-supplemental-document` class. Per its instructions, the author block and copyright are blank; the publisher adds the supplemental coversheet during production. All existing supplemental sections, tables, figures and S-prefixed labels are retained.

English: run `pdflatex -jobname=main_revised main.tex` and `pdflatex supplement.tex` until references are stable. Chinese review copies: use `xelatex -jobname=main_zh_revised main_zh.tex` and `xelatex supplement_zh.tex`. The latest main PDFs are `main_revised.pdf` and `main_zh_revised.pdf`; the earlier `main.pdf` and `main_zh.pdf` were retained because the English PDF was open in a reader and could not be replaced. The Chinese copies keep the same scientific text with a CJK font compatibility layer; they are not the journal submission language. This Windows build uses SimSun, SimHei, KaiTi, FangSong and Microsoft YaHei for Chinese, with TeX Gyre Heros for Latin sans-serif.

Author affiliations, superscripts and corresponding email were supplied in the user's reference image on 20 November 2026. The scientific title and abstract were retained from the current approved manuscript; text visible in that reference image was not substituted. Existing manual reference entries and citation keys are retained.

The user subsequently confirmed that the authors have no conflicts of interest and that GPT-6.1 sol was used for manuscript translation and code writing. These statements, a Data availability section and the publisher's Supplemental document callout were added to the English and Chinese main manuscripts. Scientific content and the supplements were not changed.

The supplied code repository is https://github.com/c08087768-cloud/oblique-thermal-opd . GitHub repository metadata was checked on 20 November 2026 and reports private visibility. The user intends to prepare a public release and requested that the manuscript not claim the code is already public. The statement records this current status, the selected derived data in Supplement 1, and the absence of the complete archived FEM nodal fields from the repository. No repository visibility change or promise to provide unpublished data has been made. Public release and access arrangements for the remaining underlying data need to be finalized before submission. No funding information has been supplied.

## Template sources and license

- Current publisher guidance: https://opg.optica.org/content/author/portal/item/templates-default
- Article archive, downloaded from publisher: https://opg.optica.org/resources/author/templates/latex-universal-template.zip
- Official supplemental template instructions: https://www.overleaf.com/latex/templates/latex-template-for-preparing-supplementary-material-for-submission-to-optica-publishing-groups-journals/tzgytvsqrftn
- The publisher's supplemental archive returned a CAPTCHA page in this session. The unmodified class was therefore retrieved from a fixed public mirror: https://github.com/shankhya/colour-lifting-repo/blob/f56ab4a7e7816f81fd0d93fb72f34ddb4afec449/manuscript/osa-supplemental-document.cls . Its published checksum and third-party notice were checked; the publisher-downloaded article class also matches that mirror's article class checksum. Direct byte comparison against a newly downloaded publisher supplemental ZIP was not possible.

The bundled publisher class/style/BibTeX support files retain their upstream copyright and LaTeX Project Public License notices. They have not been edited. Narrow `.gitattributes` rules preserve their original bytes across checkouts. Exact hashes and retrieval provenance are recorded in `template_conversion_manifest.json`.

## Export and validation

The repository retains its existing modular LaTeX layout. The delivery includes separate English main and supplemental source archives, each with one `.tex` file and all `\input` fragments expanded, following the publisher's submission guidance. Figure and numerical data are unchanged. Compilation verifies document construction and references; it does not establish scientific correctness.

The standalone main produces identical extracted PDF text to the modular build. Inlining the supplement can change LaTeX float placement; its text/number inventory and seven-page output were checked, and all seven exported pages were visually inspected.
