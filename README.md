# Pure fibres of the four-site matching map

Preprint by Ignoramus, revision 2, 11 November 2026.

**Version of record:** Zenodo, https://doi.org/10.5281/zenodo.24138298 . This repository mirrors that record so that the source and the checking program can be browsed and errors can be reported in the issue tracker.

## What the paper shows

Let V_1, ..., V_4 be three-dimensional complex vector spaces, fix nonzero c_i in V_i and a nonzero complex number tau, and let F_tau be the set of tuples of six pair tensors W_ij in V_i (x) V_j with W_12 W_34 + W_13 W_24 + W_14 W_23 = tau c_1 c_2 c_3 c_4 (the fibre of the four-site matching map over a pure tensor).

- F_tau has dimension 19 and exactly twelve irreducible components of dimension 19, each isomorphic to A^18 x G_m (Theorem 1.1, part 1).
- The open subset of F_tau where all six W_ij are nonzero has dimension 15 and exactly twelve irreducible components of dimension 15 (Theorem 1.1, part 2).
- The space of three ordered nonzero summands with sum tau c_1 c_2 c_3 c_4, one decomposable across each of 12|34, 13|24 and 14|23, has dimension 12 and exactly twelve irreducible components of dimension 12 (Theorem 1.1, part 3).
- Up to simultaneous relabelling of the four sites, every point of F_tau has exactly one of nine projected support and rank patterns, every pattern occurs, and Section 3 parametrizes all nine (Theorem 1.2).

Theorem 1.1 concerns the components of maximal dimension only; lower-dimensional components are not classified, and nothing is asserted for the zero fibre, for other local dimensions or for scheme multiplicities. The closest earlier results are credited in Section 1 and in the references of the paper.

## Contents

| Path | What it is |
|---|---|
| `paper.pdf` | The manuscript (20 pages) |
| `paper.tex` | Its LaTeX source |
| `atlas.tex` | The source of Section 3, input by `paper.tex` |
| `reproducibility/check_examples.py` | The checking program |
| `reproducibility/results.json` | Its saved output |
| `reproducibility/manifest.json` | Sizes and SHA-256 checksums of the other three files, and the Python version the program was tested with |
| `reproducibility/README.md` | How to run the program and what it does and does not establish |
| `FILES.sha256` | SHA-256 checksums of the other seven files |

## How to run the check

In the folder `reproducibility/`, run `python check_examples.py`. It uses only the Python standard library and writes `results.json` beside itself. In exact rational arithmetic it checks 20 examples with rational entries, covering all nine projected support and rank types, each under all 24 permutations of the four sites: for each, all 81 coefficients of LK + PS + QR = tau c_l c_m c_t c_z, the whole and projected ranks of the six blocks, the projected support type, and which blocks and which matching products are zero; for the empty-support examples, also that q(a) = tau and that the sum of the images of the maps B_i is totally isotropic, orthogonal to a and of dimension two, in the notation of Proposition 3.1. These are finite checks of examples. They do not prove that the nine types are exhaustive or that the parametrizations are complete, and they do not check the dimension or component counts.

The saved `results.json` was written on Windows and has CRLF line endings, and the checksum in `manifest.json` is of that file; a run on Linux writes the same lines with LF endings, hence a file with a different checksum.

## Disclosure

The author directed the project, takes responsibility for the paper, and used OpenAI Codex to produce the results, proofs, checking program and text, and Anthropic Claude to check and revise the paper; the paper's closing section "Proof artifacts and AI use" gives the details. The attribution to OpenAI Codex rests on the project's records and the author's account; which model versions were used at each step has not been verified. The author has not verified the mathematics; no human has checked the proofs, and no human peer review or formal proof-assistant verification is claimed. No AI system is an author. The author publishes under a pseudonym. All checks described in that section were made by AI agents; none of them is journal peer review, a formal verification or a certificate of originality or priority.

## Corrections

Please open an issue, or write to ignoramus110926@gmail.com.

## Licence

Creative Commons Attribution 4.0 International (CC BY 4.0), https://creativecommons.org/licenses/by/4.0/
