# Exact example checks

This directory accompanies *Pure fibres of the four-site matching map*.
It supports the computational statement in the manuscript's “Proof artifacts
and AI use” section: representative tensors and all their site
permutations satisfy the stated equations. The classification and geometry
are proved analytically in the manuscript.

Run with Python 3.9 or newer, using only its standard library:

```text
python check_examples.py
```

The program writes `results.json` beside its source and prints a summary.
The saved `results.json` was written on Windows and has CRLF line endings;
the checksum in `manifest.json` is of that file. A run on Linux writes the
same lines with LF endings, hence a file with a different checksum.
It needs no network, external data, numerical libraries, or historical research
programs. It uses `fractions.Fraction` throughout, with no floating-point
arithmetic. Failed checks raise an exception and produce a nonzero exit status;
the checks remain enabled under `python -O`.

## What is checked

The coordinate order at every site is `(c,a,b)`. Sites are ordered
`(l,m,t,z)`, and the six matrices are ordered
`(L,P,Q,R,S,K)=(W_lm,W_lt,W_lz,W_mt,W_mz,W_tz)`.
Every product uses the same actual matrices. For each example and each of the
24 site permutations, the program checks:

1. All 81 coefficients of `LK + PS + QR = tau * c_l*c_m*c_t*c_z`.
2. Each whole matrix rank and each projected rank after deleting its `c` row
   and `c` column, against explicit expected ranks.
3. The projected support graph and its rank type.
4. Each whole block's zero/nonzero status directly from its entries.
5. Each matching product's zero/nonzero status directly from its 81
   coefficients, including the cases with a whole zero block.

A permutation is recorded as `new_to_old_sites`: its value at a new site is
the old site occupying that role. Reversed endpoint order transposes the
corresponding matrix. No sign or conjugation is introduced.

The empty-support examples receive an additional exact check: extract their
constant edge vector and the eight columns of the four site maps, verify
`q(a)=tau`, verify orthogonality and pairwise isotropy for the polarized form,
and verify that their combined image has rank two. This is repeated after
every site permutation. It checks the example's image rank, not the dimension
of the empty-support stratum.

## Examples and manuscript correspondence

All unscaled examples use `tau=3`. One scaled example uses `tau=12`.
The source includes 18 original representative and boundary witnesses, the
literal triangle displayed in the manuscript, and one additional scaling of
the manuscript's empty-support example.

| Projected type | Number of examples | Relevant atlas statement |
|---|---:|---|
| Empty | 2 | `atlas:empty` |
| Single rank-one edge | 2 | `atlas:single` |
| Single rank-two edge | 2 | `atlas:single` |
| Wedge `(1,1)` | 6 | `atlas:wedge11`, all five families |
| Wedge `(2,1)` | 2 | `atlas:wedge21`, both branches |
| Wedge `(2,2)` | 1 | `atlas:wedge22` |
| Triangle `(1,1,1)` | 2 | `atlas:triangle`, including a whole-rank drop |
| Triangle `(2,1,1)` | 2 | `atlas:triangle`, `atlas:triangle-ranks` |
| Four-cycle `(1,1,1,1)` | 1 | `atlas:cycle` |

`paper_triangle211_kappa1` is the literal triangle from “Examples and scope”
at `tau=3`; its whole ranks are `(3,2,1,2,1,1)`, and all matching products
are nonzero. `allzero_dimW2` is the homogenization of the six affine edge
forms printed there. `paper_empty_scaled_s2` multiplies each of those six
matrices by `s=2`, giving target `3*s^2=12`.

Other witnesses include both zero and nonzero opposite-edge scalar branches,
free whole rank-three matrices, all five low-low wedge families, zero and
nonzero parameters in family V, correlated whole `c` components, and a
rank-one triangle whose opposite whole matrix has rank one.

Each result records all six matrices as rational strings, their expected
rank pattern, their activity, and all 24 permutation records. Rational strings
such as `3/2` are exact, not decimal approximations.

## Fresh result

The included run passed:

- 20 examples and all 9 projected support/rank types;
- 480 labelled tuples, comprising all 24 site permutations of each example;
- 38,880 exact tensor-coefficient equalities;
- 2,880 whole-rank checks and 2,880 projected-rank checks;
- 2,880 whole-block activity checks and 1,440 matching-product activity checks.

`results.json` includes the SHA-256 of the exact checker source used.
`manifest.json` records the hashes of the source, this README, and the fresh
result. The output contains no timestamps and is deterministic for unchanged
source. Re-running overwrites only `results.json`.

## Limits

These are finite rational examples inside the complex fibre. They do not
check every parameter value or every subcase of every atlas chart. They do
not establish exhaustion of the nine types, necessity of a normal form,
dimensions, irreducible components, a scheme structure, source realization,
novelty, or formal verification. The scaling example tests one rational
scaling; the manuscript's argument for arbitrary nonzero complex `tau` is
analytic. The number of checked coefficients is an accounting total and is
not a measure of theorem coverage.
