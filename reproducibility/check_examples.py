"""Exact, finite example checks for Pure fibres of the four-site matching map.

Python standard library only. No inherited research executable is imported.
Run: python check_examples.py
This checks explicit rational points, not necessity, completeness, dimensions,
irreducible components, or the validity/novelty of the manuscript's proofs.
"""

from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
EDGE_NAMES = ("L", "P", "Q", "R", "S", "K")
MATCHINGS = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))
COORDINATES = tuple(product(range(3), repeat=4))
PERMUTATIONS = tuple(permutations(range(4)))
TYPE_NAMES = (
    "empty", "single1", "single2", "wedge11", "wedge21", "wedge22",
    "triangle111", "triangle211", "cycle1111",
)
PROJECTED_RANKS = (
    (0, 0, 0, 0, 0, 0), (1, 0, 0, 0, 0, 0),
    (2, 0, 0, 0, 0, 0), (1, 1, 0, 0, 0, 0),
    (2, 1, 0, 0, 0, 0), (2, 2, 0, 0, 0, 0),
    (1, 1, 0, 1, 0, 0), (2, 1, 0, 1, 0, 0),
    (0, 1, 1, 1, 1, 0),
)
c, a, b = (F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(1))
TAU = F(3)


def require(condition, detail):
    """Checks remain enabled even when Python is run with -O."""
    if not condition:
        raise AssertionError(detail)


def vec(*terms):
    return tuple(sum((F(k) * v[i] for k, v in terms), F(0)) for i in range(3))


def tensor(*terms):
    return tuple(tuple(sum((F(k) * u[i] * v[j] for k, u, v in terms), F(0))
                       for j in range(3)) for i in range(3))


def outer(u, v):
    return tensor((1, u, v))


def scale(s, matrix):
    return tuple(tuple(F(s) * x for x in row) for row in matrix)


def transpose(matrix):
    return tuple(zip(*matrix))


def rank(matrix):
    """Gaussian elimination over Q, with a proved nonzero pivot each time."""
    rows = [list(map(F, row)) for row in matrix]
    if not rows:
        return 0
    pivot_row = 0
    for column in range(len(rows[0])):
        selected = next((i for i in range(pivot_row, len(rows)) if rows[i][column]), None)
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        pivot = rows[pivot_row][column]
        rows[pivot_row] = [x / pivot for x in rows[pivot_row]]
        for i in range(len(rows)):
            if i != pivot_row:
                factor = rows[i][column]
                rows[i] = [x - factor * y for x, y in zip(rows[i], rows[pivot_row])]
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


def quotient_rank(matrix):
    # Coordinate order is (c,a,b); quotienting by c deletes row/column zero.
    return rank([row[1:] for row in matrix[1:]])


def identify(projected_ranks):
    live = [e for e in EDGES if projected_ranks[e]]
    degree = [sum(i in edge for edge in live) for i in range(4)]
    ranks = tuple(sorted(projected_ranks[e] for e in live))
    if not live:
        return 0
    if len(live) == 1:
        require(ranks in ((1,), (2,)), ("invalid single rank", ranks))
        return 1 if ranks == (1,) else 2
    if len(live) == 2:
        require(sorted(degree) == [0, 1, 1, 2], ("not a wedge", live))
        return {(1, 1): 3, (1, 2): 4, (2, 2): 5}[ranks]
    if len(live) == 3:
        require(sorted(degree) == [0, 2, 2, 2], ("not a triangle", live))
        return {(1, 1, 1): 6, (1, 1, 2): 7}[ranks]
    require(len(live) == 4 and degree == [2, 2, 2, 2] and ranks == (1, 1, 1, 1),
            ("not a rank-one four-cycle", live, degree, ranks))
    return 8


Z, cc = tensor((0, c, c)), outer(c, c)
EXAMPLES = []


def add(name, kind, blocks, whole_ranks, atlas_reference, tau=TAU, manuscript_example=False):
    require(len(blocks) == len(whole_ranks) == 6, name)
    EXAMPLES.append({"name": name, "type": kind, "blocks": dict(zip(EDGES, blocks)),
                     "expected_whole_ranks": tuple(whole_ranks), "tau": F(tau),
                     "atlas_reference": atlas_reference, "manuscript_example": manuscript_example})


# The 18 witnesses from the fresh source-review checker are retained literally.
empty_blocks = [outer(c, vec((1, c), (-2, a))), outer(vec((1, c), (1, a)), c),
                outer(vec((1, c), (-1, a)), c), outer(vec((1, c), (1, a)), c),
                outer(vec((1, c), (1, a)), c), cc]
add("allzero_dimW2", 0, empty_blocks, (1, 1, 1, 1, 1, 1),
    "atlas:empty; paper Examples and scope, affine edge forms", manuscript_example=True)
add("single1_nonzero_kappa", 1,
    [tensor((2, c, c), (-1, a, a)), outer(a, c), cc, cc, outer(a, c), cc],
    (2, 1, 1, 1, 1, 1), "atlas:single, nonzero kappa, both products nonzero")
add("single1_zero_kappa_rank3_free_L", 1,
    [tensor((1, a, a), (1, c, b), (1, b, c)), Z, cc, tensor((TAU, c, c)), cc, Z],
    (3, 0, 1, 1, 1, 0), "atlas:single, zero kappa; atlas:crossing, PS=0")
add("single2_nonzero_kappa", 2,
    [tensor((TAU, c, c), (1, a, a), (1, b, b)), outer(a, c), outer(b, c),
     tensor((-1, b, c)), tensor((-1, a, c)), cc],
    (3, 1, 1, 1, 1, 1), "atlas:single, projected rank two, nonzero kappa")
add("single2_zero_kappa", 2,
    [tensor((1, a, a), (1, b, b)), Z, cc, tensor((TAU, c, c)), cc, Z],
    (2, 0, 1, 1, 1, 0), "atlas:single, projected rank two, zero kappa")
add("wedge11_I", 3,
    [tensor((1, a, a), (1, c, b), (1, b, c)), outer(b, a), cc, tensor((TAU, c, c)), Z, Z],
    (3, 1, 1, 1, 0, 0), "atlas:wedge11, I")
add("wedge11_II", 3,
    [tensor((TAU, c, c), (-1, a, a)), outer(b, a), outer(a, c), outer(a, c), Z, cc],
    (2, 1, 1, 1, 0, 1), "atlas:wedge11, II")
add("wedge11_III", 3,
    [outer(b, a), tensor((TAU, c, c), (-1, a, a)), outer(a, c), outer(c, a), cc, Z],
    (1, 2, 1, 1, 1, 0), "atlas:wedge11, III")
add("wedge11_IV", 3,
    [outer(a, a), tensor((1, a, a), (TAU, c, c)), outer(a, c),
     tensor((-1, a, c), (-1, c, a)), cc, cc],
    (1, 2, 1, 2, 1, 1), "atlas:wedge11, IV")
v, s, k = vec((1, a), (2, c)), vec((1, a), (1, c)), vec((1, b), (-3, c))
add("wedge11_V_nonzero_alpha_beta", 3,
    [tensor((1, v, s), (2, c, c)), tensor((-1, v, k), (-1, c, c)), tensor((2, c, c)),
     tensor((TAU / 2, c, c), (-1, c, k), (F(1, 2), s, c)), outer(s, c), outer(k, c)],
    (2, 2, 1, 2, 1, 1), "atlas:wedge11-V, alpha=2, beta=-1, eta=2")
add("wedge11_V_zero_alpha_beta", 3,
    [outer(a, a), tensor((-1, a, a)), cc, tensor((TAU, c, c)), outer(a, c), outer(a, c)],
    (1, 1, 1, 1, 1, 1), "atlas:wedge11-V, alpha=beta=0")
add("wedge21_B1", 4,
    [tensor((1, a, a), (1, b, b)), outer(a, a), cc, tensor((TAU, c, c)), Z, Z],
    (2, 1, 1, 1, 0, 0), "atlas:wedge21, S=0")
add("wedge21_B2", 4,
    [tensor((1, a, a), (1, b, b)), tensor((TAU, c, c), (1, a, a)),
     tensor((-1, a, c)), outer(c, a), cc, Z],
    (2, 2, 1, 1, 1, 0), "atlas:wedge21, sigma=1")
add("wedge22", 5,
    [tensor((1, a, a), (1, b, b)), tensor((1, a, a), (1, b, b)), cc,
     tensor((TAU, c, c)), Z, Z], (2, 2, 1, 1, 0, 0), "atlas:wedge22")
add("triangle111_two_commonlines", 6,
    [outer(a, a), tensor((-1, a, a), (1, c, b)), tensor((2, c, c)),
     tensor((TAU / 2, c, c), (-F(1, 2), a, b)), outer(a, c), outer(a, c)],
    (1, 2, 1, 2, 1, 1), "atlas:triangle, rank111 after pure-spoke relabelling")
add("triangle111_opposite_wholerank1", 6,
    [tensor((1, a, a), (1, c, c), (1, c, a)), tensor((-1, a, a), (-TAU, c, c)),
     tensor((2, c, c)), tensor((TAU / 2, c, c), (-F(1, 2), c, a),
                             (-F(1, 2), a, a), (TAU / 2, a, c)),
     outer(a, c), outer(a, c)],
    (2, 2, 1, 1, 1, 1), "atlas:triangle-ranks, rank111 whole-rank drop")
add("triangle211", 7,
    [tensor((TAU / 2, c, c), (-F(1, 2), b, a), (-F(1, 2), a, b)),
     tensor((1, a, a), (1, b, c)), outer(a, c), tensor((-1, a, a), (1, b, c)),
     outer(a, c), tensor((2, c, c))],
    (3, 2, 1, 2, 1, 1), "atlas:triangle-form, kappa=2")
u = [vec((1, a), (i + 1, c), (i, b)) for i in range(4)]
add("fourcycle_arbitrary_lifts", 8,
    [tensor((2, c, c)), outer(u[0], u[2]), outer(u[0], u[3]), outer(u[1], u[2]),
     tensor((-1, u[1], u[3])), tensor((TAU / 2, c, c))],
    (1, 1, 1, 1, 1, 1), "atlas:cycle-form, gamma=2, mixed whole lifts")

# Match both literal examples in paper.tex and check one nontrivial scaling.
add("paper_triangle211_kappa1", 7,
    [tensor((TAU, c, c), (-1, b, a), (-1, a, b)), tensor((1, a, a), (1, b, c)),
     outer(a, c), tensor((-1, a, a), (1, b, c)), outer(a, c), cc],
    (3, 2, 1, 2, 1, 1), "paper Examples and scope, literal triangle at tau=3",
    manuscript_example=True)
add("paper_empty_scaled_s2", 0, [scale(2, m) for m in empty_blocks],
    (1, 1, 1, 1, 1, 1), "paper empty-support scaling, s=2 and tau=12", tau=12,
    manuscript_example=True)


def relabel(blocks, permutation):
    """permutation[new_site] = old_site; reversed endpoint order transposes."""
    result = {}
    for i, j in EDGES:
        old_i, old_j = permutation[i], permutation[j]
        result[i, j] = (blocks[old_i, old_j] if old_i < old_j
                        else transpose(blocks[old_j, old_i]))
    return result


def product_coefficients(blocks, matching):
    (i, j), (k, n) = matching
    return tuple(blocks[i, j][cs[i]][cs[j]] * blocks[k, n][cs[k]][cs[n]]
                 for cs in COORDINATES)


def q(x):
    return x[0] * x[5] + x[1] * x[4] + x[2] * x[3]


def polar(x, y):
    return q(tuple(u + v for u, v in zip(x, y))) - q(x) - q(y)


def empty_support_data(blocks, tau):
    """Extract a and all eight site-map columns from the actual six matrices."""
    constant = tuple(blocks[e][0][0] for e in EDGES)
    columns = []
    for site in range(4):
        for coordinate in (1, 2):
            column = []
            for i, j in EDGES:
                value = blocks[i, j][coordinate][0] if site == i else (
                    blocks[i, j][0][coordinate] if site == j else F(0))
                column.append(value)
            columns.append(tuple(column))
    require(q(constant) == tau, "empty-support q(a)")
    require(all(polar(constant, col) == 0 for col in columns), "a orthogonal to image")
    require(all(polar(x, y) == 0 for x in columns for y in columns), "image isotropic")
    dimension = rank(columns)
    require(dimension == 2, ("empty-support image dimension", dimension))
    return {"constant_edge_coefficients": list(map(str, constant)),
            "site_map_columns": [list(map(str, col)) for col in columns],
            "image_rank": dimension, "orthogonal_and_isotropic": True}


def main():
    require(len(PERMUTATIONS) == 24 and len(set(PERMUTATIONS)) == 24, "permutations")
    require(len(COORDINATES) == 81, "tensor coordinates")
    records = []
    seen_types = set()
    for example in EXAMPLES:
        kind, original, tau = example["type"], example["blocks"], example["tau"]
        expected_quotient = dict(zip(EDGES, PROJECTED_RANKS[kind]))
        expected_whole = dict(zip(EDGES, example["expected_whole_ranks"]))
        permutation_records = []
        for permutation in PERMUTATIONS:
            blocks = relabel(original, permutation)
            quotient_ranks = {e: quotient_rank(blocks[e]) for e in EDGES}
            whole_ranks = {e: rank(blocks[e]) for e in EDGES}
            require(identify(quotient_ranks) == kind, (example["name"], permutation, "type"))
            for e in EDGES:
                old_edge = tuple(sorted(permutation[i] for i in e))
                require(quotient_ranks[e] == expected_quotient[old_edge],
                        (example["name"], permutation, e, "quotient rank"))
                require(whole_ranks[e] == expected_whole[old_edge],
                        (example["name"], permutation, e, "whole rank"))
                actual_nonzero = any(value != 0 for row in blocks[e] for value in row)
                require(actual_nonzero == (expected_whole[old_edge] > 0),
                        (example["name"], permutation, e, "whole activity"))
            products = [product_coefficients(blocks, matching) for matching in MATCHINGS]
            for index, cs in enumerate(COORDINATES):
                total = sum((values[index] for values in products), F(0))
                target = tau if cs == (0, 0, 0, 0) else F(0)
                require(total == target, (example["name"], permutation, cs, total, target))
            product_activity = [any(values) for values in products]
            expected_activity = [all(expected_whole[tuple(sorted(permutation[i] for i in edge))] > 0
                                     for edge in matching) for matching in MATCHINGS]
            require(product_activity == expected_activity,
                    (example["name"], permutation, "matching-product activity"))
            if kind == 0:
                empty_support_data(blocks, tau)
            permutation_records.append({"new_to_old_sites": list(permutation),
                "quotient_ranks": [quotient_ranks[e] for e in EDGES],
                "whole_ranks": [whole_ranks[e] for e in EDGES],
                "whole_activity": [whole_ranks[e] > 0 for e in EDGES],
                "matching_product_activity": product_activity,
                "all_six_whole_blocks_nonzero": all(whole_ranks.values()),
                "all_three_matching_products_nonzero": all(product_activity),
                "target_coefficient": str(tau), "tensor_coefficients_checked": 81})
        base = permutation_records[0]
        record = {"name": example["name"], "type": TYPE_NAMES[kind], "type_number": kind,
                  "tau": str(tau), "atlas_reference": example["atlas_reference"],
                  "manuscript_example": example["manuscript_example"],
                  "blocks_in_c_a_b_coordinates": {name: [[str(x) for x in row] for row in original[e]]
                                                 for name, e in zip(EDGE_NAMES, EDGES)},
                  "quotient_ranks": base["quotient_ranks"], "whole_ranks": base["whole_ranks"],
                  "whole_activity": base["whole_activity"],
                  "matching_product_activity": base["matching_product_activity"],
                  "permutations": permutation_records, "coefficient_checks": 24 * 81}
        if kind == 0:
            record["empty_support_incidence_check"] = empty_support_data(original, tau)
        records.append(record)
        seen_types.add(kind)
    require(seen_types == set(range(9)), ("missing support types", seen_types))
    result = {"status": "PASS", "arithmetic": "fractions.Fraction over Q, no floating point",
        "basis_order": ["c", "a", "b"], "site_order": ["l", "m", "t", "z"],
        "edge_order": list(EDGE_NAMES), "matching_order": ["LK", "PS", "QR"],
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "number_examples": len(records), "support_types_checked": list(TYPE_NAMES),
        "permutations_per_example": 24, "labelled_tuples_checked": 24 * len(records),
        "tensor_coefficient_equalities": 81 * 24 * len(records),
        "whole_rank_assertions": 6 * 24 * len(records),
        "quotient_rank_assertions": 6 * 24 * len(records),
        "whole_activity_assertions": 6 * 24 * len(records),
        "matching_product_activity_assertions": 3 * 24 * len(records),
        "limitations": ["Only these finite rational examples and their site permutations are checked.",
            "No exhaustive search, generic-parameter verification, completeness proof, or formal proof.",
            "No dimension, irreducible-component, scheme, source-realization, or novelty verification.",
            "The empty-support image rank is checked for these points, not the dimension of a stratum."],
        "examples": records}
    output = HERE / "results.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf8")
    print(json.dumps({key: value for key, value in result.items() if key != "examples"}, indent=2))


if __name__ == "__main__":
    main()
