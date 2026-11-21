"""Exact finite examples for the kernel, the projected coefficients and the ramp bound.

Run with python3 checks/finite_average.py. The group is C2 x C2 x C4 x C2,
and the translation subgroup is its last C2. Fourth roots of unity are
represented by their exponents. All sums use integer or rational arithmetic.
"""

from fractions import Fraction
from itertools import product


MODULI = (2, 2, 4, 2)
GROUP = tuple(product(*(range(d) for d in MODULI)))
F = len(GROUP) // 2
ZERO = (0, 0, 0, 0)
IOTA = (1, 0, 0, 0)
TRANSLATIONS = (ZERO, (0, 0, 0, 1))
ORDERS = {
    (0, 0, 0, 0): 24,
    (0, 0, 0, 1): 31,
    (0, 0, 1, 0): 30,
    (0, 0, 2, 0): 27,
    (0, 0, 3, 0): 25,
}
REAL = (1, 0, -1, 0)
IMAG = (0, 1, 0, -1)


def add(*elements):
    return tuple(sum(coords) % d for coords, d in zip(zip(*elements), MODULI))


def neg(element):
    return tuple(-x % d for x, d in zip(element, MODULI))


def phase(character, element):
    return sum((4 // d) * a * x for d, a, x in zip(MODULI, character, element)) % 4


def phase_sum(terms):
    re = im = 0
    for exponent, weight in terms:
        re += weight * REAL[exponent % 4]
        im += weight * IMAG[exponent % 4]
    return re, im


def check_kernels():
    count = 0
    for inertia in ((ZERO,), (ZERO, IOTA), TRANSLATIONS):
        # Each distinct Frobenius coset, including both sides of Phi.
        cosets = {frozenset(add(frob, i) for i in inertia) for frob in GROUP}
        for coset in cosets:
            mass = 0
            for g in GROUP:
                actual = phase_sum(
                    (phase(chi, k) - phase(chi, g), Fraction(1, F * len(inertia)))
                    for chi in GROUP for k in coset
                )
                expected = Fraction(2, len(inertia)) if g in coset else 0
                assert actual == (expected, 0)
                mass += expected
                count += 1
            assert mass == 2
            if IOTA in inertia:
                assert sum(g[0] == 1 for g in coset) == 1
    return count


def check_projection():
    count = 0
    order_sum = sum(ORDERS.values())
    # The asymmetric orders on the two odd residues detect conjugation errors.
    assert phase_sum((t[2], a) for t, a in ORDERS.items())[1] != 0
    for thick_set in (TRANSLATIONS[:1], TRANSLATIONS):
        for chi in GROUP:
            actual = phase_sum(
                (-phase(chi, add(IOTA, u, s, neg(t))),
                 Fraction(a, order_sum * len(TRANSLATIONS) * len(thick_set)))
                for t, a in ORDERS.items() for u in TRANSLATIONS for s in thick_set
            )
            expected = (0, 0)
            if chi[3] == 0:
                expected = phase_sum(
                    (phase(chi, IOTA) + phase(chi, t), Fraction(a, order_sum))
                    for t, a in ORDERS.items()
                )
            assert actual == expected
            count += 1
        alpha = Fraction(len(thick_set), len(TRANSLATIONS))
        p_max = Fraction(max(sum(t[2] == q for t in ORDERS) for q in range(4)), len(ORDERS))
        for g in GROUP:
            hits = sum(
                g == add(IOTA, u, s, neg(t))
                for t in ORDERS for u in TRANSLATIONS for s in thick_set
            )
            assert Fraction(hits, len(ORDERS) * len(TRANSLATIONS)) <= alpha * p_max
    # Retain one lift of each power on the cyclic level of order four.
    discarded = [chi for chi in GROUP if chi[3] == 0 and chi[:2] != (0, 0)]
    assert sum(chi[2] == 0 for chi in discarded) == 3
    assert sum(chi[2] != 0 for chi in discarded) == 9
    return count


def check_ramps():
    x, b, h = 7, 4, Fraction(3, 2)
    count = 0
    for a in ORDERS.values():
        v = Fraction(a, x)
        for w in (b, b + Fraction(1, 2), b + h):
            for n in range(1, 15):
                difference = max(Fraction(a, n) - w, 0) - max(Fraction(a, n) - v, 0)
                if w >= v:
                    assert -(h + Fraction(b, x)) * (n <= x) <= difference <= 0
                else:
                    assert difference >= 0
                    assert v / w <= 1 + Fraction(1, x)
                count += 1
    return count


if __name__ == "__main__":
    print(f"{check_kernels()} local character-kernel identities passed")
    print(f"{check_projection()} weighted subgroup projections passed")
    print(f"{check_ramps()} directional ramp comparisons passed")
