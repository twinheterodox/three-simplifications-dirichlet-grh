"""Exact normalization and conjugation checks on one finite abelian host."""
from fractions import Fraction
from itertools import product


MODULI = (2, 2, 4, 2)
GROUP = tuple(product(*(range(d) for d in MODULI)))
ZERO = (0, 0, 0, 0)
IOTA = (1, 0, 0, 0)
A0 = (ZERO, (0, 0, 0, 1))
THICK = (A0[1],)
TOP_WEIGHTS = {
    (0, 0, 0, 0): 97,
    (0, 0, 0, 1): 103,
    (0, 0, 1, 0): 101,
    (0, 0, 2, 1): 95,
    (0, 0, 3, 1): 107,
}
FOURTH_ROOTS = ((1, 0), (0, 1), (-1, 0), (0, -1))


def add(*points):
    return tuple(sum(coords) % m for coords, m in zip(zip(*points), MODULI))


def neg(point):
    return tuple(-g % m for g, m in zip(point, MODULI))


def exponent(char, point):
    """The value of char at point is i to this power."""
    return sum((4 // m) * a * g for m, a, g in zip(MODULI, char, point)) % 4


def gaussian_sum(terms):
    real = imag = 0
    for weight, power in terms:
        r, i = FOURTH_ROOTS[power % 4]
        real += weight * r
        imag += weight * i
    return real, imag


def check_local_kernel(coset):
    for g in GROUP:
        result = gaussian_sum((1, exponent(char, z) - exponent(char, g))
                              for char in GROUP for z in coset)
        assert result == (len(GROUP) * (g in coset), 0)
    for t in TOP_WEIGHTS:
        hits = sum(add(IOTA, s, u, neg(t)) in coset for u in A0 for s in THICK)
        target = {add(IOTA, neg(t), u) for u in A0}
        assert hits == len(THICK) * len(target & coset)


def check_projector():
    for char in GROUP:
        actual = gaussian_sum((weight, -exponent(char, add(IOTA, s, u, neg(t))))
                              for t, weight in TOP_WEIGHTS.items()
                              for u in A0 for s in THICK)
        if all(exponent(char, u) == 0 for u in A0):
            expected = gaussian_sum((len(A0) * len(THICK) * weight,
                                     exponent(char, IOTA) + exponent(char, t))
                                    for t, weight in TOP_WEIGHTS.items())
        else:
            expected = (0, 0)
        assert actual == expected
    # The nonsymmetric top distribution detects reversal of the top phase.
    assert gaussian_sum((w, t[2]) for t, w in TOP_WEIGHTS.items())[1] != 0
    for j in range(4):
        ideal = gaussian_sum((Fraction(3 + (k == 0), 13), j * k) for k in range(4))
        assert ideal == (1 if j == 0 else Fraction(1, 13), 0)


if __name__ == "__main__":
    # Ordinary Frobenius, distinguished ramification, and a split ramified prime.
    local_data = (
        ((1, 0, 1, 1), (ZERO,)),
        ((0, 0, 3, 0), (ZERO, IOTA)),
        ((0, 1, 1, 0), (ZERO, (0, 1, 0, 0))),
    )
    for frob, inertia in local_data:
        for base in (frob, add(frob, frob)):
            check_local_kernel({add(base, h) for h in inertia})
    check_projector()
    print("Passed 192 local kernels, 30 translate counts, 32 complex Fourier "
          "coefficients and 4 ideal shares, all in exact arithmetic.")
