"""Reproducible exact certificate: sharp m=29 result for the two-periodic
GF(64) parity-check family correcting either one information-symbol
substitution or one adjacent information-symbol swap.

All arithmetic uses integers, bit operations and Euclid's algorithm in F_2[x].
No packages or external services are needed.
"""

from collections import Counter

MODULUS = 0x43  # X^6 + X + 1
Q_MINUS_1 = 63
# Reduction mod 2 of the 63rd cyclotomic polynomial.
PHI_63 = sum(1 << i for i in (36, 33, 27, 24, 18, 12, 9, 3, 0))


def gf64_mul(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        a <<= 1
        if a & 64:
            a ^= MODULUS
        b >>= 1
    return out


def gf64_pow(a, exponent):
    out = 1
    for _ in range(exponent):
        out = gf64_mul(out, a)
    return out


def discrete_log_table(alpha):
    table = {}
    v = 1
    for i in range(63):
        table[v] = i
        v = gf64_mul(v, alpha)
    return table


def gf2_poly_rem(a, b):
    assert b
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def gf2_poly_gcd(a, b):
    while b:
        a, b = b, gf2_poly_rem(a, b)
    return a


def collision_free(b, D, m):
    """b=log_alpha(beta), D=log_alpha((1+alpha)/(1+beta))."""
    vertex_slopes = [(i - (b if i % 2 == 0 else 0)) % 63
                     for i in range(1, m + 1)]
    edge_slopes = [(i + D) % 63 for i in range(1, m)]
    values = vertex_slopes + edge_slopes
    return len(values) == len(set(values))


def all_formal_candidates(m):
    return [(b, D) for b in range(2, 63) for D in range(63)
            if collision_free(b, D, m)]


def main():
    candidates = all_formal_candidates(30)
    expected = ([(2, d) for d in range(29, 34)]
                + [(4, d) for d in range(29, 32)]
                + [(6, 29), (59, 34)]
                + [(61, d) for d in range(32, 35)])
    assert candidates == sorted(expected)
    print('All combinatorial candidates for m=30:', candidates)
    print('Candidate count:', len(candidates))

    for b, D in candidates:
        # In characteristic two, alpha**D*(1+alpha**b)=1+alpha.
        f = (1 << (D + b)) ^ (1 << D) ^ (1 << 1) ^ 1
        assert gf2_poly_gcd(PHI_63, f) == 1
    print('All 13 cyclotomic gcds equal 1: NO primitive alpha can reach m=30.')

    primitives = [a for a in range(2, 64)
                  if len(discrete_log_table(a)) == Q_MINUS_1]
    assert len(primitives) == 36
    distribution = Counter()
    winners = []
    for alpha in primitives:
        logs = discrete_log_table(alpha)
        for beta in range(2, 64):
            b = logs[beta]
            D = (logs[alpha ^ 1] - logs[beta ^ 1]) % 63
            best = max((m for m in range(2, 34)
                        if (b != 1 if m % 2 == 0 else b != 62)
                        and collision_free(b, D, m)), default=0)
            distribution[best] += 1
            if best == 29:
                winners.append((alpha, beta, b, D))
    assert len(primitives) * 62 == sum(distribution.values()) == 2232
    assert max(distribution) == 29
    assert len(winners) == 12
    print('Primitive alpha values:', len(primitives))
    print('Primitive-alpha / beta != 1 parameter pairs evaluated:', sum(distribution.values()))
    print('Maximum m across all pairs:', max(distribution))
    print('Number of maximizing parameter pairs:', len(winners))

    alpha, beta = 7, 51
    logs = discrete_log_table(alpha)
    assert len(logs) == 63
    assert gf64_pow(alpha, 6) == beta
    b = logs[beta]
    D = (logs[alpha ^ 1] - logs[beta ^ 1]) % 63
    assert (b, D) == (6, 30)
    assert collision_free(b, D, 29)
    assert not collision_free(b, D, 30)
    assert gf64_mul(gf64_pow(alpha, D), beta ^ 1) == alpha ^ 1
    print('Explicit maximizing example: alpha=7, beta=alpha^6=51, b=6, D=30')
    print('Finite-field representation: GF(2)[x]/(x^6+x+1); alpha=x^2+x+1')
    print('The minimal polynomial of alpha is x^6+x^4+x^3+x+1.')
    print('29 information symbols + 2 check symbols: success (odd m).')
    print('28 information symbols + 2 check symbols: success (even m).')
    print('30 information symbols + 2 check symbols: impossible in this family.')


if __name__ == '__main__':
    main()
