#!/usr/bin/env python3
"""
Exact-arithmetic verification for
  "Adjacent-Transposition Syndromes and Sharp Length Bounds for
   Two-Periodic Mixed-Error Codes over Finite Fields"

Everything is computed with exact binary-polynomial arithmetic; no floating
point or randomization is used (except in the threshold check of part 6,
which only evaluates the explicit real-valued function J(q)).

Parts
  0. GF(7): the worked example of Appendix A (general locator for odd and
     even positions, failure of the characteristic-two formula, and the
     collision showing that n - 2 <= L cannot be dropped).
  1. GF(64): achievability (alpha = t^2+t+1, beta = alpha^6, m = 29) and the
     13 residue pairs with gcd(Phi_63, x^{D+b}+x^D+x+1) = 1 (m = 30 impossible).
  2. GF(64): the fixed parameters of Beaugris (2025) (alpha = t, beta = alpha^2)
     admit exactly m = 6 under joint decoding; with alpha = t kept and beta
     optimized the maximum is 22; jointly optimizing alpha and beta gives 29.
  3. GF(64): the explicit 29-symbol encoding, the two decoded errors, and the
     codon-aligned DNA representation.
  4. Exhaustive optima for 2 <= k <= K (default K = 10) over all primitive
     alpha (Frobenius-orbit representatives) and all beta != 1:
       m_max(q)      -- check symbols error-free (Table 1 of the paper);
       m_max^all(q)  -- errors allowed in every coordinate (Table 2).
  5. Cross-check of the projective criterion (Lemma 4.1 / Corollary 4.2)
     against brute-force syndrome injectivity for small k.
  6. The threshold k >= 29 in the near-half lower bound.
  7. Numerical check of the Jacobi-sum expansion (30) of the counting
     function T, with the conjugated coefficients used in the paper.

Part 4 also asserts Proposition 5.9: m_max^all(q) = m_max(q) - 2 for k >= 4.

Usage:  python3 Two_Periodic_Exhaustive_Verification.py [K]
        (the article uses K = 12; standard library only, Python >= 3.8)

All finite-field and polynomial computations are exact.  The only floating
point is in Part 6 (evaluating the real threshold 3J(q) <= q/2 - 2) and in
the optional spot check of Part 7; neither is used in any proof.
"""
import math
import sys

# ----------------------------------------------------------------------------
# GF(2^k) arithmetic in a polynomial basis; elements are integers (bit i <-> t^i)
# ----------------------------------------------------------------------------

def gf_mul(a, b, mod, k):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if (a >> k) & 1:
            a ^= mod
    return r


def is_primitive_poly(mod, k):
    """True iff t has multiplicative order 2^k - 1 modulo mod."""
    Q = (1 << k) - 1
    x = 1
    for i in range(1, Q + 1):
        x = gf_mul(x, 2, mod, k)
        if x == 1:
            return i == Q
    return False


def first_primitive_poly(k):
    for mod in range((1 << k) + 1, 1 << (k + 1), 2):
        if is_primitive_poly(mod, k):
            return mod
    raise ValueError(k)


class Field:
    def __init__(self, k, mod=None):
        self.k = k
        self.mod = mod if mod is not None else first_primitive_poly(k)
        self.Q = (1 << k) - 1
        self.exp = [1] * self.Q            # exp[i] = t^i
        for i in range(1, self.Q):
            self.exp[i] = gf_mul(self.exp[i - 1], 2, self.mod, k)
        self.log = {v: i for i, v in enumerate(self.exp)}

    def mul(self, a, b):
        return gf_mul(a, b, self.mod, self.k)

    def inv(self, a):
        return self.exp[(-self.log[a]) % self.Q]

    def log_base(self, v, a_exp):
        """log_alpha(v) where alpha = t^{a_exp}, gcd(a_exp, Q) = 1."""
        return (self.log[v] * pow(a_exp, -1, self.Q)) % self.Q


# ----------------------------------------------------------------------------
# Projective exponents of the substitution and swap directions
# ----------------------------------------------------------------------------

def directions(m, b, D, Q):
    """V_m(b) followed by E_m(D), all reduced mod Q."""
    V = [j % Q if j % 2 else (j - b) % Q for j in range(1, m + 1)]
    E = [(i + D) % Q for i in range(1, m)]
    return V + E


def separated_prefix(b, D, Q, limit):
    """Largest M <= limit such that the 2M-1 directions of length M are distinct."""
    seen = set()
    M = 0
    for m in range(1, limit + 1):
        v = m % Q if m % 2 else (m - b) % Q
        if v in seen:
            break
        seen.add(v)
        if m >= 2:
            w = (m - 1 + D) % Q
            if w in seen:
                break
            seen.add(w)
        M = m
    return M


def info_encodable(m, b, Q):
    """Check columns m+1, m+2 independent: beta != alpha (m even), alpha*beta != 1 (m odd)."""
    return (b % Q != 1) if m % 2 == 0 else ((b + 1) % Q != 0)


def last_two_independent(N, b, Q):
    """Columns N-1, N independent: alpha != beta (N even), alpha*beta != 1 (N odd)."""
    return (b % Q != 1) if N % 2 == 0 else ((b + 1) % Q != 0)


# ----------------------------------------------------------------------------
# Part 1-3: GF(64)
# ----------------------------------------------------------------------------

def poly_mod(a, b):
    db = b.bit_length()
    while a and a.bit_length() >= db:
        a ^= b << (a.bit_length() - db)
    return a


def poly_gcd(a, b):
    while b:
        a, b = b, poly_mod(a, b)
    return a


def part_gf7():
    print("== Part 0: GF(7) example (Appendix A) ==")
    p, P = 7, 3
    W = lambda j: 1 if j % 2 else P * P % p
    syn = lambda r: (sum(x * W(j) for j, x in enumerate(r, 1)) % p,
                     sum(x * pow(P, j, p) for j, x in enumerate(r, 1)) % p)
    a = (1, 2, 4, 6, 1, 3)
    print(f"codeword {a}: syndromes {syn(a)}")
    for i in (1, 2, 3):
        r = list(a); r[i - 1], r[i] = r[i], r[i - 1]
        e, f = syn(r)
        L = e * pow(P, i, p) % p
        R = (-1) ** (i + 1) * f * (1 + P) % p
        loc = [j for j in (1, 2, 3) if e * pow(P, j, p) % p == (-1) ** (j + 1) * f * (1 + P) % p]
        print(f"swap at i={i}: r={tuple(r)}, (e,f)=({e},{f}), (L,R)=({L},{R}), locator solutions {loc}")
    e, f = syn((1, 4, 2, 6, 1, 3))
    lhs = e * pow(1 + P * P, -1, p) * P * P % p
    rhs = f * pow(1 + P, -1, p) % p
    print(f"characteristic-two formula at i=2: {lhs} vs {rhs} (fails: {lhs != rhs})")
    a1, b1 = (0, 1, 0, 0, 1, 0, 5, 3), (1, 0, 0, 1, 0, 0, 5, 3)
    r1 = (1, 0, 0, 0, 1, 0, 5, 3)
    print(f"collision: syn(a')={syn(a1)}, syn(b')={syn(b1)}, syn(r')={syn(r1)}; "
          f"order of -3 mod 7 = {next(t for t in range(1, 7) if pow(4, t, 7) == 1)}\n")


def part_gf64():
    F = Field(6, 0b1000011)          # t^6 + t + 1
    Q = F.Q
    print("== Part 1: GF(64) sharp optimum ==")
    alpha = 0b111                    # t^2 + t + 1
    a_exp = F.log[alpha]
    assert math.gcd(a_exp, Q) == 1, "alpha not primitive"
    b = 6
    beta = F.exp[(a_exp * b) % Q]
    D = F.log_base(F.mul(1 ^ alpha, F.inv(1 ^ beta)), a_exp)
    L29 = directions(29, b, D, Q)
    print(f"alpha=7 (t^2+t+1), beta={beta}, b={b}, D={D}")
    print(f"m=29: {len(L29)} directions, {len(set(L29))} distinct; "
          f"encodable (odd m): {info_encodable(29, b, Q)}")
    L30 = directions(30, b, D, Q)
    print(f"m=30 with same parameters: {len(set(L30))} distinct of {len(L30)}")

    pairs = []
    for bb in range(2, Q):
        for DD in range(Q):
            if len(set(directions(30, bb, DD, Q))) == 59:
                pairs.append((bb, DD))
    phi63 = sum(1 << e for e in (36, 33, 27, 24, 18, 12, 9, 3, 0))
    gcds = [poly_gcd(phi63, (1 << (DD + bb)) | (1 << DD) | 0b11) for bb, DD in pairs]
    print(f"residue pairs with 59 distinct exponents at m=30: {len(pairs)}")
    print(f"  {pairs}")
    print(f"  all gcd(Phi_63, F_(b,D)) == 1: {all(g == 1 for g in gcds)}")
    # Phi_63 vanishes at all 36 primitive elements
    def ev(poly, x):
        r, p = 0, 1
        for e in range(poly.bit_length()):
            if (poly >> e) & 1:
                r ^= p
            p = F.mul(p, x)
        return r
    roots = sum(1 for i in range(Q) if math.gcd(i, Q) == 1 and ev(phi63, F.exp[i]) == 0)
    print(f"  Phi_63 vanishes at {roots}/36 primitive elements")
    # all-positions model: N = 29, columns 28, 29 independent
    print(f"all-positions model, N=29: directions distinct, last two columns "
          f"independent: {last_two_independent(29, b, Q)}  -> 27 information symbols")

    print("\n== Part 2: fixed parameters of Beaugris (2025) ==")
    a0 = 1                            # alpha = t
    b0 = 2
    beta0 = F.exp[2]
    D0 = F.log_base(F.mul(1 ^ 2, F.inv(1 ^ beta0)), a0)
    M0 = separated_prefix(b0, D0, Q, Q // 2 + 1)
    print(f"alpha=t, beta=alpha^2: D={D0}; largest separated length M={M0}")
    print(f"  V_6 = {sorted(directions(6, b0, D0, Q)[:6])}, E_6 = {directions(6, b0, D0, Q)[6:]}")
    best_fixed = (0, [])
    for bb in range(1, Q):
        Db = F.log_base(F.mul(1 ^ 2, F.inv(1 ^ F.exp[bb])), a0)
        Mb = separated_prefix(bb, Db, Q, (Q + 1) // 2)
        mb = next((mm for mm in range(Mb, 1, -1) if info_encodable(mm, bb, Q)), 0)
        if mb > best_fixed[0]:
            best_fixed = (mb, [bb])
        elif mb == best_fixed[0]:
            best_fixed[1].append(bb)
    print(f"alpha=t kept, beta=alpha^b optimized: max m = {best_fixed[0]}, attained for b in {best_fixed[1]}")

    print("\n== Part 3: explicit encoding, decoding and DNA representation ==")
    m = 29
    u = lambda j: 1 if j % 2 else beta
    apow = lambda j: F.exp[(a_exp * j) % Q]
    info = list(range(1, 30))
    S0 = S1 = 0
    for j, cj in enumerate(info, 1):
        S0 ^= F.mul(u(j), cj)
        S1 ^= F.mul(apow(j), cj)
    sols = [(x, y) for x in range(64) for y in range(64)
            if F.mul(u(30), x) ^ F.mul(u(31), y) == S0
            and F.mul(apow(30), x) ^ F.mul(apow(31), y) == S1]
    assert len(sols) == 1
    c30, c31 = sols[0]
    code = info + [c30, c31]
    s0 = s1 = 0
    for j, cj in enumerate(code, 1):
        s0 ^= F.mul(u(j), cj)
        s1 ^= F.mul(apow(j), cj)
    print(f"S0={S0}, S1={S1}, checks c30={c30}, c31={c31}, full syndrome=({s0},{s1})")

    def decode(r):
        s0 = s1 = 0
        for j, rj in enumerate(r, 1):
            s0 ^= F.mul(u(j), rj)
            s1 ^= F.mul(apow(j), rj)
        if (s0, s1) == (0, 0):
            return ("none", r, (s0, s1))
        if s0 == 0 or s1 == 0:                      # step (2): reject
            return ("detected", None, (s0, s1))
        ell = F.log_base(F.mul(s1, F.inv(s0)), a_exp)
        table = {}
        for j in range(1, m + 1):
            table[j % Q if j % 2 else (j - b) % Q] = ("sub", j)
        for i in range(1, m):
            table[(i + D) % Q] = ("swap", i)
        if ell not in table:                        # step (3): reject
            return ("detected", None, (s0, s1))
        kind, pos = table[ell]
        out = list(r)
        if kind == "sub":
            out[pos - 1] ^= F.mul(s0, F.inv(u(pos)))
        else:
            out[pos - 1], out[pos] = out[pos], out[pos - 1]
        t0 = t1 = 0                                 # step (4): re-check
        for j, oj in enumerate(out, 1):
            t0 ^= F.mul(u(j), oj)
            t1 ^= F.mul(apow(j), oj)
        if (t0, t1) != (0, 0):
            return ("detected", None, (s0, s1))
        return ((kind, pos, ell), out, (s0, s1))

    r1 = list(code); r1[4] ^= 17
    t1, o1, sy1 = decode(r1)
    r2 = list(code); r2[27], r2[28] = r2[28], r2[27]
    t2, o2, sy2 = decode(r2)
    print(f"substitution (magnitude 17, position 5): syndrome {sy1}, decoded {t1}, ok={o1 == code}")
    print(f"swap (positions 28,29): syndrome {sy2}, decoded {t2}, ok={o2 == code}")
    # decoder never outputs a non-codeword: run it on many arbitrary words
    import random
    rng = random.Random(7)
    outs = {"none": 0, "detected": 0, "corrected": 0}
    for _ in range(20000):
        w = [rng.randrange(64) for _ in range(31)]
        t, o, _ = decode(w)
        outs["corrected" if isinstance(t, tuple) else t] += 1
    print(f"decoder on 20000 random words: {outs} (every output passed the re-check)")

    nt = "acgt"
    codon = lambda v: "".join(nt[(v >> s) & 3] for s in (4, 2, 0))
    print("information codons :", " ".join(codon(v) for v in info))
    print("check codons       :", codon(c30), codon(c31))
    print(f"strand length      : {3 * len(code)} nucleotides")
    print(f"codon 5 after the substitution: {codon(5)} -> {codon(5 ^ 17)}")
    print(f"codons 28, 29: {codon(28)} {codon(29)}")


# ----------------------------------------------------------------------------
# Part 4: exhaustive optima
# ----------------------------------------------------------------------------

def frobenius_reps(Q, k):
    seen, reps = set(), []
    for a in range(1, Q):
        if math.gcd(a, Q) != 1 or a in seen:
            continue
        x = a
        for _ in range(k):
            seen.add(x)
            x = (2 * x) % Q
        reps.append(a)
    return reps


def optima(k):
    """Return (m_max, m_max_all) over all primitive alpha and beta != 1."""
    F = Field(k)
    Q = F.Q
    best_info, best_all = 0, 0
    for a in frobenius_reps(Q, k):      # Frobenius maps C_{a,b} to C_{a^2,b^2}
        alpha = F.exp[a]
        l1a = F.log_base(1 ^ alpha, a)
        for b in range(1, Q):
            beta = F.exp[(a * b) % Q]
            D = (l1a - F.log_base(1 ^ beta, a)) % Q
            M = separated_prefix(b, D, Q, (Q + 1) // 2)   # m <= q/2
            for mm in range(M, 1, -1):
                if info_encodable(mm, b, Q):
                    best_info = max(best_info, mm)
                    break
            for N in range(M, 3, -1):
                if last_two_independent(N, b, Q):
                    best_all = max(best_all, N - 2)
                    break
    return best_info, best_all


# ----------------------------------------------------------------------------
# Part 5: brute-force syndrome injectivity
# ----------------------------------------------------------------------------

def brute_check(k):
    """Compare the direction criterion with syndrome injectivity over all error vectors."""
    F = Field(k)
    Q = F.Q
    mismatches, tested = 0, 0
    for a in range(1, Q):
        if math.gcd(a, Q) != 1:
            continue
        ap = lambda j: F.exp[(a * j) % Q]
        for b in range(1, Q):
            beta = ap(b)
            alpha = ap(1)
            D = (F.log_base(1 ^ alpha, a) - F.log_base(1 ^ beta, a)) % Q
            for N in range(2, (Q + 1) // 2 + 1):
                u = lambda j: 1 if j % 2 else beta
                H = [(u(j), ap(j)) for j in range(1, N + 1)]
                errs = H + [(H[i][0] ^ H[i + 1][0], H[i][1] ^ H[i + 1][1]) for i in range(N - 1)]
                syn, inj = set(), True
                for (x, y) in errs:
                    for d in range(1, Q + 1):
                        s = (F.mul(d, x), F.mul(d, y))
                        if s == (0, 0) or s in syn:
                            inj = False
                        syn.add(s)
                crit = len(set(directions(N, b, D, Q))) == 2 * N - 1
                tested += 1
                mismatches += (crit != inj)
    return tested, mismatches


# ----------------------------------------------------------------------------
# Part 6: threshold of the near-half lower bound
# ----------------------------------------------------------------------------

def threshold():
    out = []
    for k in range(26, 33):
        q = 2 ** k
        J = math.ceil(q ** 0.75 * (2 + math.log(q - 1))) + 1
        out.append((k, 3 * J <= q // 2 - 2))
    return out


def jacobi_expansion_check(k, trials=3, seed=1):
    """Compare T = #{(x,y) in A x B : x + y = 1} with
       Q^-2 * sum_{chi,psi} Ahat(chi) Bhat(psi) Jac(chi,psi),
       Ahat(chi) = sum_{a in A} conj(chi(a)).  Floating point is used only here."""
    import cmath
    import random
    rng = random.Random(seed)
    F = Field(k)
    Q = F.Q
    chi = lambda j, x: 0 if x == 0 else cmath.exp(2j * math.pi * j * F.log[x] / Q)
    jac = {(j, l): sum(chi(j, x) * chi(l, 1 ^ x) for x in range(2, Q + 1))
           for j in range(Q) for l in range(Q)}
    worst = 0.0
    for _ in range(trials):
        A = rng.sample(range(1, Q + 1), rng.randint(3, Q // 2))
        B = rng.sample(range(1, Q + 1), rng.randint(3, Q // 2))
        T = sum(1 for x in A for y in B if x ^ y == 1)
        hatA = [sum(chi(j, a).conjugate() for a in A) for j in range(Q)]
        hatB = [sum(chi(l, y).conjugate() for y in B) for l in range(Q)]
        val = sum(hatA[j] * hatB[l] * jac[j, l] for j in range(Q) for l in range(Q)) / Q ** 2
        worst = max(worst, abs(val - T))
    return worst


if __name__ == "__main__":
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    part_gf7()
    part_gf64()
    print("\n== Part 4: exhaustive optima (Tables 1 and 2) ==")
    print(" k     q   m_max   m_max_all")
    for k in range(2, K + 1):
        mi, ma = optima(k)
        if k >= 4:
            assert ma == mi - 2, "Proposition 5.9 violated"
        print(f"{k:2d} {2**k:5d} {mi:7d} {ma:11d}", flush=True)
    print("Proposition 5.9 (m_max^all = m_max - 2 for k >= 4): confirmed")
    print("\n== Part 5: criterion vs brute-force syndrome injectivity ==")
    for k in (3, 4, 5):
        t, mm = brute_check(k)
        print(f"k={k}: {t} codes tested, {mm} mismatches")
    print("\n== Part 6: condition 3J(q) <= q/2 - 2 ==")
    print(threshold())
    print("\n== Part 7: Jacobi-sum expansion (30) ==")
    for k in (4, 5, 6):
        print(f"k={k}: max |formula - T| over random A, B = {jacobi_expansion_check(k):.2e}")
