#!/usr/bin/env python3
"""
Cantor's algorithm for genus-2 hyperelliptic curve in imaginary form:
    y^2 + h(x) y = f(x), deg f = 5, deg h <= 3
Polynomials represented as lists [a_0, a_1, ...] (low-to-high) mod p.

Approach: We transform the original curve y^2 = F(x) (deg F = 6) by
substituting y = z + c x^3 where c = sqrt(a_6), getting
z^2 + 2c x^3 z = F(x) - a_6 x^6 (degree 5 in x).
"""

# Global modulus, set by user
P = None

def setp(p): 
    global P
    P = p

def poly_norm(a):
    a = list(a)
    while len(a) > 1 and a[-1] % P == 0:
        a.pop()
    return [x % P for x in a]

def poly_deg(a):
    a = poly_norm(a)
    if len(a) == 1 and a[0] == 0:
        return -1
    return len(a) - 1

def poly_add(a, b):
    n = max(len(a), len(b))
    r = [0]*n
    for i in range(n):
        if i < len(a): r[i] = (r[i] + a[i]) % P
        if i < len(b): r[i] = (r[i] + b[i]) % P
    return poly_norm(r)

def poly_sub(a, b):
    n = max(len(a), len(b))
    r = [0]*n
    for i in range(n):
        if i < len(a): r[i] = (r[i] + a[i]) % P
        if i < len(b): r[i] = (r[i] - b[i]) % P
    return poly_norm(r)

def poly_mul(a, b):
    if not a or not b: return [0]
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b):
            r[i+j] = (r[i+j] + x*y) % P
    return poly_norm(r)

def poly_scalar(a, c):
    return poly_norm([(x*c) % P for x in a])

def poly_neg(a):
    return [(-x) % P for x in a]

def poly_divmod(a, b):
    """Return (q, r) such that a = q*b + r and deg(r) < deg(b)."""
    a = poly_norm(a)
    b = poly_norm(b)
    if poly_deg(b) < 0:
        raise ZeroDivisionError("divide by zero polynomial")
    lb = b[-1]
    lb_inv = pow(lb, -1, P)
    q = [0] * max(1, len(a) - len(b) + 1)
    r = list(a)
    while poly_deg(r) >= poly_deg(b):
        d = poly_deg(r) - poly_deg(b)
        coef = (r[-1] * lb_inv) % P
        if d < len(q):
            q[d] = coef
        # subtract coef * x^d * b from r
        for j in range(len(b)):
            r[d+j] = (r[d+j] - coef * b[j]) % P
        r = poly_norm(r)
        if poly_deg(r) < 0:
            break
    return poly_norm(q), poly_norm(r)

def poly_mod(a, b):
    return poly_divmod(a, b)[1]

def poly_make_monic(a):
    """Returns (monic_a, leading_coefficient)."""
    a = poly_norm(a)
    if poly_deg(a) < 0:
        return a, 1
    lc = a[-1]
    inv = pow(lc, -1, P)
    return poly_norm([(x*inv) % P for x in a]), lc

def poly_xgcd(a, b):
    """Return (g, s, t) with g = s*a + t*b, g monic."""
    a, b = poly_norm(a), poly_norm(b)
    old_r, r = a, b
    old_s, s = [1], [0]
    old_t, t = [0], [1]
    while poly_deg(r) >= 0:
        q, _ = poly_divmod(old_r, r)
        old_r, r = r, poly_sub(old_r, poly_mul(q, r))
        old_s, s = s, poly_sub(old_s, poly_mul(q, s))
        old_t, t = t, poly_sub(old_t, poly_mul(q, t))
    g = old_r
    # Make g monic
    g_monic, lc = poly_make_monic(g)
    if poly_deg(g) >= 0:
        inv = pow(lc, -1, P)
        old_s = poly_scalar(old_s, inv)
        old_t = poly_scalar(old_t, inv)
    return g_monic, old_s, old_t

def poly_gcd(a, b):
    g, _, _ = poly_xgcd(a, b)
    return g

# ------------- Cantor's algorithm for y^2 + h y = f ----------------

class HEC:
    def __init__(self, f, h, g=2):
        """g is genus (2 for us)."""
        self.f = poly_norm(f)
        self.h = poly_norm(h)
        self.g = g
        self.IDENTITY = ([1], [0])  # u = 1, v = 0

    def is_identity(self, D):
        u, v = D
        return poly_deg(u) == 0 and u[0] % P == 1 and poly_deg(v) < 0

    def add(self, D1, D2):
        u1, v1 = D1
        u2, v2 = D2
        # Step 1: d1 = gcd(u1, u2), u1*e1 + u2*e2 = d1
        d1, e1, e2 = poly_xgcd(u1, u2)
        # Step 2: d = gcd(d1, v1+v2+h), d1*c3 + (v1+v2+h)*c4 = d
        s_in = poly_add(poly_add(v1, v2), self.h)
        d, c3, c4 = poly_xgcd(d1, s_in)
        # s1 = c3*e1, s2 = c3*e2, s3 = c4
        s1 = poly_mul(c3, e1)
        s2 = poly_mul(c3, e2)
        s3 = c4
        # Step 3: u = u1*u2 / d^2
        u_num = poly_mul(u1, u2)
        d2 = poly_mul(d, d)
        u, rem = poly_divmod(u_num, d2)
        if poly_deg(rem) >= 0:
            raise ValueError("Composition: u not divisible by d^2")
        # Step 4: v = (s1*u1*v2 + s2*u2*v1 + s3*(v1*v2 + f - h*v1)) / d, then mod u
        # Actually: v = [s1 u1 v2 + s2 u2 v1 + s3 (v1 v2 + f) - s3 h v1] / d
        # Wait, derivation: v_new must satisfy v_new ≡ v1 mod (u1/d), v_new ≡ v2 mod (u2/d), 
        # and v_new^2 + h*v_new ≡ f mod u. The formula is:
        #   v_new = (s1*u1*v2 + s2*u2*v1 + s3*(v1*v2 + f)) / d mod u
        # for h=0. For general h, replace f with f - h*v_? Let me think...
        # From Cantor's original paper, for y^2 + h y = f:
        # The formula uses: v = (s1*u1*v2 + s2*u2*v1 + s3*(v1*v2 + f)) / d mod u
        # but the GCD step uses v1+v2+h.
        v_num = poly_add(
            poly_add(poly_mul(s1, poly_mul(u1, v2)),
                     poly_mul(s2, poly_mul(u2, v1))),
            poly_mul(s3, poly_add(poly_mul(v1, v2), self.f))
        )
        v_quot, rem = poly_divmod(v_num, d)
        if poly_deg(rem) >= 0:
            raise ValueError(f"Composition: v_num not divisible by d. rem={rem}")
        v = poly_mod(v_quot, u)
        # Make u monic (should already be if u1, u2 are)
        u, _ = poly_make_monic(u)
        # Reduce
        return self.reduce((u, v))

    def reduce(self, D):
        u, v = D
        # iterate while deg(u) > g
        while poly_deg(u) > self.g:
            # u' = (f - h*v - v^2) / u
            num = poly_sub(poly_sub(self.f, poly_mul(self.h, v)), poly_mul(v, v))
            u_new, rem = poly_divmod(num, u)
            if poly_deg(rem) >= 0:
                raise ValueError(f"Reduction: (f - h v - v^2) not divisible by u. rem deg = {poly_deg(rem)}")
            # Make u_new monic
            u_new, _ = poly_make_monic(u_new)
            # v' = (-h - v) mod u_new
            v_new = poly_mod(poly_sub(poly_neg(self.h), v), u_new)
            u, v = u_new, v_new
        return (u, v)

    def double(self, D):
        return self.add(D, D)

    def neg(self, D):
        u, v = D
        # negation: (u, -h - v mod u)
        v_new = poly_mod(poly_sub(poly_neg(self.h), v), u)
        return (u, v_new)

    def scalar_mul(self, k, D):
        if k == 0:
            return self.IDENTITY
        if k < 0:
            return self.scalar_mul(-k, self.neg(D))
        result = self.IDENTITY
        base = D
        while k > 0:
            if k & 1:
                result = self.add(result, base)
            base = self.double(base)
            k >>= 1
        return result
