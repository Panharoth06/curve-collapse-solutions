#!/usr/bin/env python3
"""
End-to-end exploit for the "Curve Collapse" challenge.

Attack outline (see README.md for the full story):
  1. The curve y^2 = f(x) has deg(f) = 6. Substitute y = z + c*x^3 with
     c = sqrt(a_6) to get the easier-to-handle form z^2 + h(x) z = g(x)
     with deg g = 5.
  2. Implement Cantor's algorithm for the Jacobian (cantor.py).
  3. p + 1 is fully smooth and (p+1) * G = identity, so apply
     Pohlig-Hellman across each prime-power factor of (p+1).
  4. Decrypt enc_flag with sha256(str(k)).
"""
import hashlib
import time
from cantor import HEC, poly_mod, poly_norm, poly_sub, setp
from sympy import factorint
from sympy.ntheory.modular import crt as sympy_crt

# --------------------------------------------------------------------------
# Challenge data (from val.txt)
# --------------------------------------------------------------------------
p = 129403459552990578380563458675806698255602319995627987262273876063027199999999
setp(p)

f_coeffs = [
    87455262955769204408909693706467098277950190590892613056321965035180446006909,
    12974562908961912291194866717212639606874236186841895510497190838007409517645,
    11783716142539985302405554361639449205645147839326353007313482278494373873961,
    55538572054380843320095276970494894739360361643073391911629387500799664701622,
    124693689608554093001160935345506274464356592648782752624438608741195842443294,
    52421364818382902628746436339763596377408277031987489475057857088827865195813,
    50724784947260982182351215897978953782056750224573008740629192419901238915128,
]
G_u_raw = [95640493847532285274015733349271558012724241405617918614689663966283911276425, 1]
G_v_raw = [23400917335266251424562394829509514520732985938931801439527671091919836508525]
Q_u_raw = [
    34277069903919260496311859860543966319397387795368332332841962946806971944007,
    343503204040841221074922908076232301549085995886639625441980830955087919004,
    1,
]
Q_v_raw = [
    102912018107558878490777762211244852581725648344091143891953689351031146217393,
    65726604025436600725921245450121844689064814125373504369631968173219177046384,
]
enc_flag = bytes.fromhex(
    "f9d31f988581d7f9f06239bf26513851d32e73e7ca713aae437ce2e7419a46"
)

# --------------------------------------------------------------------------
# Step 1: Reshape  y^2 = f(x)  into  z^2 + h(x) z = g(x)  with deg g = 5
# --------------------------------------------------------------------------
a6 = f_coeffs[6]
# p ≡ 3 (mod 4), so square root is just a^((p+1)/4)
c = pow(a6, (p + 1) // 4, p)
assert (c * c) % p == a6 % p, "square root failed"

h_poly = [0, 0, 0, (2 * c) % p]   # h(x) = 2c * x^3
f_new = list(f_coeffs[:6])         # drop the a_6 * x^6 term

# Move G and Q into the new coordinates: v_new = v_old - c*x^3 (mod u)
c_x3 = [0, 0, 0, c]
G_u = poly_norm(G_u_raw)
G_v = poly_mod(poly_sub(poly_norm(G_v_raw), c_x3), G_u)
Q_u = poly_norm(Q_u_raw)
Q_v = poly_mod(poly_sub(poly_norm(Q_v_raw), c_x3), Q_u)

hec = HEC(f_new, h_poly, g=2)
G = (G_u, G_v)
Q = (Q_u, Q_v)


# --------------------------------------------------------------------------
# Step 2: Robust equality. The transformed curve still has 2 infinity
# points, so two different (u, v) pairs can describe the same Jacobian
# element. Subtract and test for identity instead of comparing tuples.
# --------------------------------------------------------------------------
def D_eq(D1, D2):
    return hec.is_identity(hec.add(D1, hec.neg(D2)))


# --------------------------------------------------------------------------
# Step 3: Pohlig-Hellman
# --------------------------------------------------------------------------
def subgroup_order(D, ell, e_max):
    """Smallest k in [0, e_max] with ell^k * D = identity."""
    cur = D
    if hec.is_identity(cur):
        return 0
    for k in range(1, e_max + 1):
        cur = hec.scalar_mul(ell, cur)
        if hec.is_identity(cur):
            return k
    raise RuntimeError("order exceeds bound")


def bf_dlog(base, target, ord_bound):
    """Brute-force DL of `target` in <base> assuming order <= ord_bound."""
    cur = hec.IDENTITY
    if D_eq(cur, target):
        return 0
    for i in range(1, ord_bound):
        cur = hec.add(cur, base)
        if D_eq(cur, target):
            return i
    raise RuntimeError(f"no DL found within {ord_bound}")


def php_prime_power(Gpp, Qpp, ell, e):
    """Solve k * Gpp = Qpp where ord(Gpp) divides ell^e."""
    e_actual = subgroup_order(Gpp, ell, e)
    if e_actual == 0:
        return 0, 1
    G_top = hec.scalar_mul(ell ** (e_actual - 1), Gpp)
    k = 0
    for i in range(e_actual):
        Q_i = hec.add(Qpp, hec.neg(hec.scalar_mul(k, Gpp)))
        h_i = hec.scalar_mul(ell ** (e_actual - 1 - i), Q_i)
        k_i = bf_dlog(G_top, h_i, ell)
        k += k_i * (ell ** i)
    return k, ell ** e_actual


def pohlig_hellman(N):
    remainders, moduli = [], []
    for ell, e in sorted(factorint(N).items()):
        cofactor = N // (ell ** e)
        Gp = hec.scalar_mul(cofactor, G)
        Qp = hec.scalar_mul(cofactor, Q)
        t0 = time.time()
        k_pp, m = php_prime_power(Gp, Qp, ell, e)
        elapsed = time.time() - t0
        print(f"  {ell}^{e:<3}: k = {k_pp:<14} (mod {m}) [{elapsed:.2f}s]")
        if m > 1:
            remainders.append(k_pp)
            moduli.append(m)
    k, mod = sympy_crt(moduli, remainders)
    return int(k), int(mod)


# --------------------------------------------------------------------------
# Run the attack
# --------------------------------------------------------------------------
def main():
    print("[*] Verifying group order: (p+1) * G should be identity")
    assert hec.is_identity(hec.scalar_mul(p + 1, G))
    print("    ok\n")

    print("[*] Pohlig-Hellman on each prime-power factor of (p+1):")
    t0 = time.time()
    k, mod = pohlig_hellman(p + 1)
    print(f"\n[*] Recovered k in {time.time() - t0:.2f}s")
    print(f"    k = {k}")
    print(f"    k.bit_length() = {k.bit_length()}")

    print("\n[*] Sanity check: k * G == Q ?")
    assert D_eq(hec.scalar_mul(k, G), Q)
    print("    ok\n")

    print("[*] Decrypting flag with sha256(str(k))...")
    key = hashlib.sha256(str(k).encode()).digest()
    flag = bytes(a ^ b for a, b in zip(enc_flag, key))
    print(f"\n[+] Flag: {flag.decode()}")


if __name__ == "__main__":
    main()
