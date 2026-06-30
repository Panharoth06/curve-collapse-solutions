#!/usr/bin/env python3
"""
Factor the curve polynomial f(x) over F_p.

If f(x) splits into low-degree pieces, the Jacobian has lots of 2-torsion
and is a strong candidate for a Pohlig-Hellman attack. See README §2.
"""
from sympy import GF, Poly, symbols

p = 129403459552990578380563458675806698255602319995627987262273876063027199999999

f_coeffs = [
    87455262955769204408909693706467098277950190590892613056321965035180446006909,
    12974562908961912291194866717212639606874236186841895510497190838007409517645,
    11783716142539985302405554361639449205645147839326353007313482278494373873961,
    55538572054380843320095276970494894739360361643073391911629387500799664701622,
    124693689608554093001160935345506274464356592648782752624438608741195842443294,
    52421364818382902628746436339763596377408277031987489475057857088827865195813,
    50724784947260982182351215897978953782056750224573008740629192419901238915128,
]

x = symbols("x")
f_expr = sum(c * x**i for i, c in enumerate(f_coeffs))
f_poly = Poly(f_expr, x, domain=GF(p))

leading, factors = f_poly.factor_list()
print(f"f(x) has degree {f_poly.degree()} over F_p")
print(f"Leading coefficient: {leading}")
print(f"Number of irreducible factors: {len(factors)}")
print()
print("Factors (each one over F_p):")
for fac, mult in factors:
    coeffs = [int(c) % p for c in fac.all_coeffs()]
    print(f"  degree {fac.degree()} (multiplicity {mult}): {coeffs}")
