#!/usr/bin/env python3
"""
Stare at the prime. Look for structure (smoothness, special form, etc.).

This is the first reconnaissance step described in README §2.
"""
from sympy import factorint, isprime

p = 129403459552990578380563458675806698255602319995627987262273876063027199999999

print(f"p              = {p}")
print(f"p.bit_length() = {p.bit_length()}")
print(f"is prime?      = {isprime(p)}")
print(f"p mod 4        = {p % 4}      # 3 means sqrt is just a^((p+1)/4)")
print()

print("Factor p + 1 (smoothness check):")
fp1 = factorint(p + 1)
for q, e in sorted(fp1.items()):
    print(f"  {q}^{e}")
print()

print("Factor p - 1 (for comparison):")
fp_1 = factorint(p - 1, limit=10**6)  # partial; bigger factor left for clarity
for q, e in sorted(fp_1.items()):
    print(f"  {q}^{e}")
print()

prod = 1
for q, e in fp1.items():
    prod *= q ** e
print(f"All of p+1 explained by small primes? {prod == p + 1}")
print(f"Largest prime dividing p+1: {max(fp1)}")
