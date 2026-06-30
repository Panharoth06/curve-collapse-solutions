# Curve Collapse — A Genus-2 Story

> *"My friend invented 'military-grade' curve crypto after 3 coffees and 1 wiki page.
> Break it, decrypt it, and humble him."* — Cyrus

**Flag:** `MPTC{my_fri3nd_15_n0t_4_g3n1u5}`

---

## Quick Start

Just want the flag? Clone the repo and run:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python solve.py
# → [+] Flag: MPTC{my_fri3nd_15_n0t_4_g3n1u5}
```

The whole attack finishes in about 4 seconds. Read on for the full story of how it works.

---

## Table of Contents

- [Quick Start](#quick-start)
- [0. The Cast](#0-the-cast)
- [1. What on Earth Is a Hyperelliptic Curve?](#1-what-on-earth-is-a-hyperelliptic-curve)
- [2. First Move: Look at the Numbers](#2-first-move-look-at-the-numbers)
- [3. The Plan](#3-the-plan)
- [4. Reshaping the Curve](#4-reshaping-the-curve)
- [5. Cantor's Algorithm in Two Paragraphs](#5-cantors-algorithm-in-two-paragraphs)
- [6. Confirming the Group Order](#6-confirming-the-group-order)
- [7. The Pothole: Two Representations of the Same Point](#7-the-pothole-two-representations-of-the-same-point)
- [8. Pohlig–Hellman in One Picture](#8-pohlighellman-in-one-picture)
- [9. Cracking the Flag](#9-cracking-the-flag)
- [10. Running the Proof-of-Concept](#10-running-the-proof-of-concept)
- [11. Lessons Learned](#11-lessons-learned)

---

## 0. The Cast

Before we start, meet everyone in our story.

| Item | What it is |
|---|---|
| `p` | A 257-bit prime |
| `f(x)` | A degree-6 polynomial. Defines the curve `y² = f(x)` |
| `G = (G_u, G_v)` | A "generator" divisor on the Jacobian (Mumford representation) |
| `Q = (Q_u, Q_v)` | Public point such that `Q = k·G` for a secret integer `k` |
| `enc_flag` | The flag, encrypted using a key derived from `k` |

The cryptosystem is the **Hyperelliptic Curve Discrete Log Problem (HECDLP)**: given `G` and `Q = k·G`, recover `k`. If `k` is hard to find, the system is secure. Our job is to show it is not.

The values themselves live in `val.txt`. We treat them as the inputs to our adventure.

---

## 1. What on Earth Is a Hyperelliptic Curve?

For elliptic curves (the cousins everyone knows), points satisfy `y² = x³ + ax + b`. You can "add" two points to get a third, and "multiplying" a point by an integer is just adding it to itself.

A **genus-2 hyperelliptic curve** generalises this. The equation becomes

```
y² = f(x)
```

where `f(x)` has degree `2g+1` or `2g+2`. Here `g = 2`, so `deg f ∈ {5, 6}`. Our `f` has degree 6.

Two new complications appear:

1. **Single points don't form a group on their own.** Instead, you work with *formal sums of points* called **divisors**. A divisor of weight ≤ 2 (i.e. up to two points "added together") is the typical citizen on the Jacobian.

2. **Mumford representation.** Each divisor of weight ≤ 2 is stored as a pair of polynomials `(u(x), v(x))`:
   - `u(x)` is monic, `deg u ≤ 2`. Its roots tell you the x-coordinates of the points.
   - `v(x)` has `deg v < deg u`. Plugging an x-root of `u` into `v` gives the y-coordinate.
   - A consistency condition holds: `v² ≡ f (mod u)`.

So when the challenge says

```python
G_u = [95640..., 1]   # u(x) = x + 95640...
G_v = [23400...]      # v(x) = 23400... (a constant)
```

it means `G` is one point: `x = -95640... mod p`, `y = 23400... mod p`.

For `Q`:

```python
Q_u = [342770..., 343503..., 1]   # u(x) = x² + 343503·x + 342770
Q_v = [102912..., 65726...]       # v(x) = 65726·x + 102912
```

`Q` represents two points (the x-coordinates are roots of that quadratic).

The "group operation" on these `(u, v)` pairs is **Cantor's algorithm**. Think of it as a recipe for the `+` button on a hyperelliptic-curve calculator.

---

## 2. First Move: Look at the Numbers

Before bashing on crypto with brute force, look at the inputs. Mathematicians call this *staring at the problem*.

```python
# check.py — what does the prime look like?
p = 129403459552990578380563458675806698255602319995627987262273876063027199999999
print(p.bit_length())            # 257
print(p % 4)                     # 3  → handy for square roots
print(factorint(p + 1))          # smoothness check
print(factorint(p - 1))
```

Output:

```
p + 1 = 2^23 · 3^14 · 5^8 · 7^4 · 11^10 · 13^10 · 17^9 · 19^6 · 23^5 · 29 · 31^4
p - 1 = 2 · 113177 · 34591397 · (one huge prime)
```

**Whoa**. `p + 1` is *insanely* smooth — every prime factor is at most 31. `p - 1`, in contrast, contains a giant prime. That asymmetry is a big hint: somebody designed `p` specifically so that `p + 1` is super smooth.

In curve crypto land, the order of the group of points is usually close to `p + 1` (or `(p + 1)²` for genus 2). When that number factors into tiny primes, an old attack named **Pohlig–Hellman** crushes the discrete log. The fact that `p + 1` is fully smooth is therefore screaming "I am breakable."

Next, let's check the polynomial `f(x)`.

```python
# analyze.py — does f(x) split nicely?
f = Poly(sum(c*x**i for i, c in enumerate(f_coeffs)), x, domain=GF(p))
print(f.factor_list())
```

We discover that `f(x)` factors into **three irreducible quadratics** over `F_p`:

```
f(x) = a₆ · q₁(x) · q₂(x) · q₃(x)
```

A "random" degree-6 polynomial almost never factors this neatly. Like the smooth `p + 1`, this is intentional — and it forces the Jacobian to have plenty of 2-torsion (elements `D` with `2D = 0`).

So our suspicions: the Jacobian is **supersingular** with `#J(F_p) = (p+1)²`, and the Pohlig–Hellman attack will eat it for breakfast.

---

## 3. The Plan

Here is the road map:

1. **Re-shape the curve.** The "deg 6" form is awkward to compute on. We change variables to get an equivalent "deg 5" form `z² + h(x)·z = g(x)`, which Cantor's algorithm prefers.
2. **Write Cantor's algorithm.** Pure-Python polynomial arithmetic modulo `p`.
3. **Confirm the group order.** Verify that `(p+1)² · G = identity` so we know the order divides `(p+1)²`.
4. **Pohlig–Hellman.** Solve the discrete log in each tiny prime-power subgroup, then glue the answers with CRT.
5. **Decrypt the flag.** The key is derived from `k`.

Let's go.

---

## 4. Reshaping the Curve

Our curve is `y² = f(x)` with `deg f = 6`. The leading coefficient `a₆` happens to be a *quadratic residue* modulo `p` (you can check `a₆^((p-1)/2) ≡ 1`). That means there's a square root `c ∈ F_p` with `c² = a₆`.

**Trick:** substitute `y = z + c·x³`. Squaring gives

```
y² = z² + 2c·x³·z + c²·x⁶ = z² + 2c·x³·z + a₆·x⁶
```

Now `y² = f(x)` becomes

```
z² + 2c·x³·z = f(x) − a₆·x⁶
```

The right-hand side **loses the `x⁶` term**, so it now has degree 5. The new curve has the shape

```
z² + h(x)·z = g(x),    h(x) = 2c·x³,    deg g = 5
```

which behaves more like the cosy "imaginary" model that Cantor's algorithm is built for.

We also have to translate `G` and `Q` into the new `(u, z-poly)` coordinates. The x-coordinates stay the same; we only update the `v` polynomial:

```
v_new(x) = v_old(x) − c·x³   (mod u(x))
```

Code:

```python
a6 = f_coeffs[6]
c = pow(a6, (p + 1) // 4, p)        # √a₆ mod p (p ≡ 3 mod 4)
h_poly     = [0, 0, 0, (2*c) % p]   # 2c·x³
f_new_poly = f_coeffs[:6]           # drop a₆·x⁶

c_x3 = [0, 0, 0, c]                  # c·x³
G_v_new = poly_mod(poly_sub(G_v_old, c_x3), G_u)
Q_v_new = poly_mod(poly_sub(Q_v_old, c_x3), Q_u)
```

That's our new arena.

---

## 5. Cantor's Algorithm in Two Paragraphs

Given two divisors `D₁ = (u₁, v₁)` and `D₂ = (u₂, v₂)`, the sum `D₁ + D₂` is computed in two phases:

**Compose.** Find common roots and stitch the polynomials together:

```
d₁ = gcd(u₁, u₂)            (with cofactors e₁, e₂)
d  = gcd(d₁, v₁ + v₂ + h)   (with cofactors c₃, c₄)
u  = u₁·u₂ / d²
v  = (s₁·u₁·v₂ + s₂·u₂·v₁ + s₃·(v₁·v₂ + f)) / d   (mod u)
```

The composed `u` may have degree up to 4, which is too big.

**Reduce.** While `deg u > 2`:

```
u′ = (f − h·v − v²) / u    (then make monic)
v′ = (−h − v) mod u′
```

Each reduction step shrinks `deg u` until we're back to weight ≤ 2.

For scalar multiplication `k·G`, do good old square-and-multiply: at each bit of `k`, double the current divisor, and add `G` whenever the bit is 1.

The implementation lives in `cantor.py`. The polynomial arithmetic is all hand-rolled lists like `[a₀, a₁, a₂, ...]` with operations mod `p`.

```python
from cantor import setp, HEC
setp(p)
hec = HEC(f_new_poly, h_poly, g=2)
G   = (G_u, G_v_new)
Q   = (Q_u, Q_v_new)
two_G = hec.double(G)           # 2·G
big_G = hec.scalar_mul(123, G)  # 123·G
```

---

## 6. Confirming the Group Order

We hypothesised that the Jacobian is supersingular with `#J(F_p) = (p+1)²`. Let's *test* the hypothesis instead of trusting it:

```python
N = (p + 1) ** 2
result = hec.scalar_mul(N, G)
print(hec.is_identity(result))   # → True ✓
```

Identity. So `ord(G)` divides `(p+1)²`. Even better, when we try the smaller candidate `N = p + 1`:

```python
result = hec.scalar_mul(p + 1, G)
print(hec.is_identity(result))   # → True ✓
```

`ord(G)` actually divides `p + 1`. (The group is structured like `Z/(p+1) × Z/(p+1)`, not a single cyclic `Z/(p+1)²`.) Great — the smaller, fully smooth modulus is all we need.

---

## 7. The Pothole: Two Representations of the Same Point

While testing, I hit a wild bug. The element

```
G_top := 2²² · G_pp
```

should generate a subgroup of order 2. Its sibling

```
h_1 := 2²¹ · Q_pp
```

should be either the identity or `G_top`. But Python claimed `h_1 ≠ G_top`, while `G_top − h_1 = identity` algebraically.

The reason: even after our reshape, the curve still has **two points at infinity**, so the Mumford `(u, v)` pair is not a canonical name — two different `(u, v)` pairs can describe the same Jacobian element.

The fix is simple and robust:

```python
def D_eq(D1, D2):
    return hec.is_identity(hec.add(D1, hec.neg(D2)))
```

We compare by subtracting and checking if the result is identity. A bit slower per comparison, but bullet-proof. From this moment on, all equality tests go through `D_eq`.

This is the kind of bug that turns a five-minute solve into a two-hour solve. Lesson: when the math says equal but your code says not equal, the bug is in your representation, not in the math.

---

## 8. Pohlig–Hellman in One Picture

Pohlig–Hellman is a divide-and-conquer attack on the discrete log.

If `ord(G) = ℓ₁^{e₁} · ℓ₂^{e₂} · … · ℓ_n^{e_n}` factors into small prime powers, you solve a tiny DLP for each `ℓ_i^{e_i}` and then glue the answers with the Chinese Remainder Theorem.

The recipe for one prime power `ℓ^e`:

1. **Project into the ℓ-component.**

   ```
   G_pp = (N / ℓ^e) · G        Q_pp = (N / ℓ^e) · Q
   ```

   Now `G_pp` has order dividing `ℓ^e`.

2. **Build a tiny generator of order ℓ.**

   ```
   G_top = ℓ^(e-1) · G_pp
   ```

3. **Find `k mod ℓ^e` one digit at a time.** Write `k = k₀ + k₁·ℓ + k₂·ℓ² + …`. At step `i`:

   ```
   Q_i = Q_pp − (current k partial) · G_pp
   h_i = ℓ^(e-1-i) · Q_i           # forced into the order-ℓ subgroup
   k_i = brute-force DL of h_i in base G_top   # at most ℓ tries
   ```

For our `(p+1)` the biggest prime is 31, so each "brute force" only tries 31 possibilities. The most expensive factor was `13¹⁰` ≈ 2³⁷, but Pohlig–Hellman turns it into about `10·13 ≈ 130` group operations. Trivial.

The full PHP solver is in `solve.py`. Output:

```
ell=2^23:  k ≡ 5374250         (mod 8388608)        [0.19s]
ell=3^14:  k ≡ 1886420         (mod 4782969)        [0.16s]
ell=5^8:   k ≡ 188545          (mod 390625)         [0.08s]
ell=7^4:   k ≡ 1665            (mod 2401)           [0.04s]
ell=11^10: k ≡ 4250443025      (mod 25937424601)    [0.20s]
ell=13^10: k ≡ 104022560557    (mod 137858491849)   [0.24s]
ell=17^9:  k ≡ 23450484380     (mod 118587876497)   [0.20s]
ell=19^6:  k ≡ 41475237        (mod 47045881)       [0.12s]
ell=23^5:  k ≡ 4255640         (mod 6436343)        [0.10s]
ell=29^1:  k ≡ 10              (mod 29)             [0.01s]
ell=31^4:  k ≡ 337189          (mod 923521)         [0.07s]

k ≡ 91527621348541142496688581834442276703691715094599257862319082414424378704170
   (mod p+1)
k·G == Q  → True ✓
```

Less than four seconds. The whole "military-grade curve" was a paper tiger.

---

## 9. Cracking the Flag

The encrypted flag is

```
enc_flag = f9d31f988581d7f9f06239bf26513851d32e73e7ca713aae437ce2e7419a46
```

That's 31 bytes. Now we just guess a small list of typical key derivations and XOR. The first try that yields printable ASCII wins.

```python
import hashlib
k = 91527621348541142496688581834442276703691715094599257862319082414424378704170
enc = bytes.fromhex("f9d31f988581d7f9f06239bf26513851d32e73e7ca713aae437ce2e7419a46")

for label, key in [
    ("sha256(str(k))",  hashlib.sha256(str(k).encode()).digest()),
    ("sha256(k_bytes)", hashlib.sha256(k.to_bytes(32, 'big')).digest()),
    ("sha1(str(k))",    hashlib.sha1(str(k).encode()).digest()),
    ("md5(str(k))",     hashlib.md5(str(k).encode()).digest()),
]:
    out = bytes(a ^ b for a, b in zip(enc, key))
    print(label, "→", out)
```

Output:

```
sha256(str(k))  → b'MPTC{my_fri3nd_15_n0t_4_g3n1u5}'
```

There she blows.

```
🚩 MPTC{my_fri3nd_15_n0t_4_g3n1u5}
```

---

## 10. Running the Proof-of-Concept

The repository ships with four files that do the work:

- `check.py` — peek at the prime and factor `p ± 1`.
- `analyze.py` — factor `f(x)` over `F_p`.
- `cantor.py` — pure-Python Cantor's algorithm on a genus-2 Jacobian.
- `solve.py` — end-to-end attack: reshape → Pohlig–Hellman → decrypt.

Setup and run:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python solve.py
```

Expected output (roughly four seconds on a laptop):

```
[*] Verifying group order: (p+1) * G should be identity
    ok

[*] Pohlig-Hellman on each prime-power factor of (p+1):
  2^23 : k = 5374250        (mod 8388608) [0.18s]
  3^14 : k = 1886420        (mod 4782969) [0.15s]
  5^8  : k = 188545         (mod 390625) [0.08s]
  7^4  : k = 1665           (mod 2401) [0.04s]
  11^10: k = 4250443025     (mod 25937424601) [0.20s]
  13^10: k = 104022560557   (mod 137858491849) [0.24s]
  17^9 : k = 23450484380    (mod 118587876497) [0.20s]
  19^6 : k = 41475237       (mod 47045881) [0.12s]
  23^5 : k = 4255640        (mod 6436343) [0.10s]
  29^1 : k = 10             (mod 29) [0.01s]
  31^4 : k = 337189         (mod 923521) [0.07s]

[*] Recovered k in 3.82s
    k = 91527621348541142496688581834442276703691715094599257862319082414424378704170

[*] Sanity check: k * G == Q ?
    ok

[*] Decrypting flag with sha256(str(k))...

[+] Flag: MPTC{my_fri3nd_15_n0t_4_g3n1u5}
```

---

## 11. Lessons Learned

- **Always look at the numbers.** A 30-second factor check (`factorint(p+1)`) told us the curve was broken before we wrote a single line of group code.
- **Smooth group orders = death.** When the group order is built only from small primes, Pohlig–Hellman wins.
- **Special-shaped polynomials are flags.** `f(x)` splitting into three quadratics over `F_p` was a strong tell that the Jacobian had rich 2-torsion / supersingular structure.
- **Representation ≠ identity.** Two distinct `(u, v)` pairs can describe the same Jacobian element on a "real model" curve. Compare via subtraction, not tuple equality.
- **Test as you go.** "Compute `(p+1)·G`" and "compute `2·G + (−2G)`" are cheap sanity checks that catch implementation bugs before they snowball.

The "military-grade" curve fell to a textbook attack from the 1970s. Sometimes the best crypto attack is just *reading the parameters*.
