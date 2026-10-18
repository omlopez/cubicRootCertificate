# Certified root separation for a cubic polynomial

Interval-arithmetic certificate of a root-separation property of

$$
G(t;x,y,z)=2it^3+(x-iy+2)t^2-2zt-(x+iy).
$$

The statement and definitions are in the docstring of [`certify.py`](certify.py).

```sh
python3 certify.py
```

Expected output:

```
PASS: root-separation bounds certified.
Certified: 343037 Outside: 21733 Maximum certified depth: 28
```

The initial box covers the stated domain, and each bisection retains both
closed halves, preserving coverage. A box is excluded only when its interval
lower bound for Q exceeds 1. Every other terminal box is certified using
outward-rounded interval arithmetic to satisfy |G'(zeta)|^2 >= 1/100 and
Q(lambda_eta, beta_eta) >= 2 for both rival roots.

`PASS` is printed only after all boxes have been processed. An unresolved box
at depth 60 raises an exception before that message; this means the run did
not certify the claim, not necessarily that the claim is false.
