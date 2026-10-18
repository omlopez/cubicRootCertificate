"""Certify root separation for the cubic polynomial

G(t; x,y,z) = 2i*t**3 + (x-i*y+2)*t**2 - 2*z*t - (x+i*y).

At (x,y,z) = forward(u,v,beta), zeta = u+i*v is a prescribed root.
On Q(lambda,beta) = 900*(lambda-3/20)**2+100*(beta-1/50)**2 <= 1,
with lambda = u*u+v*v, certify |G'(zeta)|**2 >= 1/100 and
Q(lambda_eta,beta_eta) >= 2 for both roots of the remaining quadratic
factor. The recovered parameters are defined in recovered_q().

Run: python3 certify.py
Every interval endpoint is an integer in units of 2**(-80).
The subdivision uses boxes in (u,v,beta).
"""

from math import isqrt

S = 2**80


# Outward-rounded real interval arithmetic. Intervals are pairs (lo, hi).
def ceildiv(n, d):
    return -((-n) // d)


def rational(n, d=1):
    return n*S // d, ceildiv(n*S, d)


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def neg(a):
    return -a[1], -a[0]


def sub(a, b):
    return add(a, neg(b))


def mul(a, b):
    p = (a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1])
    return min(p)//S, ceildiv(max(p), S)


def scale(a, n, d=1):
    p = (n*a[0], n*a[1])
    return min(p)//d, ceildiv(max(p), d)


def div(a, b):
    if b[0] <= 0 <= b[1]:
        raise ArithmeticError("division interval contains zero")
    return mul(a, (S*S//b[1], ceildiv(S*S, b[0])))


def square(a):
    p, q = a[0]**2, a[1]**2
    lo = 0 if a[0] <= 0 <= a[1] else min(p, q)
    return lo//S, ceildiv(max(p, q), S)


def sqrt(a):
    # Used only for expressions known mathematically to be nonnegative.
    if a[1] < 0:
        raise ArithmeticError("negative square-root interval")
    lo, hi = isqrt(max(0, a[0])*S), isqrt(a[1]*S)
    return lo, hi + (hi*hi < a[1]*S)


ZERO, ONE = rational(0), rational(1)
LAMBDA0, BETA0 = rational(3, 20), rational(1, 50)


# Complex rectangles are pairs of real intervals.
def cadd(a, b):
    return add(a[0], b[0]), add(a[1], b[1])


def cscale(a, n, d=1):
    return scale(a[0], n, d), scale(a[1], n, d)


def cmul(a, b):
    return (sub(mul(a[0], b[0]), mul(a[1], b[1])),
            add(mul(a[0], b[1]), mul(a[1], b[0])))


def norm2(a):
    return add(square(a[0]), square(a[1]))


def csqrt(a):
    # The two signs of the returned rectangle enclose both square roots.
    modulus = sqrt(norm2(a))
    imaginary = sqrt(scale(sub(modulus, a[0]), 1, 2))
    if imaginary[0] > 0:
        return div(a[1], scale(imaginary, 2)), imaginary
    real = sqrt(scale(add(modulus, a[0]), 1, 2))
    if real[0] > 0:
        return real, div(a[1], scale(real, 2))
    bound = sqrt(modulus)[1]
    return (-bound, bound), (-bound, bound)


def q(lam, beta):
    # Q(lambda, beta) = 900*(lambda-3/20)^2 + 100*(beta-1/50)^2.
    return add(scale(square(sub(lam, LAMBDA0)), 900),
               scale(square(sub(beta, BETA0)), 100))


def forward(u, v, beta):
    # Rational map X(u,v,beta) at which G(u+i*v; X) = 0.
    lam = add(square(u), square(v))
    d = add(ONE, lam)
    gamma = div(sub(lam, ONE), d)
    kappa_coef = div(scale(lam, 4), square(d))
    uv = mul(u, v)
    k = add(sub(beta, mul(add(gamma, kappa_coef), u)),
            mul(add(scale(gamma, 2), kappa_coef), uv))
    ell = add(v, sub(square(u), square(v)))
    x = div(scale(sub(mul(u, k), mul(v, ell)), 2), d)
    y = div(scale(add(mul(v, k), mul(u, ell)), 2), d)
    z = add(mul(gamma, beta),
            mul(kappa_coef, sub(mul(sub(ONE, gamma), u),
                                mul(sub(rational(2), gamma), uv))))
    return x, y, z


def recovered_q(root, point):
    # Q at the recovered parameters (lambda_eta, beta_eta) of a root eta.
    u, v = root
    x, y, z = point
    lam = norm2(root)
    d = add(ONE, lam)
    numerator = add(scale(add(mul(u, x), mul(v, y)), 2),
                    mul(sub(lam, ONE), z))
    correction = div(scale(mul(mul(lam, u), sub(ONE, v)), 4), square(d))
    beta = add(div(numerator, d), correction)
    return q(lam, beta)


def bounds(u, v, beta):
    point = forward(u, v, beta)
    x, y, z = point
    root = u, v
    # D = 6i*zeta^2 + 2*(conj(w)+2)*zeta - 2*z.
    derivative = cadd(
        cmul((ZERO, rational(6)), cmul(root, root)),
        cadd(cmul((add(x, rational(2)), neg(y)), cscale(root, 2)),
             (scale(z, -2), ZERO)))
    d_lower = norm2(derivative)[0]
    # G(Z)/(2i) = (Z-zeta)*(Z^2 + L*Z + C), exactly.
    L = cadd(root, (scale(y, -1, 2),
                    scale(add(x, rational(2)), -1, 2)))
    C = cadd(cmul(root, L), (ZERO, z))
    discriminant = cadd(cmul(L, L), cscale(C, -4))
    radical = csqrt(discriminant)
    rival_roots = [
        cscale(cadd(cscale(L, -1), cscale(radical, sign)), 1, 2)
        for sign in (-1, 1)]
    q_lower = min(recovered_q(eta, point)[0] for eta in rival_roots)
    return d_lower, q_lower


def certify():
    # Q <= 1 implies |u|,|v| < 1/2 and -2/25 <= beta <= 3/25.
    # The slightly outward-rounded rectangle contains the entire closed set.
    box = ((-S//2, S//2), (-S//2, S//2),
           (rational(-2, 25)[0], rational(3, 25)[1]))
    stack = [(box, 0)]
    certified = outside = deepest = 0
    while stack:
        box, depth = stack.pop()
        u, v, beta = box
        if q(add(square(u), square(v)), beta)[0] > S:
            outside += 1
            continue
        d_lower, q_lower = bounds(u, v, beta)
        if 100*d_lower >= S and q_lower >= 2*S:
            certified += 1
            deepest = max(deepest, depth)
            continue
        if depth >= 60:
            raise ArithmeticError("unresolved box: " + str(box))
        axis = depth % 3
        lo, hi = box[axis]
        mid = (lo+hi)//2
        left, right = list(box), list(box)
        left[axis], right[axis] = (lo, mid), (mid, hi)
        # Both resulting closed boxes are retained: coverage follows by bisection.
        stack.extend(((tuple(right), depth+1), (tuple(left), depth+1)))
    print("PASS: root-separation bounds certified.")
    print("Certified:", certified, "Outside:", outside,
          "Maximum certified depth:", deepest)


if __name__ == "__main__":
    certify()
