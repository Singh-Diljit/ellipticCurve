# ellipticCurve

## Overview

This project implements an **Elliptic Curve Class** with group operations for points on an elliptic curve and includes functionality to compute the order of points using **Mazur's Theorem**. Additionally, there is a **Point Class** representing points on the elliptic curve. All mathematical operations are performed using Python’s **fractions module** to ensure exact results without floating-point approximations.
## Classes

### 1. `EllipticCurve`
Represents an elliptic curve (over $\mathbb{Q}$) defined by the equation: $y^2 = x^3 + ax + b$.
### 2. `Coordinate` and `Point`
Represents either a point on an elliptic curve (`Point` class) or a a point in $\mathbb{Q}^2$ (`Coordinate` class), the point in either class may be the 'point at infinity'.
## Mazur's Theorem and Nagell-Lutz

Nagell-Lutz tells us if a point is a candidate for having finite order without having to perform any group operations. Mazur's Theorem restricts the possible orders of torsion points on elliptic curves over the rationals (it does more by telling us the group structure but that is not used). For the purpose of the `order` method the maximum finite order of a point is 12.
## Usage Example

```python
# Define the elliptic curve y^2 = x^3 + 8 over the rationals
curve = EllipticCurve(0, 8)

# Define a point on the curve
P = Point(curve, 2, 4)
Q = Point(curve, 1, 3)

print(4*Q)
>>> (31073/2704, 5491823/140608)

A = (P + 2*Q)
print(A)
>>> (2, -4)

print(A.order)
>>> inf

Z = Point(e, -2, 0)
print(Z.order)
>>> 2
```

## Dependencies

- Python 3.8+
- Uses the `fractions` module for exact arithmetic.
