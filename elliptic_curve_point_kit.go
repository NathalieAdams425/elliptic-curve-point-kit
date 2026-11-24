// Package elliptic_curve_point_kit provides affine and projective point
// arithmetic for short Weierstrass elliptic curves over the rational numbers.
//
// The curve is described by y^2 = x^3 + a*x + b. All coordinates are float64,
// which makes the package useful for teaching and for sketching, but it is not
// suitable for cryptography: floating point cannot represent field elements
// exactly.
package elliptic_curve_point_kit

import "fmt"

// Curve is a short Weierstrass curve y^2 = x^3 + a*x + b over the reals.
type Curve struct {
	A, B float64
}

// NewCurve returns a curve with the given coefficients.
func NewCurve(a, b float64) *Curve {
	return &Curve{A: a, B: b}
}

// Discriminant returns 4*a^3 + 27*b^2. A non-zero discriminant means the
// curve is non-singular (no cusps or self-intersections).
func (c *Curve) Discriminant() float64 {
	return 4*c.A*c.A*c.A + 27*c.B*c.B
}

// IsNonSingular reports whether the curve is smooth.
func (c *Curve) IsNonSingular() bool {
	return c.Discriminant() != 0
}

// AffinePoint is a point in affine coordinates. The zero value (0, 0) is not
// automatically the identity; use the Infinity flag to mark the point at
// infinity. This explicit flag avoids relying on a sentinel coordinate.
type AffinePoint struct {
	X, Y     float64
	Infinity bool
}

// AffineInfinity returns the identity element for affine arithmetic.
func AffineInfinity() AffinePoint {
	return AffinePoint{Infinity: true}
}

// IsInfinity reports whether p is the point at infinity.
func (p AffinePoint) IsInfinity() bool {
	return p.Infinity
}

// OnCurve reports whether p lies on c. The point at infinity is considered to
// lie on every curve.
func (p AffinePoint) OnCurve(c *Curve) bool {
	if p.Infinity {
		return true
	}
	lhs := p.Y * p.Y
	rhs := p.X*p.X*p.X + c.A*p.X + c.B
	// Use an absolute tolerance because the coordinates are floats.
	const eps = 1e-9
	d := lhs - rhs
	if d < 0 {
		d = -d
	}
	return d < eps
}

// AddAffine returns p + q on c using the chord-and-tangent rule.
//
// Edge cases handled:
//   - either point is infinity: return the other
//   - p and q are inverses (same x, opposite y): return infinity
//   - p == q: use the tangent slope, and return infinity when y == 0
//   - p and q are distinct but share the same x (vertical line): infinity
func (c *Curve) AddAffine(p, q AffinePoint) AffinePoint {
	if p.Infinity {
		return q
	}
	if q.Infinity {
		return p
	}
	if p.X == q.X && p.Y == -q.Y {
		return AffineInfinity()
	}

	var lambda float64
	if p == q {
		if p.Y == 0 {
			return AffineInfinity()
		}
		lambda = (3*p.X*p.X + c.A) / (2 * p.Y)
	} else {
		lambda = (q.Y - p.Y) / (q.X - p.X)
	}

	x3 := lambda*lambda - p.X - q.X
	y3 := lambda*(p.X-x3) - p.Y
	return AffinePoint{X: x3, Y: y3}
}

// ScalarMulAffine returns k*p using double-and-add.
//
// Negative k is supported by negating the result of |k|*p. The point at
// infinity is returned for k == 0 or when p is infinity.
func (c *Curve) ScalarMulAffine(p AffinePoint, k int) AffinePoint {
	if k == 0 || p.Infinity {
		return AffineInfinity()
	}
	if k < 0 {
		neg := AffinePoint{X: p.X, Y: -p.Y}
		return c.ScalarMulAffine(neg, -k)
	}

	result := AffineInfinity()
	addend := p
	for k > 0 {
		if k&1 == 1 {
			result = c.AddAffine(result, addend)
		}
		addend = c.AddAffine(addend, addend)
		k >>= 1
	}
	return result
}

// ProjectivePoint is a point in Jacobian coordinates (X, Y, Z). The affine
// point is (X/Z^2, Y/Z^3). The identity is represented by Z == 0; in that case
// X and Y are ignored.
type ProjectivePoint struct {
	X, Y, Z float64
}

// ProjectiveInfinity returns the identity in Jacobian coordinates.
func ProjectiveInfinity() ProjectivePoint {
	return ProjectivePoint{X: 0, Y: 1, Z: 0}
}

// IsInfinity reports whether p is the identity.
func (p ProjectivePoint) IsInfinity() bool {
	return p.Z == 0
}

// ToAffine converts p to affine coordinates. The point at infinity maps to
// AffineInfinity.
func (p ProjectivePoint) ToAffine() AffinePoint {
	if p.IsInfinity() {
		return AffineInfinity()
	}
	zInv := 1 / p.Z
	zInv2 := zInv * zInv
	zInv3 := zInv2 * zInv
	return AffinePoint{X: p.X * zInv2, Y: p.Y * zInv3}
}

// FromAffine converts an affine point to Jacobian coordinates. The point at
// infinity maps to ProjectiveInfinity.
func FromAffine(p AffinePoint) ProjectivePoint {
	if p.Infinity {
		return ProjectiveInfinity()
	}
	return ProjectivePoint{X: p.X, Y: p.Y, Z: 1}
}

// AddProjective returns p + q on c using Jacobian coordinates.
//
// The formulas are the standard ones from Cohen-Miyaji-Ono. They avoid
// inversion during the whole computation, which is the point of using
// projective coordinates in the first place.
func (c *Curve) AddProjective(p, q ProjectivePoint) ProjectivePoint {
	if p.IsInfinity() {
		return q
	}
	if q.IsInfinity() {
		return p
	}

	z1z1 := p.Z * p.Z
	z2z2 := q.Z * q.Z
	u1 := p.X * z2z2
	u2 := q.X * z1z1
	s1 := p.Y * q.Z * z2z2
	s2 := q.Y * p.Z * z1z1

	if u1 == u2 {
		if s1 != s2 {
			// p == -q, so the sum is the identity.
			return ProjectiveInfinity()
		}
		// p == q, use doubling.
		return c.doubleProjective(p)
	}

	h := u2 - u1
	i := (2 * h) * (2 * h)
	j := h * i
	r := 2 * (s2 - s1)
	v := u1 * i

	x3 := r*r - j - 2*v
	y3 := r*(v-x3) - 2*s1*j
	z3 := ((p.Z + q.Z) * (p.Z + q.Z) - z1z1 - z2z2) * h

	return ProjectivePoint{X: x3, Y: y3, Z: z3}
}

// doubleProjective doubles p. It is unexported because callers should use
// AddProjective with p == q; keeping it private avoids two public entry points
// with different edge-case behaviour.
func (c *Curve) doubleProjective(p ProjectivePoint) ProjectivePoint {
	if p.IsInfinity() || p.Y == 0 {
		return ProjectiveInfinity()
	}

	ysq := p.Y * p.Y
	s := 4 * p.X * ysq
	m := 3*p.X*p.X + c.A*p.Z*p.Z*p.Z*p.Z
	t := m*m - 2*s
	x3 := t
	y3 := m*(s-t) - 8*ysq*ysq
	z3 := 2 * p.Y * p.Z
	return ProjectivePoint{X: x3, Y: y3, Z: z3}
}

// ScalarMulProjective returns k*p using double-and-add in Jacobian
// coordinates. Negative k is handled by negating p first.
func (c *Curve) ScalarMulProjective(p ProjectivePoint, k int) ProjectivePoint {
	if k == 0 || p.IsInfinity() {
		return ProjectiveInfinity()
	}
	if k < 0 {
		neg := ProjectivePoint{X: p.X, Y: -p.Y, Z: p.Z}
		return c.ScalarMulProjective(neg, -k)
	}

	result := ProjectiveInfinity()
	addend := p
	for k > 0 {
		if k&1 == 1 {
			result = c.AddProjective(result, addend)
		}
		addend = c.AddProjective(addend, addend)
		k >>= 1
	}
	return result
}

// String formats an affine point for debugging.
func (p AffinePoint) String() string {
	if p.Infinity {
		return "Infinity"
	}
	return fmt.Sprintf("(%g, %g)", p.X, p.Y)
}

// String formats a projective point for debugging.
func (p ProjectivePoint) String() string {
	if p.IsInfinity() {
		return "Infinity"
	}
	return fmt.Sprintf("(%g, %g, %g)", p.X, p.Y, p.Z)
}
