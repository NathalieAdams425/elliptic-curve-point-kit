# elliptic-curve-point-kit

Affine and projective point addition and scalar multiplication on short Weierstrass curves `y^2 = x^3 + a*x + b` over `float64`.

```go
package main

import (
	"fmt"

	ecp "github.com/example/elliptic-curve-point-kit"
)

func main() {
	c := ecp.NewCurve(-1, 1)
	p := ecp.AffinePoint{X: 0, Y: 1}
	q := ecp.AffinePoint{X: 1, Y: 1}

	sum := c.AddAffine(p, q)
	fmt.Println(sum) // (-1, -1)

	triple := c.ScalarMulAffine(p, 3)
	fmt.Println(triple) // (56, 419)

	// The same computation in Jacobian coordinates, which avoids the field
	// inversion that affine addition needs.
	pj := ecp.FromAffine(p)
	qj := ecp.FromAffine(q)
	fmt.Println(c.AddProjective(pj, qj).ToAffine()) // (-1, -1)
}
```

## Why this exists

The affine chord-and-tangent formulas are short and easy to read, but every addition performs a division, which is the expensive operation in a finite field. Projective (Jacobian) coordinates remove that division: the coordinates are scaled so that the affine point is `(X/Z^2, Y/Z^3)`, and the group law becomes a sequence of multiplications and additions. `AddProjective` and `ScalarMulProjective` implement that form, and `ToAffine`/`FromAffine` convert between the two views.

The trade-off is that projective arithmetic is harder to read and has its own special cases (the identity is `Z == 0`, and doubling a point with `Y == 0` returns the identity). This library keeps both representations so the affine code can be used as a reference for the projective code, and the tests compare the two on the same inputs.

## The awkward edge

The `AffinePoint` zero value is `(0, 0)`, which is a perfectly ordinary point on many curves. It is not the identity. The identity is marked by the `Infinity` field; use `AffineInfinity()` to construct it and `IsInfinity()` to test for it. If you build an `AffinePoint` literal and forget the flag, the arithmetic will treat `(0, 0)` as a real point.

Coordinates are `float64`, so this package is for exploration and teaching, not for cryptography. Equality of points is approximate; the tests use an absolute tolerance of `1e-9`.

## API

- `type Curve struct { A, B float64 }`, `NewCurve(a, b float64) *Curve`
- `(*Curve) Discriminant() float64`, `(*Curve) IsNonSingular() bool`
- `type AffinePoint struct { X, Y float64; Infinity bool }`, `AffineInfinity() AffinePoint`
- `(AffinePoint) IsInfinity() bool`, `(AffinePoint) OnCurve(*Curve) bool`, `(AffinePoint) String() string`
- `(*Curve) AddAffine(p, q AffinePoint) AffinePoint`, `(*Curve) ScalarMulAffine(p AffinePoint, k int) AffinePoint`
- `type ProjectivePoint struct { X, Y, Z float64 }`, `ProjectiveInfinity() ProjectivePoint`
- `(ProjectivePoint) IsInfinity() bool`, `(ProjectivePoint) ToAffine() AffinePoint`, `(ProjectivePoint) String() string`
- `FromAffine(p AffinePoint) ProjectivePoint`
- `(*Curve) AddProjective(p, q ProjectivePoint) ProjectivePoint`, `(*Curve) ScalarMulProjective(p ProjectivePoint, k int) ProjectivePoint`
