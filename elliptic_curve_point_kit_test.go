package elliptic_curve_point_kit

import (
	"math"
	"testing"
)

const eps = 1e-9

func almostEqual(a, b float64) bool {
	return math.Abs(a-b) < eps
}

func pointsAlmostEqual(a, b AffinePoint) bool {
	if a.Infinity || b.Infinity {
		return a.Infinity == b.Infinity
	}
	return almostEqual(a.X, b.X) && almostEqual(a.Y, b.Y)
}

func projAlmostEqual(a, b ProjectivePoint) bool {
	if a.IsInfinity() || b.IsInfinity() {
		return a.IsInfinity() == b.IsInfinity()
	}
	return pointsAlmostEqual(a.ToAffine(), b.ToAffine())
}

// curve y^2 = x^3 - x + 1 has discriminant 4*(-1)^3 + 27*1 = 23, non-singular.
func testCurve() *Curve {
	return NewCurve(-1, 1)
}

// A point that satisfies y^2 = x^3 - x + 1: x=0, y=1.
func p0() AffinePoint { return AffinePoint{X: 0, Y: 1} }

// x=1, y=1: 1 = 1 - 1 + 1 = 1.
func p1() AffinePoint { return AffinePoint{X: 1, Y: 1} }

func TestCurveDiscriminant(t *testing.T) {
	c := testCurve()
	if !almostEqual(c.Discriminant(), 23) {
		t.Fatalf("discriminant = %v, want 23", c.Discriminant())
	}
	if !c.IsNonSingular() {
		t.Fatal("curve should be non-singular")
	}
}

func TestSingularCurve(t *testing.T) {
	// y^2 = x^3 has a cusp at the origin.
	c := NewCurve(0, 0)
	if c.IsNonSingular() {
		t.Fatal("y^2 = x^3 should be singular")
	}
}

func TestOnCurve(t *testing.T) {
	c := testCurve()
	if !p0().OnCurve(c) {
		t.Fatal("p0 should be on the curve")
	}
	if !p1().OnCurve(c) {
		t.Fatal("p1 should be on the curve")
	}
	bad := AffinePoint{X: 0, Y: 2}
	if bad.OnCurve(c) {
		t.Fatal("point (0,2) should not be on the curve")
	}
	if !AffineInfinity().OnCurve(c) {
		t.Fatal("infinity should be on the curve")
	}
}

func TestAddAffineIdentity(t *testing.T) {
	c := testCurve()
	got := c.AddAffine(p0(), AffineInfinity())
	if !pointsAlmostEqual(got, p0()) {
		t.Fatalf("p0 + O = %v, want %v", got, p0())
	}
	got = c.AddAffine(AffineInfinity(), p1())
	if !pointsAlmostEqual(got, p1()) {
		t.Fatalf("O + p1 = %v, want %v", got, p1())
	}
	got = c.AddAffine(AffineInfinity(), AffineInfinity())
	if !got.Infinity {
		t.Fatalf("O + O = %v, want infinity", got)
	}
}

func TestAddAffineInverse(t *testing.T) {
	c := testCurve()
	neg := AffinePoint{X: p0().X, Y: -p0().Y}
	got := c.AddAffine(p0(), neg)
	if !got.Infinity {
		t.Fatalf("p + (-p) = %v, want infinity", got)
	}
}

func TestAddAffineDoubling(t *testing.T) {
	c := testCurve()
	got := c.AddAffine(p0(), p0())
	// Tangent at (0,1): slope = (3*0^2 - 1) / (2*1) = -0.5.
	// x3 = 0.25, y3 = -0.5*(0 - 0.25) - 1 = 0.125 - 1 = -0.875.
	want := AffinePoint{X: 0.25, Y: -0.875}
	if !pointsAlmostEqual(got, want) {
		t.Fatalf("2*p0 = %v, want %v", got, want)
	}
	if !got.OnCurve(c) {
		t.Fatalf("2*p0 = %v is not on the curve", got)
	}
}

func TestAddAffineDistinct(t *testing.T) {
	c := testCurve()
	got := c.AddAffine(p0(), p1())
	// slope = (1-1)/(1-0) = 0; x3 = -1, y3 = 0*(0+1) - 1 = -1.
	want := AffinePoint{X: -1, Y: -1}
	if !pointsAlmostEqual(got, want) {
		t.Fatalf("p0 + p1 = %v, want %v", got, want)
	}
	if !got.OnCurve(c) {
		t.Fatalf("p0 + p1 = %v is not on the curve", got)
	}
}

func TestAddAffineVertical(t *testing.T) {
	c := testCurve()
	neg := AffinePoint{X: p1().X, Y: -p1().Y}
	got := c.AddAffine(p1(), neg)
	if !got.Infinity {
		t.Fatalf("p1 + (-p1) = %v, want infinity", got)
	}
}

func TestScalarMulAffine(t *testing.T) {
	c := testCurve()
	p := p0()

	got := c.ScalarMulAffine(p, 0)
	if !got.Infinity {
		t.Fatalf("0*p = %v, want infinity", got)
	}
	got = c.ScalarMulAffine(p, 1)
	if !pointsAlmostEqual(got, p) {
		t.Fatalf("1*p = %v, want %v", got, p)
	}
	got = c.ScalarMulAffine(p, 2)
	want := AffinePoint{X: 0.25, Y: -0.875}
	if !pointsAlmostEqual(got, want) {
		t.Fatalf("2*p = %v, want %v", got, want)
	}
	got = c.ScalarMulAffine(p, 3)
	// 3*p = 2*p + p. 2*p = (0.25, -0.875). slope = (-0.875-1)/(0.25-0) = -7.5.
	// x3 = 56.25 - 0.25 = 56; y3 = -7.5*(0 - 56) - 1 = 420 - 1 = 419.
	want = AffinePoint{X: 56, Y: 419}
	if !pointsAlmostEqual(got, want) {
		t.Fatalf("3*p = %v, want %v", got, want)
	}
}

func TestScalarMulAffineNegative(t *testing.T) {
	c := testCurve()
	p := p0()
	got := c.ScalarMulAffine(p, -1)
	want := AffinePoint{X: 0, Y: -1}
	if !pointsAlmostEqual(got, want) {
		t.Fatalf("-1*p = %v, want %v", got, want)
	}
}

func TestProjectiveRoundTrip(t *testing.T) {
	p := p0()
	proj := FromAffine(p)
	back := proj.ToAffine()
	if !pointsAlmostEqual(back, p) {
		t.Fatalf("round trip = %v, want %v", back, p)
	}
	if !proj.ToAffine().OnCurve(testCurve()) {
		t.Fatal("round-tripped point not on curve")
	}
}

func TestAddProjectiveIdentity(t *testing.T) {
	c := testCurve()
	p := FromAffine(p0())
	got := c.AddProjective(p, ProjectiveInfinity())
	if !projAlmostEqual(got, p) {
		t.Fatalf("p + O = %v, want %v", got, p)
	}
	got = c.AddProjective(ProjectiveInfinity(), p)
	if !projAlmostEqual(got, p) {
		t.Fatalf("O + p = %v, want %v", got, p)
	}
}

func TestAddProjectiveMatchesAffine(t *testing.T) {
	c := testCurve()
	pairs := [][2]AffinePoint{
		{p0(), p0()},
		{p0(), p1()},
		{p1(), p1()},
	}
	for _, pair := range pairs {
		a := c.AddAffine(pair[0], pair[1])
		b := c.AddProjective(FromAffine(pair[0]), FromAffine(pair[1])).ToAffine()
		if !pointsAlmostEqual(a, b) {
			t.Fatalf("affine %v vs projective %v for %v + %v", a, b, pair[0], pair[1])
		}
	}
}

func TestScalarMulProjectiveMatchesAffine(t *testing.T) {
	c := testCurve()
	p := p0()
	for k := 0; k <= 5; k++ {
		a := c.ScalarMulAffine(p, k)
		b := c.ScalarMulProjective(FromAffine(p), k).ToAffine()
		if !pointsAlmostEqual(a, b) {
			t.Fatalf("k=%d affine %v vs projective %v", k, a, b)
		}
	}
}

func TestAddProjectiveInverse(t *testing.T) {
	c := testCurve()
	p := FromAffine(p0())
	neg := ProjectivePoint{X: p.X, Y: -p.Y, Z: p.Z}
	got := c.AddProjective(p, neg)
	if !got.IsInfinity() {
		t.Fatalf("p + (-p) = %v, want infinity", got)
	}
}

func TestDoubleProjectiveYZero(t *testing.T) {
	c := testCurve()
	// A point with y = 0 on this curve would need x^3 - x + 1 = 0, which has a
	// real root near -1.3247. Find it numerically to exercise the y == 0 branch.
	x := -1.324717957244746
	p := FromAffine(AffinePoint{X: x, Y: 0})
	got := c.AddProjective(p, p)
	if !got.IsInfinity() {
		t.Fatalf("doubling a y=0 point = %v, want infinity", got)
	}
}

func TestStringFormatting(t *testing.T) {
	p := AffinePoint{X: 0, Y: 1}
	if p.String() != "(0, 1)" {
		t.Fatalf("affine string = %q", p.String())
	}
	if AffineInfinity().String() != "Infinity" {
		t.Fatalf("affine infinity string = %q", AffineInfinity().String())
	}
	pp := ProjectivePoint{X: 0, Y: 1, Z: 1}
	if pp.String() != "(0, 1, 1)" {
		t.Fatalf("projective string = %q", pp.String())
	}
	if ProjectiveInfinity().String() != "Infinity" {
		t.Fatalf("projective infinity string = %q", ProjectiveInfinity().String())
	}
}
