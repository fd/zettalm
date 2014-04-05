package lsq

// copyright (c) 2014 by Jason E. Aten

import (
	"errors"
	"fmt"
	"math"
	"strings"
)

//  Project Zettalm, source file: lsq.go
//
//  * fits linear models to data *
//
//  This does basic linear regression, otherwise know as least-squares fiting.
//
//  * but: handles really big data *
//
//  This particular solution to fitting linear models can handle infinite data,
//  and that's the point.
//  Again: there a tons of model fitting packages, but the point of interest here
//  is that the online QR decomposition algorithm can handle Zetta-Bytes of observations,
//  *actually unlimited* rows of data. It needs space only proportional to
//  O(p^2) for p variables. It only ever has to look at each row once. The original
//  Fortran90 did not handle more than one dependent variable, but the current
//  code does. Again, observations of x and y variables are not stored in memory.
//  Weighted observations are supported, including negative weights to delete
//  cases.
//
//  * is still quite fast *
//
//  You are certainly most likely to be IO-bound. For CPU time, I benchmarked
//  the current code at 10x the performance of R's linear model solver.
//  Asymptodically, for p variables (both x and y variables), each additional
//  row of observations incurs O(p^2) time. To process n rows will take
//  you time O(n*p*p) assuming n is > p.
//
//  * enhancements over the fortran code *
//
//  I've added the ability to handle multiple y-variables at once, as well
//  as online computation of means and standard deviations for all x and y variables.
//  Also I've fixed bugs that were present when using the downdating functionality
//  to remove rows from the set.
//
//  * online, but can be a moving window too *
//
//  Even though it only ever looks at any row once, rows are weighted, and
//  hence rows can be deleted, if you wish, from the data set considered
//  by adding the same row again but using a weight of -1 (instead of the default +1 weight)
//  on the second call to Includ().

//
//  Origins:
//
//  The structure and algorithm at the heart of this file is derived from
//  Alan Miller's AS 274 Algorithm and his Fortran90 source code for that publication.
//  AS 274 in turn builds upon W. Morven Gentleman's AG 75 algorithm.
//  http://www.jstor.org/stable/2347147
//
//  Additional information and the Fortran90 source is available at
//  the mirror, http://jblevins.org/mirror/amiller/lsq.f90.
//  The web page there indicates that Alan Miller's version of the code
//  has been placed in the public domain. Nonetheless, extensive credit is due to him
//  for this rigorous and numerically exacting algorithm. As a tribute, we
//  name the central struct the MillerLSQ struct. LSQ is short for least
//  squares, the method of regression implemented here.
//
//
//  Relevant commentary from lsq.f90 at http://jblevins.org/mirror/amiller/lsq.f90
//  follows, and has been edited for relevance to the current package.
//
//
//  Module for unconstrained linear least-squares calculations.
//  The algorithm is suitable for updating LS calculations as more
//  data are added.   This is sometimes called recursive estimation.
//
//  Based upon Applied Statistics algorithm AS 274.
//
//  The PUBLIC variables of struct MillerLSQ:
//
//  Nobs    = the number of observations processed to date.
//  Ncol    = the total number of dependent (x) variables,
//            including one for the constant. In this version of the code,
//            a constant or intercept term is always fitted.
//
//  R_dim   = the dimension of the upper triangular matrix R = Ncol*(Ncol-1)/2
//  Vorder  = an integer vector storing the current order of the variables
//            in the QR-factorization.   The initial order is 0, 1, 2, ...
//
//  Initialized = a logical variable which indicates whether space has
//                been allocated for various arrays.
//  Tol_set = a logical variable which is set when subroutine Tolset() has
//            been called to calculate tolerances for use in testing for
//            singularities.
//  Rss_set = a logical variable indicating whether residual sums of squares
//            are available and usable.
//  D       = array of row multipliers for the Cholesky factorization.
//            The factorization is X = Q * sqrt(D) * R where Q is an ortho-
//            normal matrix which is NOT stored, D is a diagonal matrix
//            whose diagonal elements are stored in array d, and R is an
//            upper-triangular matrix with 1's as its diagonal elements.
//  Rhs     = vector of RHS projections (after scaling by sqrt(D)).
//            Thus Q'y = sqrt(D) * Rhs
//
//  R       = the upper-triangular matrix R.   The upper triangle only,
//            excluding the implicit 1's on the diagonal, are stored by
//            rows. This is where the R of the QR decomposition is maintained.
//            R is also called the Orthogonal reduction in AS 75.1, or
//            the Cholesky factorization.
//
//  Tol     = array of tolerances used in testing for singularities.
//  Rss     = array of residual sums of squares, one for each Y variable.
//            Rss[wycol][i] is the residual
//            sum of squares with the first i variables in the model for the wycol-th
//            y-variable.
//            By changing the order of variables, the residual sums of
//            squares can be found for all possible subsets of the variables.
//
//            The residual sum of squares with NO variables in the model,
//            that is the total sum of squares of the first column of y-values, can be
//            calculated as Rss[0][0] + D[0]*Rhs[0][0]^2.
//            Since we always fit a constant, by itself Rss[0][0] is the sum of squares of
//            (y - ybar) where ybar is the average value of y.
//
//  Sserr   = residual sum of squares with all of the variables included, one
//             for each y-variable.
//  Row_ptr = array of indices of first elements in each row of R.
//
//--------------------------------------------------------------------------

type NanHandling uint32

const (
	NAN_TO_ZERO  NanHandling = 0
	NAN_OMIT_ROW             = 1
)

// change assignments from 1 (Fortan) based to 0 (Go) based by doing -Adj
// inside many of [] access
const Adj = 1

type MillerLSQ struct {
	XStats SdTracker
	YStats SdTracker

	Nobs           int64   // number of observations included
	AccumWeightSum float64 // accumulate negative weights too, to track downdating
	Ncol           int
	R_dim          int
	Nxvar          int // not counting Constant (intercept) coefficient. Nxvar = Ncol -1.
	Nyvar          int // number of y-target variables

	UseMeanSd bool // Apply x' = (x - mean)/sd normalization
	Xmean     []float64
	Xsd       []float64
	Ymean     []float64
	Ysd       []float64

	Initialized bool
	Tol_set     bool
	Rss_set     bool

	Vsmall float64
	Sserr  []float64
	Toly   float64

	Vorder  []int
	Row_ptr []int

	R []float64 // row-major, upper triangular without the diagonal (implicit 1's)

	D   []float64   // the diagonal
	Rhs [][]float64 // related to the beta, a list of column vectors, each vector a Rhs for a different y-target.
	Tol []float64
	Rss [][]float64 // one for each y

	Rinv     []float64
	Dim_rinv int

	// Curxrow: copy incomming xrow here, and we add the 1 in the constant
	//          position so users don't need to worry about it.
	Curxrow []float64
	Curyrow []float64
}

func (m *MillerLSQ) SetMeanSd(xmean []float64, xsd []float64, ymean []float64, ysd []float64) {
	if m.UseMeanSd {
		panic("can only call SetMeanSd once.")
	}
	if m.Nobs > 0 {
		panic("can only call SetMeanSd *before* calling Includ()")
	}
	if len(ymean) != m.Nyvar {
		panic(fmt.Sprintf("ymean(%v) did not agree with m.Nyvar(%v)", ymean, m.Nyvar))
	}
	if len(xmean) != m.Nxvar {
		panic(fmt.Sprintf("xmean(%v) did not agree with m.Nxvar(%v)", xmean, m.Nxvar))
	}
	copy(m.Xmean, xmean)
	copy(m.Xsd, xsd)
	copy(m.Ymean, ymean)
	copy(m.Ysd, ysd)
	m.UseMeanSd = true
}

func DeepCopy(a []float64) []float64 {
	res := make([]float64, len(a))
	copy(res, a)
	return res
}

// NewMillerLSQ(): the constructor
//
// nxvar should include the count of all x, but not the Constant 1 column.
// nyvar should count the number of y-target variables (dependent variables).
//  while beginning usage will typically only have one y-target, supporting
//  many targets without having to re-run through a large data set is a boon.
//
// there is always an intercept fit, so don't count it
// nxvar = number of x columns (indep variables)
// nyvar = number of y columns (dependent variable(s))
//
//
func NewMillerLSQ(nxvar int, nyvar int) *MillerLSQ {
	//fmt.Printf("NewMillerLSQ called with nxvar=%v   and nyvar=%v\n", nxvar, nyvar)
	m := &MillerLSQ{}

	m.Nobs = 0
	m.AccumWeightSum = 0
	m.Ncol = nxvar + 1 // 1 for the constant term (intercept)
	m.Nxvar = nxvar
	m.Nyvar = nyvar

	m.XStats.ZeroTracker(m.Nxvar)
	m.YStats.ZeroTracker(m.Nyvar)

	m.R_dim = m.Ncol * (m.Ncol - 1) / 2

	m.UseMeanSd = false
	m.Xmean = make([]float64, m.Nxvar)
	m.Xsd = make([]float64, m.Nxvar)
	m.Ymean = make([]float64, m.Nyvar)
	m.Ysd = make([]float64, m.Nyvar)

	m.Initialized = true
	m.Tol_set = false
	m.Rss_set = false

	m.Vsmall = 1e-12
	m.Sserr = make([]float64, m.Nyvar)
	m.Toly = 0.0

	m.Vorder = make([]int, m.Ncol)
	m.Row_ptr = make([]int, m.Ncol)

	m.R = make([]float64, m.R_dim)

	m.D = make([]float64, m.Ncol)
	m.Rhs = make([][]float64, m.Nyvar)
	for i := range m.Rhs {
		m.Rhs[i] = make([]float64, m.Ncol)
	}
	m.Tol = make([]float64, m.Ncol)
	m.Rss = make([][]float64, m.Nyvar)
	for i := range m.Rss {
		m.Rss[i] = make([]float64, m.Ncol)
	}

	for i := 1; i <= m.Ncol; i++ {
		m.Vorder[i-Adj] = i - 1
	}

	// m.Row_ptr[i] is the position of element R(i,i+1) in array m.R.
	// Since the diagonal is implicit 1's at R(i,i), this is the
	// first relevant/interesting entry.

	// Row_ptr entries are accurately 0-based or 1-based depending on Adj
	m.Row_ptr[1-Adj] = 1 - Adj
	for i := 2; i <= m.Ncol-1; i++ {
		m.Row_ptr[i-Adj] = m.Row_ptr[(i-1)-Adj] + m.Ncol - i + 1
	}

	// the m.Row_ptr[m.Ncol-Adj] entry should never be used; it doesn't exist since
	// there is only a 1 on the diagnonal in the R(Ncol, Ncol) 1-based position.
	// so here we make sure this last entry is out-of-bounds for m.R
	m.Row_ptr[m.Ncol-Adj] = m.R_dim + 1

	m.Curxrow = make([]float64, m.Ncol)
	m.Curyrow = make([]float64, m.Nyvar)

	return m
} // end startup

func StdFormat6dec(x float64) string {
	return fmt.Sprintf(" %.6f", x)
}
func StdFormat5dec(x float64) string {
	return fmt.Sprintf(" %.5f", x)
}
func StdFormat2dec(x float64) string {
	return fmt.Sprintf(" %.2f", x)
}
func StdFormat1dec(x float64) string {
	return fmt.Sprintf(" %.1f", x)
}
func StdFormat0dec(x float64) string {
	return fmt.Sprintf(" %.0f", x)
}
func StdFormatE(x float64) string {
	return fmt.Sprintf(" %e", x)
}

// linesep defaults to "\n" unless supplied differently
func UpperTriToString(dat []float64, ncol int, linesep string, formatFunc func(float64) string) string {
	if linesep == "" {
		linesep = "\n"
	}
	var i, j, pos int
	s := "["
	pos = 0
	for i = 1; i <= ncol; i++ {

		for j = 1; j <= ncol; j++ {
			if j < i {
				s += formatFunc(0.0)
			} else if j == i {
				s += formatFunc(1.0)
			} else {
				s += formatFunc(dat[pos])
				pos++
			}
		}

		if i < ncol {
			s += " " + linesep
		} else {
			s += "]"
		}
	}

	return s
}

func fEquals(a, b float64) bool {
	if math.Abs(a-b) < 1e-10 {
		return true
	}
	return false
}

func EpsEquals(a, b, eps float64) bool {
	if math.Abs(a-b) < eps {
		return true
	}
	return false
}

func fSliceEqual(a, b []float64) bool {
	if len(a) != len(b) {
		return false
	}
	if len(a) == 0 {
		return true
	}
	for i := range a {
		if !fEquals(a[i], b[i]) {
			return false
		}
	}
	return true
}

func EpsSliceEqual(a, b []float64, eps float64) bool {
	if len(a) != len(b) {
		return false
	}
	if len(a) == 0 {
		return true
	}
	for i := range a {
		if !EpsEquals(a[i], b[i], eps) {
			return false
		}
	}
	return true
}

func Int64SliceEqual(a, b []int64, eps float64) bool {
	if len(a) != len(b) {
		return false
	}
	if len(a) == 0 {
		return true
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func IntSliceEqual(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	if len(a) == 0 {
		return true
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func StringSliceEqual(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	if len(a) == 0 {
		return true
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// Includ(): the main method for adding new rows of data and thus updating the
//  QR factorization.
//
//  reference: Golub and Van Loan, Matrix Computations, Section 5.1.13 on Fast Givens Transformations
//             page 219 in the 3rd edition.
//
//   consider: http://en.wikipedia.org/wiki/Givens_rotation #Stability section, at some point.
//
//     ALGORITHM AS75.1  APPL. STATIST. (1974) VOL.23, NO. 3
//
//     Calling this routine updates D, R, RHS and SSERR by the
//     inclusion of xrow, yelem = each of yrow member in turn, with the specified weight.
//
//   iff row was ommitted due to nanapproach == NAN_OMIT_ROW, we return true
func (m *MillerLSQ) Includ(weight float64, xrow []float64, yrow []float64, nanapproach NanHandling) (rowOmitted bool) {

	// xi having NaN is messing us up. Zero them out. This takes care of
	//  nanapproach == NAN_TO_ZERO
	xHasNaN := NanToZero(xrow)
	yHasNaN := NanToZero(yrow)

	// but in order to match with R regression (default) linear model,
	// we also need to be able to omit rows with NaN in them.
	// Not yet implemented: replace missing value with xmean?
	if nanapproach == NAN_OMIT_ROW {
		if xHasNaN || yHasNaN {
			// returning early omits the row
			return true
		}
	}

	//fmt.Printf("Includ() called with xrow = %v   and yrow = %v\n", xrow, yrow)

	if len(yrow) != m.Nyvar {
		panic(fmt.Sprintf("len(yrow) == %v did not match m.Nyvar == %v", len(yrow), m.Nyvar))
	}

	if len(xrow) != m.Nxvar {
		panic(fmt.Sprintf("len(xrow) == %v did not match m.Nxvar == %v", len(xrow), m.Nxvar))
	}

	m.Curxrow[0] = 1.0
	copy(m.Curxrow[1:], xrow) // Curxrow is 1 longer than xrow
	copy(m.Curyrow, yrow)

	var i, k, nextr int
	var w, y, xi, di, wxi, dpi, cbar, sbar, xk float64

	w = weight

	// don't bother adding a zero-weighted x vector
	if w == 0 || math.Abs(w) < m.Vsmall {
		return
	} else {
		m.AccumWeightSum += weight
	}

	// track our mean and sd
	m.XStats.AddObs(xrow, weight)
	m.YStats.AddObs(yrow, weight)

	if m.UseMeanSd {
		for j := range m.Curyrow {
			if m.Ysd[j] != 0.0 {
				m.Curyrow[j] = (m.Curyrow[j] - m.Ymean[j]) / m.Ysd[j]
			}
		}

		// use bump to skip past the 1.0 constant in xrow
		bump := 1
		for j := range xrow {
			if m.Xsd[j] != 0.0 {
				m.Curxrow[j+bump] = (m.Curxrow[j+bump] - m.Xmean[j]) / m.Xsd[j]
			}
		}
	}

	//y = m.Curyrow[0]

	m.Nobs = m.Nobs + 1
	m.Rss_set = false
	nextr = 1 - Adj
	for i = 1; i <= m.Ncol; i++ {

		//     Skip unnecessary transformations.   Test on exact zeroes must be
		//     used or stability can be destroyed.

		xi = m.Curxrow[i-Adj]

		if math.Abs(xi) < m.Vsmall {
			nextr = nextr + m.Ncol - i // okay, no Adj needed.
		} else {
			di = m.D[i-Adj]
			wxi = w * xi
			dpi = di + wxi*xi

			if math.IsNaN(dpi) {
				msg := "detected IsNaN(dpi) in Includ()"
				panic(msg)
			}

			// can't divide by dpi when dpi is 0
			if fEquals(dpi, 0.0) {
				cbar = 0.0
				sbar = 0.0
			} else {
				cbar = di / dpi
				sbar = wxi / dpi
			}
			w = cbar * w
			m.D[i-Adj] = dpi
			for k = i + 1; k <= m.Ncol; k++ {
				xk = m.Curxrow[k-Adj]
				m.Curxrow[k-Adj] = xk - xi*m.R[nextr]  // no Adj needed for [nextr]
				m.R[nextr] = cbar*m.R[nextr] + sbar*xk // no Adj needed for [nextr]
				//fmt.Printf("R is %s\n", m.RtoString(2))

				nextr = nextr + 1
			}
			// update all the m.Rhs[]
			for yi := range m.Curyrow {
				y = m.Curyrow[yi]
				xk = y
				y = xk - xi*m.Rhs[yi][i-Adj]
				m.Rhs[yi][i-Adj] = cbar*m.Rhs[yi][i-Adj] + sbar*xk

				m.Curyrow[yi] = y
			}
		}
		//fmt.Printf("during Includ(), pass %d: m.Rhs is %v\n", i, m.Rhs)
	} // end i over 1..Ncol

	//     Y * math.Sqrt(W) is now equal to the Brown, Durbin & Evans recursive
	//     residual.

	for yi := range m.Curyrow {
		y = m.Curyrow[yi]
		m.Sserr[yi] = m.Sserr[yi] + w*y*y
	}

	return
} // end Includ()

//    Regcf(): This returns the least-squares regression coefficients in array beta.
//
//     wxcol = which of the x-variables you want in your model. The constant/intercept
//      term is always included and is the 0-th term. So the contents of wxcol should
//      not include zero but should start at 1, 2, 3, ... for instance, if a model
//      using the first three x-variables is desired.
//
//     wycol = which single y column you want to regress on. 0-based.
//
//    If singularities are detected, one or more of the regression coefficients
//    will be returned as zero, and non-nil error will be returned.
//
//    fotran notes:
//     Modified version of AS75.4 to calculate regression coefficients
//     for the first NREQ variables, given an orthogonal reduction from
//     AS75.1.

func (m *MillerLSQ) Regcf(wxcol []int, wycol int) (err error, beta []float64) {

	var nreq int = len(wxcol) + 1
	beta = make([]float64, nreq)

	if wycol < 0 {
		panic("wcol cannot be negative")
	}
	if wycol >= len(m.Rhs) {
		panic(fmt.Sprintf("out of bounds: wycol==%d was >= len(m.Rhs)==%d", wycol, m.Rhs))
	}

	// put the [xcol] right after the Constant
	/* // problem: this is changing the beta output, when it shouldn't, should it?!?.
	// to be fixed later
		var err error
		if len(wxcol) > 0 {
			err = m.Reorder(wxcol, 1)
			if err != nil {
				panic(err)
			}
		}
	*/

	var i, j, nextr int

	//     Some checks.

	var res error = nil
	if nreq < 1 || nreq > m.Ncol {
		return errors.New(fmt.Sprintf("Regcf() call error: nreq==%v out of range.", nreq)), beta
	}

	if !m.Tol_set {
		m.Tolset(1e-12)
	}

	//fmt.Printf("D = %v\n", m.D)
	//Rstr := "\n" + UpperTriToString(m.R, m.Ncol, "\n", StdFormat6dec)
	//fmt.Printf("R = %s\n", Rstr)

	//fmt.Printf("m.Rhs = %v\n", m.Rhs)

	badList := make([]int, 0)
	for i = nreq; i >= 1; i-- {
		if math.Sqrt(m.D[i-Adj]) < m.Tol[i-Adj] {
			beta[i-Adj] = 0.0
			m.D[i-Adj] = 0.0

			badList = append(badList, i-Adj)
			res = errors.New(fmt.Sprintf("Regcf() error, singular variable(s) at locations: %v.", badList))
		} else {
			beta[i-Adj] = m.Rhs[wycol][i-Adj]
			nextr = m.Row_ptr[i-Adj]
			for j = i + 1; j <= nreq; j++ {
				beta[i-Adj] = beta[i-Adj] - m.R[nextr]*beta[j-Adj]
				nextr = nextr + 1
			}
		}
	}

	return res, beta
}

var EpsilonFloat64 float64 = math.Nextafter(1.0, 2.0) - 1.0 // 2.220446049250313e-16

func (m *MillerLSQ) Tolset(eps float64) {

	//     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

	//     Sets up array TOL for testing for zeroes in an orthogonal
	//     reduction formed using AS75.1.

	// REAL (dp), INTENT(IN), OPTIONAL :: eps

	//     Unless the argument eps is set, it is assumed that the input data are
	//     recorded to full machine accuracy.   This is often not the case.
	//     If, for instance, the data are recorded to `single precision' of about
	//     6-7 significant decimal digits, then singularities will not be detected.
	//     It is suggested that in this case eps should be set equal to
	//     10.0 * EPSILON(1.0)
	//     If the data are recorded to say 4 significant decimals, then eps should
	//     be set to 1.0E-03
	//     The above comments apply to the predictor variables, not to the
	//     dependent variable.

	//     Correction - 19 August 2002
	//     When negative weights are used, it is possible for an alement of D
	//     to be negative.

	var col, row, pos int
	var eps1, total float64
	work := make([]float64, m.Ncol)

	//     EPS is a machine-dependent constant.
	/*
		// EPSILON(ten) == 2.22E-16, so ten * EPSILON(ten) = 2.22e-15

			   if (PRESENT(eps)) {
			     eps1 = MAX(math.Abs(eps), ten * EPSILON(ten))
			   } else {
			     eps1 = ten * EPSILON(ten)
			   }
	*/
	eps1 = 10 * EpsilonFloat64
	if eps > 0 && eps > eps1 {
		eps1 = eps
	}

	//     Set tol[i] = sum of absolute values in column I of R after
	//     scaling each element by the square root of its row multiplier,
	//     multiplied by EPS1.

	for k := range work {
		work[k] = math.Sqrt(math.Abs(m.D[k])) // 0-based already, no Adj needed.
	}
	for col = 1; col <= m.Ncol; col++ {
		pos = (col - 1) - Adj
		total = work[col-Adj]
		// for row loop does nothing when col == 1, note the condition that row <= col-1
		//  and therefore R[0] when pos is 0 never gets indexed.
		for row = 1; row <= col-1; row++ {
			total = total + math.Abs(m.R[pos])*work[row-Adj] // R[pos] is correct, doesn't need Adj
			pos = pos + m.Ncol - row - 1
		}
		m.Tol[col-Adj] = eps1 * total
	}

	m.Tol_set = true
} // end Tolset

func assign(to []float64, from []float64) {
	if len(to) != len(from) {
		panic(fmt.Sprintf("to(%d) and from(%d) lengths must match in assign()", len(to), len(from)))
	}
	for i := range from {
		to[i] = from[i]
	}
}

func zero_out(x []float64) {
	for i := range x {
		x[i] = 0.0
	}
}

// -----------------------------------------------------------------
// SingularCheck()
//
//     Checks for singularities, reports, and adjusts orthogonal
//     reductions produced by AS75.1.
//
//     Correction - 19 August 2002
//     When negative weights are used, it is possible for an alement of D
//     to be negative.
//
//     Auxiliary routines called: Includ(), Tolset()
//
//     returns number of linearly dependent variables found
//      and fills in lindep with true or false
// -----------------------------------------------------------------

func (m *MillerLSQ) SingularCheck(lindep *[]bool, wycol int) int {

	if len(*lindep) != m.Ncol {
		panic(fmt.Sprintf("len(*lindep)==%d did not match m.Ncol==%d", len(*lindep), m.Ncol))
	}

	var temp, y, weight float64
	var pos, row, pos2 int
	x := make([]float64, m.Ncol)
	work := make([]float64, m.Ncol)

	var ifault int = 0

	for k := range work {
		work[k] = math.Sqrt(math.Abs(m.D[k]))
	}
	if !m.Tol_set {
		m.Tolset(1e-12)
	}

	for row = 1; row <= m.Ncol; row++ {
		temp = m.Tol[row-Adj]
		pos = m.Row_ptr[row-Adj] // pos = location of first element in row

		//     If diagonal element is near zero, set it to zero, set appropriate
		//     element of LINDEP, and use INCLUDE to augment the projections in
		//     the lower rows of the orthogonalization.

		(*lindep)[row-Adj] = false

		// 0 <= 0 is making this fire when we've seen no data. try 0 < 0 strictly.
		//if work[row-Adj] <= temp {
		if work[row-Adj] <= temp {
			(*lindep)[row-Adj] = true
			ifault++
			if row < m.Ncol {
				pos2 = pos + m.Ncol - row - 1
				zero_out(x[:row+1-Adj])
				// x[row+1:m.Ncol] = r[pos:pos2] // F90
				// assign(x[row+1-Adj:m.Ncol], m.R[pos-Adj:pos2+1-Adj]) // crashing
				slcTo := x[row+1-Adj : m.Ncol+1-Adj]
				slcFrom := m.R[pos : pos2+1] // no Adj on pos-based access to m.R[]
				//fmt.Printf("len(m.R) = %v,  pos=%d, pos2=%d\n", len(m.R), pos, pos2)
				//fmt.Printf("len(slcTo)=%d   len(slcFrom)=%d\n", len(slcTo), len(slcFrom))
				assign(slcTo, slcFrom)
				y = m.Rhs[wycol][row-Adj]
				weight = m.D[row-Adj]
				zero_out(m.R[pos : pos2+1]) // no Adj on pos-based access to m.R[]
				m.D[row-Adj] = 0.0
				m.Rhs[wycol][row-Adj] = 0.0

				//fmt.Printf("weight=%v   x=%v   y=%v\n", weight, x, y)
				m.Includ(weight, x[1:], []float64{y}, NAN_TO_ZERO)
				// INCLUD automatically increases the number
				// of cases each time it is called. compensate by decreasing m.Nobs
				m.Nobs = m.Nobs - 1
			} else {
				m.Sserr[0] = m.Sserr[0] + m.D[row-Adj]*m.Rhs[wycol][row-Adj]*m.Rhs[wycol][row-Adj]
			}
		}
	}

	return ifault
} // end sing

//
// SS()
//     Calculates partial residual sums of squares from an orthogonal
//     reduction from AS75.1.
//
func (m *MillerLSQ) SS(wycol int) {

	var i int
	var total float64

	total = m.Sserr[wycol]
	m.Rss[wycol][m.Ncol-Adj] = m.Sserr[wycol]
	for i = m.Ncol; i >= 2; i-- {
		total = total + m.D[i-Adj]*m.Rhs[wycol][i-Adj]*m.Rhs[wycol][i-Adj]
		m.Rss[wycol][(i-1)-Adj] = total
	}

	m.Rss_set = true
}

//--------------------------------------------------------------------------
//   Cov()
//
//     Calculate covariance matrix for regression coefficients for the
//     first nreq variables, from an orthogonal reduction produced from
//     AS75.1.
//
//     Auxiliary routine called: Inv()
//
//     Fills in covmat, sterr as the major intended effect.
//
//     nreq, as above, is the number of regression coefficients to be
//     accounted for, include a count of one for the constant term.
//--------------------------------------------------------------------------

func (m *MillerLSQ) Cov(nreq int, covmat []float64, sterr []float64, wycol int) (err error, Variance float64) {

	var pos, row, start, pos2, col, pos1, k int
	var total float64

	//     Check that dimension of array covmat is adequate.

	needs := nreq * (nreq + 1) / 2
	if len(covmat) < needs {
		return errors.New(fmt.Sprintf("Cov() error: covmat length(%v) is "+
			"inadquequate for nreq==%v. should be >= nreq*(nreq+1)/2 == %v", len(covmat), nreq, needs)), 0
	}

	//     Check for small or zero multipliers on the diagonal.

	for row = 1; row <= nreq; row++ {
		if math.Abs(m.D[row-Adj]) < m.Vsmall {
			return errors.New(fmt.Sprintf("Cov() error: small or zero multiplier "+
				"on the diagonal(row-Adj=%v)", row-Adj)), 0
		}
	}

	//     Calculate estimate of the residual variance.

	if m.Nobs > int64(nreq) {
		if !m.Rss_set {
			m.SS(wycol)
		}
		Variance = m.Rss[wycol][nreq-Adj] / float64(m.Nobs-int64(nreq))
	} else {
		return errors.New(fmt.Sprintf("Cov() error: m.Nobs(%v) <= nreq(%v)", m.Nobs, nreq)), 0
	}

	m.Dim_rinv = nreq * (nreq - 1) / 2
	m.Rinv = make([]float64, m.Dim_rinv)
	m.Inv(nreq, m.Rinv)

	pos = 1
	start = 1
	for row = 1; row <= nreq; row++ {
		pos2 = start
		for col = row; col <= nreq; col++ {
			pos1 = start + col - row
			if row == col {
				total = 1.0 / m.D[col-Adj]
			} else {
				total = m.Rinv[(pos1-1)-Adj] / m.D[col-Adj]
			}
			for k = col + 1; k <= nreq; k++ {
				total = total + m.Rinv[pos1-Adj]*m.Rinv[pos2-Adj]/m.D[k-Adj]
				pos1 = pos1 + 1
				pos2 = pos2 + 1
			}
			covmat[pos-Adj] = total * Variance
			if row == col {
				sterr[row-Adj] = math.Sqrt(covmat[pos-Adj])
			}
			pos = pos + 1
		}
		start = start + nreq - row
	}

	return nil, Variance
} // end Cov

//--------------------------------------------------------------------------
//  Inv()
//     Invert first nreq rows and columns of Cholesky factorization
//     produced by AS 75.1.
//
//--------------------------------------------------------------------------

func (m *MillerLSQ) Inv(nreq int, rinv []float64 /*out*/) {

	var pos, row, col, start, k, pos1, pos2 int
	var total float64

	//     Invert R ignoring row multipliers, from the bottom up.

	pos = nreq * (nreq - 1) / 2
	for row = nreq - 1; row >= 1; row-- {
		start = m.Row_ptr[row-Adj]
		for col = nreq; col >= row+1; col-- {
			pos1 = start
			pos2 = pos
			total = 0.0
			for k = row + 1; k <= col-1; k++ {
				pos2 = pos2 + nreq - k
				total = total - m.R[pos1]*m.Rinv[pos2-Adj]
				pos1 = pos1 + 1
			}
			m.Rinv[pos-Adj] = total - m.R[pos1]
			pos = pos - 1
		}
	}

} // end inv

func (m *MillerLSQ) Partial_corr(In int, cormat []float64 /*out*/, dimc int, ycorr []float64 /*out*/, wycol int) (ifault int) {

	// Fortran comments:
	//
	//     Replaces subroutines PCORR and COR of:
	//     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
	//
	//     Calculate partial correlations after the variables in rows
	//     1, 2, ..., IN have been forced into the regression.
	//     If IN = 1, and the first row of R represents a constant in the
	//     model, then the usual simple correlations are returned.
	//
	//     If IN = 0, the value returned in array CORMAT for the correlation
	//     of variables Xi & Xj is:
	//       sum ( Xi.Xj ) / Sqrt ( sum (Xi^2) . sum (Xj^2) )
	//
	//     On return, array CORMAT contains the upper triangle of the matrix of
	//     partial correlations stored by rows, excluding the 1's on the diagonal.
	//     e.g. if IN = 2, the consecutive elements returned are:
	//     (3,4) (3,5) ... (3,m.Ncol), (4,5) (4,6) ... (4,m.Ncol), etc.
	//     Array YCORR stores the partial correlations with the Y-variable
	//     starting with YCORR(IN+1) = partial correlation with the variable in
	//     position (IN+1).
	//
	//--------------------------------------------------------------------------

	var base_pos, pos, row, col, col1, col2, pos1, pos2, k int

	var sumxx, sumxy, sumyy, tmp, tmp1, tmp2 float64
	rmsOffset := In
	rms := make([]float64, m.Ncol-In)

	workOffset := In
	work := make([]float64, m.Ncol-In)

	//     Some checks.

	ifault = 0
	if In < 0 || In > m.Ncol-1 {
		ifault = ifault + 4
	}
	if dimc < (m.Ncol-In)*(m.Ncol-In-1)/2 {
		ifault = ifault + 8
	}
	if ifault != 0 {
		return ifault
	}

	//     Base position for calculating positions of elements in row (IN+1) of R.

	base_pos = In*m.Ncol - (In+1)*(In+2)/2 - Adj

	//     Calculate 1/RMS of elements in columns from IN to (m.Ncol-1).

	if m.D[In+1-Adj] > 0.0 {
		rms[0] = 1.0 / math.Sqrt(m.D[In+1-Adj])
	}

	for col = In + 2; col <= m.Ncol; col++ {
		pos = base_pos + col
		sumxx = m.D[col-Adj]
		for row = In + 1; row <= col-1; row++ {
			//fmt.Printf("row is %d, col is %d, In is %d, m.Ncol is %d, pos is %d\n", row, col, In, m.Ncol, pos)
			tmp = m.R[pos]
			sumxx = sumxx + m.D[row-Adj]*tmp*tmp
			pos = pos + m.Ncol - row - 1
		}
		if sumxx > 0.0 {
			rms[col-rmsOffset-Adj] = 1.0 / math.Sqrt(sumxx)
		} else {
			rms[col-rmsOffset-Adj] = 0.0
			ifault = -col
		}
	}

	//     Calculate 1/RMS for the Y-variable

	sumyy = m.Sserr[wycol]
	for row = In + 1; row <= m.Ncol; row++ {
		sumyy = sumyy + m.D[row-Adj]*m.Rhs[wycol][row-Adj]*m.Rhs[wycol][row-Adj]
	}
	if sumyy > 0.0 {
		sumyy = 1.0 / math.Sqrt(sumyy)
	}

	//     Calculate sums of cross-products.
	//     These are obtained by taking dot products of pairs of columns of R,
	//     but with the product for each row multiplied by the row multiplier
	//     in array D.

	pos = 1
	for col1 = In + 1; col1 <= m.Ncol; col1++ {
		sumxy = 0.0
		// work[col1+1:m.Ncol] = 0.0
		for i := range work {
			work[i] = 0.0
		}
		pos1 = base_pos + col1
		for row = In + 1; row <= col1-1; row++ {
			pos2 = pos1 + 1
			for col2 = col1 + 1; col2 <= m.Ncol; col2++ {
				tmp1 = m.R[pos1]
				tmp2 = m.R[pos2]
				work[col2-workOffset-Adj] = work[col2-workOffset-Adj] + m.D[row-Adj]*tmp1*tmp2
				pos2 = pos2 + 1
			}
			sumxy = sumxy + m.D[row-Adj]*m.R[pos1]*m.Rhs[wycol][row-Adj]
			pos1 = pos1 + m.Ncol - row - 1
		}

		//     Row COL1 has an implicit 1 as its first element (in column COL1)

		pos2 = pos1 + 1
		for col2 = col1 + 1; col2 <= m.Ncol; col2++ {
			work[col2-workOffset-Adj] = work[col2-workOffset-Adj] + m.D[col1-Adj]*m.R[pos2]
			pos2 = pos2 + 1
			//fmt.Printf("col1-rmsOffset-Adj = %d\n", col1-rmsOffset-Adj)
			tmp1 = rms[col1-rmsOffset-Adj]
			tmp2 = rms[col2-rmsOffset-Adj]
			cormat[pos-Adj] = work[col2-workOffset-Adj] * tmp1 * tmp2
			pos = pos + 1
		}
		sumxy = sumxy + m.D[col1-Adj]*m.Rhs[wycol][col1-Adj]
		ycorr[col1-Adj] = sumxy * rms[col1-rmsOffset-Adj] * sumyy
	}

	for k = 1; k <= In; k++ {
		ycorr[k-Adj] = 0.0
	}

	return 0
} // end partial_corr

// ---------------------------------------------------------------
// Vmove()
//     Move variable from position FROM to position TO in an
//     orthogonal reduction produced by AS75.1.
//
// ---------------------------------------------------------------

func (m *MillerLSQ) Vmove(from int, to int) error {

	var d1, d2, x, d1new, d2new, cbar, sbar, y float64
	var i, first, last, inc, m1, m2, mp1, col, pos, row int

	//     Check input parameters

	if from < 1 || from > m.Ncol {
		return errors.New(fmt.Sprintf("Vmove() error: from(%d) was out of [1, %d] range", from, m.Ncol))
	}
	if to < 1 || to > m.Ncol {
		return errors.New(fmt.Sprintf("Vmove() error: to(%d) was out of [1, %d] range", to, m.Ncol))
	}

	if from == to {
		return errors.New(fmt.Sprintf("Vmove() error: to(%d) is equal from(%d)", to, from))
	}

	if !m.Rss_set {
		for wycol := range m.Rhs {
			m.SS(wycol)
		}
	}

	if from < to {
		first = from
		last = to - 1
		inc = 1
	} else {
		first = from - 1
		last = to
		inc = -1
	}

	// jea: Apologies for the ugly gotos below. This is a tricky and complex
	// routine, and I didn't try to digest it nor to improve it.
	// To maintain correctness, it is translated
	// directly from the fortran structure.
	//
	for i = first; i != (last + inc); i += inc {

		//     Find addresses of first elements of R in rows i and (i+1).

		m1 = m.Row_ptr[i-Adj]
		m2 = m.Row_ptr[i+1-Adj]
		mp1 = i + 1
		d1 = m.D[i-Adj]
		d2 = m.D[mp1-Adj]

		//     Special cases.

		if d1 < m.Vsmall && d2 < m.Vsmall {
			goto forty40
		}
		x = m.R[m1]
		if math.Abs(x)*math.Sqrt(d1) < m.Tol[mp1-Adj] {
			x = 0.0
		}
		if d1 < m.Vsmall || math.Abs(x) < m.Vsmall {
			m.D[i-Adj] = d2
			m.D[mp1-Adj] = d1
			m.R[m1] = 0.0
			for col = i + 2; col <= m.Ncol; col++ {
				m1 = m1 + 1
				x = m.R[m1]
				m.R[m1] = m.R[m2]
				m.R[m2] = x
				m2 = m2 + 1
			}
			// generalized with wycol
			for wycol := range m.Rhs {
				x = m.Rhs[wycol][i-Adj]
				m.Rhs[wycol][i-Adj] = m.Rhs[wycol][mp1-Adj]
				m.Rhs[wycol][mp1-Adj] = x
			}
			goto forty40
		} else if d2 < m.Vsmall {
			m.D[i-Adj] = d1 * x * x
			m.R[m1] = 1.0 / x
			// r(m1+1:m1+m.Ncol-i-1) = r(m1+1:m1+m.Ncol-i-1) / x
			for k := m1 + 1; k <= m1+m.Ncol-i-1; k++ {
				m.R[k] = m.R[k] / x
			}
			// generalized with wycol
			for wycol := range m.Rhs {
				m.Rhs[wycol][i-Adj] = m.Rhs[wycol][i-Adj] / x
			}
			goto forty40
		}

		//     Planar rotation in regular case.

		d1new = d2 + d1*x*x
		cbar = d2 / d1new
		sbar = x * d1 / d1new
		d2new = d1 * cbar
		m.D[i-Adj] = d1new
		m.D[mp1-Adj] = d2new
		m.R[m1] = sbar
		for col = i + 2; col <= m.Ncol; col++ {
			m1 = m1 + 1
			y = m.R[m1]
			m.R[m1] = cbar*m.R[m2] + sbar*y
			m.R[m2] = y - x*m.R[m2]
			m2 = m2 + 1
		}
		// generalized with wycol
		for wycol := range m.Rhs {
			y = m.Rhs[wycol][i-Adj]
			m.Rhs[wycol][i-Adj] = cbar*m.Rhs[wycol][mp1-Adj] + sbar*y
			m.Rhs[wycol][mp1-Adj] = y - x*m.Rhs[wycol][mp1-Adj]
		}
		//     Swap columns i and (i+1) down to row (i-1).

	forty40:
		pos = i - Adj
		for row = 1; row <= i-1; row++ {
			x = m.R[pos]
			m.R[pos] = m.R[pos-1]
			m.R[pos-1] = x
			pos = pos + m.Ncol - row - 1
		}

		//     Adjust variable order (VORDER), the tolerances (TOL) and
		//     the vector of residual sums of squares (RSS).

		m1 = m.Vorder[i-Adj]
		m.Vorder[i-Adj] = m.Vorder[mp1-Adj]
		m.Vorder[mp1-Adj] = m1
		x = m.Tol[i-Adj]
		m.Tol[i-Adj] = m.Tol[mp1-Adj]
		m.Tol[mp1-Adj] = x
		// generalized with wycol
		for wycol := range m.Rhs {
			m.Rss[wycol][i-Adj] = m.Rss[wycol][mp1-Adj] + m.D[mp1-Adj]*m.Rhs[wycol][mp1-Adj]*m.Rhs[wycol][mp1-Adj]
		}
	} //

	return nil
} // end Vmove

//--------------------------------------------------------------------------
// Reorder()
//     Re-order the variables in an orthogonal reduction produced by
//     AS75.1 so that the N variables in LIST start at position POS1,
//     though will not necessarily be in the same order as in LIST.
//     Any variables in Vorder() before position POS1 are not moved.
//
//     Auxiliary routine called: Vmove()
//
//     on call into the Go-version here, pos1 is zero-based
//--------------------------------------------------------------------------

func (m *MillerLSQ) Reorder(list []int, pos1 int) error {
	if len(list) == 0 {
		return nil
	}
	n := len(list)
	pos1++ // internally, pos1 is still 1-based

	var next, i, l, j int

	//     Check N.

	if n < 1 || n > m.Ncol+1-pos1 {
		return errors.New(fmt.Sprintf("Reorder(): %d == len(list) > m.Ncol+1-pos1 == %d", n, m.Ncol+1-pos1))
	}

	//     Work through VORDER finding variables which are in LIST.

	next = pos1
	i = pos1

ten10:

	l = m.Vorder[i-Adj]
	for j = 1; j <= n; j++ {
		if l == list[j-Adj] {
			goto forty40
		}
	}

thirty30:
	i = i + 1
	if i <= m.Ncol {
		goto ten10
	}

	//     If this point is reached, one or more variables in LIST has not
	//     been found.

	return errors.New(fmt.Sprintf("Reorder() error: could not find one or more variables in list=%v", list))

	//     Variable L is in LIST; move it up to position NEXT if it is not
	//     already there.

forty40:
	if i > next {
		err := m.Vmove(i, next)

		if err != nil {
			panic(fmt.Sprintf("error on call to m.Vmove(i=%v, next=%v) call: %s", i, next, err))
		}
	}
	next = next + 1
	if next < n+pos1 {
		goto thirty30
	}

	return nil
} // end reordr

func (m *MillerLSQ) Hdiag(xrow []float64, nreq int) (hii float64, err error) {

	//     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
	//
	//                         -1           -1
	// The hat matrix H = x(X'X) x' = x(R'DR) x' = z'Dz

	//              -1
	// where z = x'R

	// Here we only calculate the diagonal element hii corresponding to one
	// row (xrow).   The variance of the i-th least-squares residual is (1 - hii).
	//--------------------------------------------------------------------------

	var col, row, pos int
	var total float64
	wk := make([]float64, m.Ncol)

	//     Some checks

	if nreq > m.Ncol {
		return 0, errors.New(fmt.Sprintf("Hdiag() error: nreq=%v was > m.Ncol=%v", nreq, m.Ncol))
	}

	//     The elements of xrow.inv(R).sqrt(D) are calculated and stored in WK.

	hii = 0.0

	for col = 1; col <= nreq; col++ {
		if math.Sqrt(m.D[col-Adj]) <= m.Tol[col-Adj] {
			wk[col-Adj] = 0.0
		} else {
			pos = col - 1 - Adj
			total = xrow[col-Adj]
			for row = 1; row <= col-1; row++ {
				total = total - wk[row-Adj]*m.R[pos]
				pos = pos + m.Ncol - row - 1
			}
			wk[col-Adj] = total
			//hii = hii + total*total/m.D[col]
			hii += total * total / m.D[col-Adj]

		}
	} // ! col = 1, nreq

	return hii, nil
} // end Hdiag

// ---------------------------------------------------------------------
//   Varprd()
//     Calculate the variance of x'b where b consists of the first nreq
//     least-squares regression coefficients.
//
// ---------------------------------------------------------------------

func (m *MillerLSQ) Varprd(x []float64, nreq int, wycol int) float64 {

	var ifault, row int
	wk := make([]float64, nreq)
	var Variance float64

	//     Check input parameter values

	var fn_val float64
	ifault = 0
	if nreq < 1 || nreq > m.Ncol {
		ifault = ifault + 4
	}
	if m.Nobs <= int64(nreq) {
		ifault = ifault + 8
	}
	if ifault != 0 {
		panic(fmt.Sprintf("Error in function Varprd: ifault = %v\n", ifault))
	}

	//     Calculate the residual variance estimate.

	Variance = m.Sserr[wycol] / float64(m.Nobs-int64(nreq))

	//     Variance of x'b = var.x'(inv R)(inv D)(inv R')x
	//     First call BKSUB2 to calculate (inv R')x by back-substitution.

	m.Bksub2(x, wk, nreq)
	for row = 1; row <= nreq; row++ {
		if m.D[row-Adj] > m.Tol[row-Adj] {
			fn_val = fn_val + wk[row-Adj]*wk[row-Adj]/m.D[row-Adj]
		}
	}

	fn_val = fn_val * Variance

	return fn_val
}

//--------------------------------------------------------------------------
// Bksub2()
//     Solve x = R'b for b given x, using only the first nreq rows and
//     columns of R, and only the first nreq elements of R.
//
//--------------------------------------------------------------------------

func (m *MillerLSQ) Bksub2(x []float64, b []float64 /*out*/, nreq int) {

	var pos, row, col int
	var temp float64

	//     Solve by back-substitution, starting from the top.

	for row = 1; row <= nreq; row++ {
		pos = row - 1 - Adj
		temp = x[row-Adj]
		for col = 1; col <= row-1; col++ {
			temp = temp - m.R[pos]*b[col-Adj]
			pos = pos + m.Ncol - col - 1
		}
		b[row-Adj] = temp
	}

} // end Bksub2

func Min(a int, b int) int {
	if a <= b {
		return a
	}
	return b
}

//--------------------------------------------------------------------------
// Printc()
//     Print (partial) correlations calculated using partial_corr to unit lout.
//     If yCorOnly, print correlations with the Y-variable only.
//
//--------------------------------------------------------------------------

func (m *MillerLSQ) Printc(In int, cormat []float64, dimc int, ycorr []float64,
	vname []string, yname string, yCorOnly bool, formatFunc func(float64) string) (string, error) {

	s := ""
	text := ""
	var nrows, j1, j2, i1, i2, row, upos, tpos, last int

	//     Check validity of arguments

	if In >= m.Ncol {
		return "", errors.New("Printc error: In was >= m.Ncol")
	}
	if m.Ncol <= 1 {
		return "", errors.New("Printc error: m.Ncol was <= 1")
	}
	nrows = m.Ncol - In
	if dimc <= nrows*(nrows-1)/2 {
		return "", errors.New("Printc error: dimc too small.")
	}

	if !yCorOnly {
		s += "Correlation matrix\n           "
		j1 = In + 1

		for {
			j2 = Min(j1+6, m.Ncol)

			i1 = j1 - In
			i2 = j2 - In
			for k := j1; k <= j2; k++ {
				s += vname[m.Vorder[k-Adj]] + " "
			}
			//     Print correlations for rows 1 to i2, columns i1 to i2.

			s += "\n"
			for row = 1; row <= i2; row++ {
				text = " " + vname[m.Vorder[row+In-Adj]]
				if i1 > row {
					upos = (row-1)*(nrows+nrows-row)/2 + (i1 - row)
					last = upos + i2 - i1
					// WRITE(text(12:74), `(7(F8.5, " "))`) cormat(upos:last)
					s += text + "  "
					for k := upos; k <= last; k++ {
						s += formatFunc(cormat[k-Adj]) + "  "
					}

				} else {
					upos = (row-1)*(nrows+nrows-row)/2 + 1
					tpos = 12 + 9*(row-i1-Adj)
					// text(tpos:tpos+8) = char1 // " 1.0 "

					text += strings.Repeat(" ", int(tpos)) + fmt.Sprintf("%-8.1f", 1.0)
					last = upos + i2 - row - 1
					s += text
					if row < i2 {
						//WRITE(text(tpos+9:74), `(6(F8.5, " "))`) cormat(upos:last)
						for k := upos; k <= last; k++ {
							//s += formatFunc(cormat[k-Adj])
							s += fmt.Sprintf("%8.5f ", cormat[k-Adj])
						}
					}
				} // end else
				//WRITE(lout, `(a)`) text
				s += "\n"
			} // end for

			//     Move onto the next block of columns.

			j1 = j2 + 1
			if j1 > m.Ncol {
				break
			}
		} // end for

	} // end if !yCorOnly

	//     Correlations with the Y-variable.

	s += "\n\nCorrelations with the dependent variable: " + yname + "\n"
	j1 = In + 1

	for {
		j2 = Min(j1+7, m.Ncol)
		for k := j1; k <= j2; k++ {
			s += " " + vname[m.Vorder[k-Adj]] + "     "
		}

		s += "\n"
		for k := j1; k <= j2; k++ {
			s += formatFunc(ycorr[k-Adj]) + "   "
		}
		j1 = j2 + 1
		if j1 > m.Ncol {
			break
		}
	}

	return s, nil
}

func PrintSliceDiff(a, b []float64) {
	if len(a) != len(b) {
		fmt.Printf("len(a) is %d, while len(b) is %d\n", len(a), len(b))
		return
	}
	if len(a) == 0 {
		return
	}
	for i := range a {
		if !fEquals(a[i], b[i]) {
			fmt.Printf("a[%d]==%v and b[%d]==%v differ by %v\n", i, a[i], i, b[i], a[i]-b[i])
		}
	}

}

func RegressionTableString(nreq int, vname []string, Vorder []int, beta []float64, sterr []float64, Tstat []float64, Rss []float64) string {

	s := ""

	s += fmt.Sprintf("Variable    Regn.coeff.   Std.error  t-value   Res.sum of sq.\n")
	for i := int(0); i < nreq; i++ {
		s += fmt.Sprintf("%v  %12.4g  %11.4g  %7.2f  %14.6g\n", vname[Vorder[i]], beta[i], sterr[i], Tstat[i], Rss[i])
	}

	return s
}

func RstyleRegressionTableString(nreq int, vname []string, Vorder []int, beta []float64, sterr []float64, Tstat []float64, Rss []float64) string {

	s := ""

	s += fmt.Sprintf("Variable      Regn.coeff.     Std.error         t-value       Res.sum of sq.\n")
	s += fmt.Sprintf("--------      -----------     ---------         -------       --------------\n")
	for i := int(0); i < nreq; i++ {
		//s += fmt.Sprintf("%v  %14.7f  %14.7f  %14.3f  %14.0f\n", vname[Vorder[i]], beta[i], sterr[i], Tstat[i], Rss[i])
		s += fmt.Sprintf("%v  %14.7e  %14.7e  %14.3e  %14.0e\n", vname[Vorder[i]], beta[i], sterr[i], Tstat[i], Rss[i])
	}

	return s
}

// generate 1, 2, ..., max
func Seq(max int) []int {
	r := make([]int, max)
	for i := range r {
		r[i] = int(i + 1)
	}
	return r
}

// return column j from m.R
func (m *MillerLSQ) RColumn(k int) []float64 {
	if k >= m.Ncol {
		panic(fmt.Sprintf("RColumn() error: k=%d is out of bounds; m.Ncol = %d", k, m.Ncol))
	}
	//fmt.Printf("RColumn request for k=%d\n", k)
	r := make([]float64, m.Ncol)
	nextr := k - 1
	var i int
	for i = 0; i < m.Ncol; i++ {
		//fmt.Printf("nextr is %d   m.Ncol is %d\n", nextr, m.Ncol)
		if i == k {
			r[i] = 1.0
			continue
		} else if i > k {
			r[i] = 0
			continue
		}

		r[i] = m.R[nextr]
		nextr = nextr + m.Ncol - i - 2
	}

	//fmt.Printf("column %d from %s\n    is %v\n", k, UpperTriToString(m.R, m.Ncol, "\n", StdFormat6dec), r)

	return r
}

// returns true if found nan
func NanToZero(s []float64) bool {
	foundNaN := false
	for i := range s {
		if math.IsNaN(s[i]) {
			s[i] = 0
			foundNaN = true
		}
	}
	return foundNaN
}
