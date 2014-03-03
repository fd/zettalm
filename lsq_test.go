package lsq

import (
	"fmt"
	"math"
	"testing"

	cv "github.com/smartystreets/goconvey/convey"
)

func TestDeepCopy(t *testing.T) {
	a := make([]float64, 3)
	for i := range a {
		a[i] = float64(i)
	}
	b := DeepCopy(a)
	a[0] = 99

	//fmt.Printf("a = %v\n", a)
	//fmt.Printf("b = %v\n", b)
	cv.Convey("When we make a deep copy of an array, modifying the original should leave the copy untouched", t, func() {
		cv.So(b[0], cv.ShouldEqual, 0)
	})

}

func TestSimplestNewLsq(t *testing.T) {

	nvar := 2
	m := NewMillerLSQ(nvar, 1)

	totalcol := nvar + 1

	cv.Convey("Given the MillerLSQ library, when we ask for a NewMillerLSQ, m", t, func() {
		cv.Convey("Then m.Nobs should start at 0", func() {
			cv.So(m.Nobs, cv.ShouldEqual, 0)
		})
		cv.Convey("Then m.R should be upper triangular sized.", func() {
			rsize := totalcol * (totalcol - 1) / 2
			cv.So(len(m.R), cv.ShouldEqual, rsize)
			cv.So(m.R_dim, cv.ShouldEqual, rsize)
		})
		cv.Convey("Then m.Ncol should account for the constant column.", func() {
			cv.So(m.Ncol, cv.ShouldEqual, totalcol)
		})
		cv.Convey("Displaying R for 2 var + const with 0 decimals and no data should yield [ 1 0 0 ; 0 1 0 ; 0 0 1].", func() {

			rString := UpperTriToString(m.R, m.Ncol, ";", StdFormat0dec)
			cv.So(rString, cv.ShouldEqual, "[ 1 0 0 ; 0 1 0 ; 0 0 1]")
		})
		cv.Convey("Displaying R for 2 var + const with 1 decimal and no data should yield [ 1.0 0.0 0.0 ; 0.0 1.0 0.0 ; 0.0 0.0 1.0]", func() {
			rString := UpperTriToString(m.R, m.Ncol, ";", StdFormat1dec)
			cv.So(rString, cv.ShouldEqual, "[ 1.0 0.0 0.0 ; 0.0 1.0 0.0 ; 0.0 0.0 1.0]")
		})
		cv.Convey("Displaying R for 2 var + const with 2 decimals and no data should yield [ 1.00 0.00 0.00 ; 0.00 1.00 0.00 ; 0.00 0.00 1.00].", func() {
			rString := UpperTriToString(m.R, m.Ncol, ";", StdFormat2dec)
			cv.So(rString, cv.ShouldEqual, "[ 1.00 0.00 0.00 ; 0.00 1.00 0.00 ; 0.00 0.00 1.00]")
		})

	})

}

func dot(a []float64, b []float64) float64 {
	prod := 0.0
	for i := range a {
		prod += a[i] * b[i]
	}
	return prod
}

func TestInclude(t *testing.T) {

	nvar := 2
	m := NewMillerLSQ(nvar, 1)

	// setup simple data, with expected beta:
	b := []float64{1, 2, 3}

	r0 := []float64{1, 2, 3}
	y0 := []float64{dot(b, r0)}
	r0 = r0[1:] // now we omit the 1
	/*
		r1 := []float64{1, 0, 3}
		y1 := dot(b, r1)

		r2 := []float64{1, 2, 0}
		y2 := dot(b, r2)
	*/

	cv.Convey("Given an empty LSQ", t, func() {
		origRmatrix := UpperTriToString(m.R, m.Ncol, ";", StdFormat2dec)

		cv.Convey("When we add and subtract an observation with Includ()", func() {
			//m.Includ(1.0, r0, len(r0)-1)
			m.Includ(1.0, r0, y0)

			afterIn := UpperTriToString(m.R, m.Ncol, ";", StdFormat6dec)
			cv.Convey("After only an include we should get [ 1.000000 2.000000 3.000000 ; 0.000000 1.000000 0.000000 ; 0.000000 0.000000 1.000000] ", func() {
				//cv.So(afterIn, cv.ShouldEqual, "[ 1.000000 2.000000 3.000000 ; 0.000000 1.000000 0.000000 ; 0.000000 0.000000 1.000000]")
				cv.So(afterIn, cv.ShouldEqual, "[ 1.000000 2.000000 3.000000 ; 0.000000 1.000000 0.000000 ; 0.000000 0.000000 1.000000]")
			})
			//m.Includ(-1.0, r0, len(r0)-1)
			m.Includ(-1.0, r0, y0)
			cv.Convey("Then we should get back the orignal identity matrix in R", func() {
				afterInOut := UpperTriToString(m.R, m.Ncol, ";", StdFormat2dec)

				cv.So(afterInOut, cv.ShouldEqual, origRmatrix)
			})
		})
	})

	/*
		cv.Convey("Given a small simple 3 row data set", t, func() {
			cv.Convey("When we add all 3 rows with Includ() we should get ...", func() {
				a := NewMillerLSQ(2, true)
				a.Includ(1.0, r0, y0)
				a.Includ(1.0, r1, y1)
				a.Includ(1.0, r2, y2)
				after3 := a.RtoString(2)

				cv.So(after3, cv.ShouldEqual, "")
			})
		})
	*/

}

// on actual data
func TestOnFuelData(t *testing.T) {

	nvar := 7
	m := NewMillerLSQ(nvar, 1)

	df, err := readData("fuelcons.dat")
	if err != nil {
		panic(err)
	}
	//fmt.Printf("in fuelconst.dat, df.Rows[0] = %v\n", df.Rows[0])
	last := df.Ncol - 1
	for i := range df.Rows {
		m.Includ(1.0, df.Rows[i][1:last], df.Rows[i][last:])
	}

	mat := "\n" + UpperTriToString(m.R, m.Ncol, "\n", StdFormat6dec)

	emat := "\n" + UpperTriToString(m.R, m.Ncol, "\n", StdFormatE)
	//fmt.Printf("R is \n%v\n", mat)

	cv.Convey("Given an LSQ with the fuelcons.dat data loaded, after 48 Includ() calls but before any other calls", t, func() {

		cv.Convey("Then the R matrix should match the lsq.f90 version output", func() {
			cv.So(mat, cv.ShouldEqual, `
[ 1.000000 4296.916667 7.668333 2362.083333 4.241833 5.565417 2253.041667 57.033333 
 0.000000 1.000000 -0.000031 0.532401 0.000053 0.000508 0.463974 -0.000456 
 0.000000 0.000000 1.000000 -87.904537 0.044881 -1.603833 -215.545385 -2.036955 
 0.000000 0.000000 0.000000 1.000000 -0.000000 -0.000707 1.270032 0.007452 
 0.000000 0.000000 0.000000 0.000000 1.000000 -1.347978 -260.251651 3.888110 
 0.000000 0.000000 0.000000 0.000000 0.000000 1.000000 54.974483 0.547643 
 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 1.000000 -0.006200 
 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 1.000000]`)

			// check more precisely for those small entries like row 4, column 5
			cv.So(emat, cv.ShouldEqual, `
[ 1.000000e+00 4.296917e+03 7.668333e+00 2.362083e+03 4.241833e+00 5.565417e+00 2.253042e+03 5.703333e+01 
 0.000000e+00 1.000000e+00 -3.139515e-05 5.324008e-01 5.295497e-05 5.081340e-04 4.639739e-01 -4.564454e-04 
 0.000000e+00 0.000000e+00 1.000000e+00 -8.790454e+01 4.488092e-02 -1.603833e+00 -2.155454e+02 -2.036955e+00 
 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 -4.463739e-09 -7.071586e-04 1.270032e+00 7.451966e-03 
 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 -1.347978e+00 -2.602517e+02 3.888110e+00 
 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 5.497448e+01 5.476433e-01 
 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 -6.199548e-03 
 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00]`)
		})

		u := m.XStats.Mean()
		sd := m.XStats.Sd()

		//fmt.Printf("u = %v\n", u)
		//fmt.Printf("sd = %v\n", sd)

		kgMean := []float64{4296.916666666666, 7.668333333333333, 2362.0833333333335, 4.241833333333332, 5.565416666666667, 2253.0416666666665, 57.033333333333324}
		kgSd := []float64{4441.10870206372, 0.9507697516051801, 2384.9691182561382, 0.5736237677697418, 3.491507166078876, 2124.270577772204, 5.547026549972453}
		cv.Convey("The mean and stddev should be what we expect", func() {
			//fmt.Printf("u = %v\n", u)
			//fmt.Printf("kgMean = %v\n", kgMean)
			cv.So(EpsSliceEqual(u, kgMean, 1e-10), cv.ShouldEqual, true)
			cv.So(EpsSliceEqual(sd, kgSd, 1e-10), cv.ShouldEqual, true)
		})

	})
}

func TestOnFuelPartialCorrelation(t *testing.T) {

	nvar := 7
	m := NewMillerLSQ(nvar, 1)
	wycol := 0

	lindep := make([]bool, m.Ncol)
	nSing := m.SingularCheck(&lindep, wycol)

	cv.Convey("The SingularCheck() call shouldn't crash on a newborn R matrix.", t, func() {
		cv.So(nSing, cv.ShouldEqual, 8)
		for i := range lindep {
			cv.So(lindep[i], cv.ShouldEqual, true)
		}
	})

	df, err := readData("fuelcons.dat")
	if err != nil {
		panic(err)
	}

	// save a version of R after each addition, so we can back then down and
	// check that we recover the previous step.
	Rversion := make([][]float64, len(df.Rows)+1) // +1 for original identity state
	Rversion[0] = DeepCopy(m.R)

	//r0 := "\n" + UpperTriToString(Rversion[0], nvar, "\n", StdFormatE)
	//fmt.Printf("\n Original Rversion[0]: %s\n", r0)

	// add dataframe's rows
	last := df.Ncol - 1
	for i := range df.Rows {
		m.Includ(1.0, df.Rows[i][1:last], df.Rows[i][last:])
		Rversion[i+1] = DeepCopy(m.R)
	}

	var i, j, N int
	const SkipBelow = 11
	cv.Convey("Downdating should recover the previous R matrix", t, func() {
		cv.Convey("When we run each of the 48 row additions backwards (but not past SkipBelow), we should get each previous matrix in turn.", func() {
			N = len(df.Rows)
			fmt.Printf("\n")
			for i = 0; i < N-SkipBelow; i++ {
				j = N - i - 1

				//fmt.Printf("Downdating data row j is %v\n", j)
				xrow := df.Rows[j]
				//m.Includ(-1.0, slc, whichY) // weight of -1 means remove this row (downdate)
				m.Includ(-1.0, xrow[1:last], []float64{xrow[last]}) // weight of -1 means remove this row (downdate)

				// we're okay down to 9
				if !fSliceEqual(Rversion[j], m.R) {
					fmt.Printf("detected difference at revision %d\n", j)

					PrintSliceDiff(Rversion[j], m.R)
				}

				cv.So(fSliceEqual(Rversion[j], m.R), cv.ShouldEqual, true)
			}
		})
	})

	// startover
	m = NewMillerLSQ(nvar, 1)

	// put back in the data
	for i := range df.Rows {
		m.Includ(1.0, df.Rows[i][1:last], []float64{df.Rows[i][last]})
	}

	lindep = make([]bool, m.Ncol)
	nSing = m.SingularCheck(&lindep, wycol)

	cv.Convey("The SingularCheck() call should find full rank.", t, func() {
		cv.So(nSing, cv.ShouldEqual, 0)
		for i := range lindep {
			cv.So(lindep[i], cv.ShouldEqual, false)
		}
	})

	In := 1
	var maxvar int = 30
	var max_cdim int = maxvar * (maxvar + 1) / 2
	cormat := make([]float64, m.R_dim)
	ycorr := make([]float64, m.Ncol)
	// compute cormat and ycorr
	m.Partial_corr(In, cormat, max_cdim, ycorr, wycol)

	vname := []string{"Constant"}
	vname = append(vname, df.Colnames[1:]...)
	//fmt.Printf("vname is %v\n", vname)

	y_name := df.Colnames[8]
	//fmt.Printf("y_name is %v\n", y_name)

	pcorStr, err := m.Printc(In, cormat, max_cdim, ycorr, vname, y_name, false, StdFormat5dec)
	if err != nil {
		panic(err)
	}
	fmt.Printf("pcor: %s\n", pcorStr)

	emat := "\n" + UpperTriToString(cormat, nvar, "\n", StdFormatE)
	//fmt.Printf("cormat is \n%v\n", emat)

	strycorr := fmt.Sprintf("%v", ycorr)
	//fmt.Printf("strycorr is %v\n", strycorr)

	cv.Convey("Given an LSQ with the fuelcons.dat data loaded, after Partial_corr()", t, func() {
		cv.Convey("The correlation matrix should match the R and Fortran code output.", func() {
			cv.So(emat, cv.ShouldEqual, `
[ 1.000000e+00 -1.466488e-01 9.913965e-01 4.099879e-01 6.463336e-01 9.700076e-01 -3.654433e-01 
 0.000000e+00 1.000000e+00 -1.796767e-01 1.266516e-02 -5.221301e-01 -2.366484e-01 -2.880372e-01 
 0.000000e+00 0.000000e+00 1.000000e+00 4.039095e-01 6.480528e-01 9.876866e-01 -2.992830e-01 
 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 5.016279e-02 3.325855e-01 1.570701e-01 
 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 7.018139e-01 -6.412950e-02 
 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 -2.864036e-01 
 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00]`)
		})
		cv.Convey("And the correlation against y should match the known good values.", func() {
			cv.So(strycorr, cv.ShouldEqual, `[0 -0.46275857811153115 -0.4512802751869866 -0.4232900127494368 -0.24486207498269863 0.019041938819586268 -0.36293909434233707 0.6989654213627073]`)
		})
	})

}

func TestReorderBetasCovarReporting(t *testing.T) {

	nvar := 7
	m := NewMillerLSQ(nvar, 1)
	wycol := 0

	df, err := readData("fuelcons.dat")
	if err != nil {
		panic(err)
	}

	// add dataframe's rows
	last := df.Ncol - 1
	for i := range df.Rows {
		m.Includ(1.0, df.Rows[i][1:last], df.Rows[i][last:])
	}

	lindep := make([]bool, m.Ncol)
	nSing := m.SingularCheck(&lindep, wycol)

	cv.Convey("The SingularCheck() call should find full rank.", t, func() {
		cv.So(nSing, cv.ShouldEqual, 0)
		for i := range lindep {
			cv.So(lindep[i], cv.ShouldEqual, false)
		}
	})

	vname := []string{"Constant"}
	vname = append(vname, df.Colnames[1:]...)
	//fmt.Printf("vname is %v\n", vname)

	//y_name := df.Colnames[8]
	//fmt.Printf("y_name is %v\n", y_name)

	// now re-order
	list := []int{2, 4, 5, 7}
	var startPosForListVarsInNewOrder int = 1

	fmt.Printf("current order of variables:\n")
	for i := range m.Vorder {
		fmt.Printf("%02d: %s   (m.Vorder[i==%d] is %v)\n", i, vname[m.Vorder[i]], i, m.Vorder[i])
	}

	err = m.Reorder(list, startPosForListVarsInNewOrder)
	if err != nil {
		panic(err)
	}

	fmt.Printf("After re-order of variables with list = %v at position %d:\n", list, startPosForListVarsInNewOrder)
	for i := range m.Vorder {
		fmt.Printf("%02d: %s   (m.Vorder[i==%d] is %v)\n", i, vname[m.Vorder[i]], i, m.Vorder[i])
	}

	cv.Convey("After m.Reorder([]int{2, 4, 5, 7}, 1) we should have m.Vorder = {0,2,4,5,7,1,3,6}.", t, func() {
		cv.So(m.Vorder[0], cv.ShouldEqual, 0)
		cv.So(m.Vorder[1], cv.ShouldEqual, 2)
		cv.So(m.Vorder[2], cv.ShouldEqual, 4)
		cv.So(m.Vorder[3], cv.ShouldEqual, 5)
		cv.So(m.Vorder[4], cv.ShouldEqual, 7)
		cv.So(m.Vorder[5], cv.ShouldEqual, 1)
		cv.So(m.Vorder[6], cv.ShouldEqual, 3)
		cv.So(m.Vorder[7], cv.ShouldEqual, 6)
	})

	//     Calculate regression coefficients of Y against the first variable, which
	//     was just a constant = 1.0, and the next 4 predictors.

	m.Tolset(1e-12) // Calculate tolerances before calling
	// subroutine regcf.
	nreq := 5 // i.e. Const, TAX, INC, ROAD, DLIC
	beta := make([]float64, df.Ncol)

	err = m.Regcf(beta, []int{1, 2, 3, 4}, 0)
	if err != nil {
		panic(err)
	}
	//fmt.Printf("beta is %v\n", beta)
	//fmt.Printf("m.D is %v\n", m.D)

	knownGoodBeta := []float64{377.29114647367385, -34.790149164432734, -66.58875178694376, -2.4258888852597504, 13.364493572600963, 0, 0, 0, 0}
	//fmt.Printf("knownGoodBeta is %v\n", knownGoodBeta)
	cv.Convey("m.Regcf() should compute the correct parameter estimates or betas for number of required betas (nreq)", t, func() {
		cv.So(EpsSliceEqual(beta, knownGoodBeta, 1e-10), cv.ShouldEqual, true)
	})

	m.SS(wycol) // Calculate residual sums of squares

	//     Calculate covariance matrix of the regression coefficients & their
	//     standard errors.

	Variance := m.Rss[wycol][nreq] / float64(m.Nobs-int64(nreq))
	//fmt.Printf("Variance is %v\n", Variance)
	//fmt.Printf("m.Nobs is %v\n", m.Nobs)
	//fmt.Printf("nreq is %v\n", nreq)
	//fmt.Printf("m.Rss is %v\n", m.Rss)

	knownGoodRss := []float64{588366.4791666667, 468543.35883752245, 434888.65355693694, 401405.2109552021, 189049.96822312832, 177344.10625162232, 171100.06850570423, 150444.92571767746}
	cv.Convey("m.SS() the residual sum-of-squares function should compute m.Rss that is correct", t, func() {
		cv.So(m.Nobs, cv.ShouldEqual, 48)
		cv.So(nreq, cv.ShouldEqual, 5)
		cv.So(EpsSliceEqual(m.Rss[wycol], knownGoodRss, 1e-10), cv.ShouldEqual, true)
		cv.So(Variance, cv.ShouldAlmostEqual, 4124.281540735403, 1e-6)
	})
	cv.Convey("m.SS() the residual sum-of-squares function should compute m.Rss that is correct", t, func() {
		//	cv.So(Variance, cv.ShouldAlmostEqual, 0, 1e-6)
	})

	covmat := make([]float64, nreq*(nreq+1)/2)
	sterr := make([]float64, nreq)
	err, Var2 := m.Cov(nreq, covmat, sterr, wycol)
	if err != nil {
		panic(err)
	}
	//fmt.Printf("covmat is %v\n", covmat)
	//fmt.Printf("sterr is %v\n", sterr)
	//fmt.Printf("Var2 is %v\n", Var2)

	knownGoodCovmat := []float64{34425.53346398021, -1875.170149827859, -657.318143861681, -331.2736516911481, -268.6604290931253, 168.22613071035346, -25.999719993904026, 25.125345490157322, 9.741839273140174, 296.58860136193107, -6.790461863024202, -6.3751315470831775, 11.486502840863489, 1.8143943977412207, 3.6978575618930836}
	knownGoodSterr := []float64{185.5411907474462, 12.97020164493804, 17.221747918313376, 3.3891743597613107, 1.9229814252595065}
	cv.Convey("m.Cov() should compute the correct Covariance of parameter estimate values", t, func() {
		cv.So(EpsSliceEqual(covmat, knownGoodCovmat, 1e-10), cv.ShouldEqual, true)
		cv.So(EpsSliceEqual(sterr, knownGoodSterr, 1e-10), cv.ShouldEqual, true)
		cv.So(Var2, cv.ShouldAlmostEqual, 4396.510888909961, 1e-10)
	})

	//     Calculate t-values
	Tstat := make([]float64, nreq)
	for k := 0; k < nreq; k++ {
		Tstat[k] = beta[k] / sterr[k]
	}

	fmt.Printf("Tstat is %v\n", Tstat)

	TstatKnownGood := []float64{2.0334630006079495, -2.6823136691949965, -3.8665501378134897, -0.7157757694799151, 6.949881780994023}
	cv.Convey("The t-values obtained should be correct; match previous implementation", t, func() {
		cv.So(EpsSliceEqual(Tstat, TstatKnownGood, 1e-10), cv.ShouldEqual, true)
	})

	//     Output regression table, residual sums of squares, and R-squared.
	rtab := RegressionTableString(nreq, vname, m.Vorder, beta, sterr, Tstat, m.Rss[wycol])
	//	fmt.Printf("Variable    Regn.coeff.   Std.error  t-value   Res.sum of sq.\n")
	//	for i := int(0); i < nreq; i++ {
	//		fmt.Printf("%v  %12.4g  %11.4g  %7.2f  %14.6g\n", vname[m.Vorder[i]], beta[i], sterr[i], Tstat[i], m.Rss[i])
	//	}
	fmt.Printf(rtab)

	cv.Convey("The regression table obtained should be correct; match previous implementation", t, func() {
		cv.So(rtab, cv.ShouldEqual,
			`Variable    Regn.coeff.   Std.error  t-value   Res.sum of sq.
Constant         377.3        185.5     2.03          588366
Tax             -34.79        12.97    -2.68          468543
Income          -66.59        17.22    -3.87          434889
RoadMls         -2.426        3.389    -0.72          401405
DLic             13.36        1.923     6.95          189050
`)
	})

	//     Output correlations of the parameter estimates

	fmt.Printf("\n\nCovariances of parameter estimates\n")
	i2 := nreq - 1
	for i := 1; i <= nreq-1; i++ {
		i1 := i2 + 1
		i2 = i2 + nreq - i
		fmt.Printf("%s", vname[m.Vorder[i]])
		for emptyCells := 1; emptyCells < i; emptyCells++ {
			fmt.Printf("          ") // 10 space empty cell
		}
		for k := i1; k <= i2; k++ {
			fmt.Printf("%10.3f", covmat[k])
		}
		fmt.Printf("\n")
	}

	//     Now delete the variable with the smallest t-value by moving it
	//     to position 5 and then repeating the calculations for the constant
	//     and the next 3 variables.

	j := 1
	var ti float64
	tmin := math.Abs(Tstat[1])
	for i := 2; i <= nreq-1; i++ {
		ti = math.Abs(Tstat[i])
		if ti < tmin {
			j = i
			tmin = ti
		}
	}
	j = j + 1 // Add 1 as the t-array started at subscript 0

	//fmt.Printf("Before Vmove(), m.Rhs[0] = %v\n", m.Rhs[0])
	//fmt.Printf("Before Vmove(), m.Sserr[0] = %v\n", m.Sserr[0])
	fmt.Printf("Before Vmove(), m.Vorder = %v\n", m.Vorder)
	kgRhs0 := []float64{576.7708333333331, -53.10629787821635, -46.65322425492822, -8.983324335629549, 13.364493572600963, -0.0066757328486958915, -0.04616781805459721, 0.10359613036358159}

	cv.Convey("Before Vmove, m.Rhs[0] should be what we expect", t, func() {
		cv.So(EpsSliceEqual(m.Rhs[0], kgRhs0, 1e-10), cv.ShouldEqual, true)
	})
	cv.Convey("Before Vmove, m.Sserr[0] should be what we expect", t, func() {
		cv.So(m.Sserr[0], cv.ShouldAlmostEqual, 150444.92571767746, 1e-10)
	})
	kgVO := []int{0, 2, 4, 5, 7, 1, 3, 6}
	cv.Convey("Before Vmove, m.Vorder should be [0 2 4 5 7 1 3 6]", t, func() {
		cv.So(IntSliceEqual(m.Vorder, kgVO), cv.ShouldEqual, true)
	})

	fmt.Printf("\n\nRemoving variable in position %d\n\n", j)
	err = m.Vmove(j, nreq)

	fmt.Printf("After Vmove(), m.Vorder = %v\n", m.Vorder)
	kgVO2 := []int{0, 2, 4, 7, 5, 1, 3, 6}
	cv.Convey("After Vmove, m.Vorder should be [0 2 4 7 5 1 3 6]", t, func() {
		cv.So(IntSliceEqual(m.Vorder, kgVO2), cv.ShouldEqual, true)
	})

	//fmt.Printf("After Vmove(), m.Rhs[0] = %v\n", m.Rhs[0])
	expected_after_move_rhs := []float64{576.7708333333331, -53.10629787821635, -46.65322425492822, 13.747684110568663, -2.4258888852597504, -0.0066757328486958915, -0.04616781805459721, 0.10359613036358159}
	cv.Convey("After Vmove, m.Rhs[0] should be what we expect", t, func() {
		cv.So(EpsSliceEqual(m.Rhs[0], expected_after_move_rhs, 1e-10), cv.ShouldEqual, true)
	})
	cv.Convey("After Vmove, m.Sserr[0] should be the same.", t, func() {
		cv.So(m.Sserr[0], cv.ShouldAlmostEqual, 150444.92571767746, 1e-10)
	})

	if err != nil {
		panic(err)
	}
	nreq = nreq - 1
	for i := range beta {
		beta[i] = 0
	}
	err = m.Regcf(beta, []int{1, 2, 3}, 0)
	if err != nil {
		panic(err)
	}
	m.SS(wycol)
	err, _ = m.Cov(nreq, covmat, sterr, wycol)
	if err != nil {
		panic(err)
	}

	for k := 0; k < nreq; k++ {
		Tstat[k] = beta[k] / sterr[k]
	}

	//     Output regression table, residual sums of squares, and R-squared.

	rtab2 := RegressionTableString(nreq, vname, m.Vorder, beta, sterr, Tstat, m.Rss[wycol])
	fmt.Printf(rtab2)

	cv.Convey("The second (after Vmove) regression table obtained should be correct; match previous implementation", t, func() {
		cv.So(rtab2, cv.ShouldEqual,
			`Variable    Regn.coeff.   Std.error  t-value   Res.sum of sq.
Constant         307.3        156.8     1.96          588366
Tax             -29.48        10.58    -2.79          468543
Income          -68.02        17.01    -4.00          434889
DLic             13.75        1.837     7.49          191302
`)
	})

	Variance = m.Rss[wycol][nreq-1] / float64(m.Nobs-int64(nreq))
	stdev_res := math.Sqrt(Variance)

	// RSS(1) is the rss if only the constant
	// is fitted; RSS(nreq) is the rss after
	// fitting the requested NREQ variables.
	r2 := 1.0 - m.Rss[wycol][nreq-1]/m.Rss[wycol][0]

	fmt.Printf("\nR^2 = %8.4f     Std. devn. of residuals = %12.4g\n", r2, stdev_res)
	fmt.Printf("\nN.B. Some statistical packages wrongly call the standard deviation\n")
	fmt.Printf("     of the residuals the standard error of prediction\n")

	/*
		fmt.Printf("\n\n")

		fmt.Printf(" after reorder and re-computation of regression:\n")
		fmt.Printf("m.Vorder is %v\n", m.Vorder)
		fmt.Printf("beta is %v\n", beta)
		fmt.Printf("sterr is %v\n", sterr)
		fmt.Printf("Tstat is %v\n", Tstat)
		fmt.Printf("m.Rss is %v\n", m.Rss)
		fmt.Printf("covmat is %v\n", covmat)
	*/

	kgVorder := []int{0, 2, 4, 7, 5, 1, 3, 6}
	kgBeta := []float64{307.3278964966638, -29.483808565617203, -68.02286155799504, 13.747684110568663, 0, 0, 0, 0, 0}
	kgSterr := []float64{156.8306699209346, 10.58357586517238, 17.009750282390495, 1.8366953982594196, 1.9229814252595065}
	kgTstat := []float64{1.9596160409925023, -2.7858078348207678, -3.9990511576419996, 7.48501037439138, 6.949881780994023}
	kgRss := []float64{588366.4791666666, 468543.3588375224, 434888.6535569369, 191302.45441913296, 189049.96822312832, 177344.10625162232, 171100.06850570423, 150444.92571767746}
	kgCovmat := []float64{24595.859027849147, -1137.7969473696087, -843.7009945158211, -213.93515407751184, 112.01207809385929, -11.022861312174904, 5.709084787771336, 289.33160466928354, -5.2437488438627735, 3.373449985987328, -6.790461863024202, -6.3751315470831775, 11.486502840863489, 1.8143943977412207, 3.6978575618930836}
	cv.Convey("After Vmove(), we should obtain the correct (from parallel/previous impl) Vorder, beta, sterr, Tstat, Rss, covmat, and R^2, stdev_res", t, func() {
		cv.So(r2, cv.ShouldAlmostEqual, 0.6749, 1e-4)
		cv.So(stdev_res, cv.ShouldAlmostEqual, 65.94, 1e-2)

		cv.So(EpsSliceEqual(beta, kgBeta, 1e-10), cv.ShouldEqual, true)
		cv.So(EpsSliceEqual(sterr, kgSterr, 1e-10), cv.ShouldEqual, true)
		cv.So(EpsSliceEqual(Tstat, kgTstat, 1e-10), cv.ShouldEqual, true)
		cv.So(EpsSliceEqual(m.Rss[wycol], kgRss, 1e-10), cv.ShouldEqual, true)
		cv.So(EpsSliceEqual(covmat, kgCovmat, 1e-10), cv.ShouldEqual, true)
		cv.So(IntSliceEqual(m.Vorder, kgVorder), cv.ShouldEqual, true)
	})

	//     Calculate residuals, hii, standardized residuals & standard errors of
	//     prediction.

	fmt.Printf("\nState     Actual    Fitted  Residual  Std.resid.  SE(prediction)\n")
	xx := make([]float64, m.Nobs)
	var fitted float64
	resid := make([]float64, m.Nobs)
	std_resid := make([]float64, m.Nobs)
	y := make([]float64, m.Nobs)
	// fill y
	for i := range df.Rows {
		y[i] = df.Rows[i][8]
	}

	var std_err_pred float64
	vStdErrPred := make([]float64, m.Nobs)
	vFitted := make([]float64, m.Nobs)
	for i := int64(0); i < m.Nobs; i++ {
		xx[0] = 1.0
		for j = 1; j < nreq; j++ {
			xx[j] = df.Rows[i][m.Vorder[j]] // N.B. Regression coefficient j is for
			//      the variable vorder(j+1)
		}
		fitted = 0
		for k := 0; k < nreq; k++ {
			fitted += beta[k] * xx[k]
		}
		vFitted[i] = fitted
		resid[i] = y[i] - fitted
		hii, err := m.Hdiag(xx, nreq)
		if err != nil {
			panic(err)
		}
		std_resid[i] = resid[i] / math.Sqrt(Variance*(1.0-hii))

		std_err_pred = math.Sqrt(m.Varprd(xx, nreq, wycol))
		vStdErrPred[i] = std_err_pred

		fmt.Printf("%s  %10.1f%10.1f%10.1f  %9.2f       %10.0f\n", df.Rownames[i], y[i], fitted, resid[i], std_resid[i], std_err_pred)
	}
	/*
		fmt.Printf("\n\n")
		fmt.Printf("y is %v\n", y)
		fmt.Printf("df.Rownames is %v\n", df.Rownames)
		fmt.Printf("vFitted is %v\n", vFitted)
		fmt.Printf("resid is %v\n", resid)
		fmt.Printf("std_resid is %v\n", std_resid)
		fmt.Printf("vStdErrPred is %v\n", vStdErrPred)
	*/

	eY := []float64{541, 524, 561, 414, 410, 457, 344, 467, 464, 498, 580, 471, 525, 508, 566, 635, 603, 714, 865, 640, 649, 540, 464, 547, 460, 566, 577, 631, 574, 534, 571, 554, 577, 628, 487, 644, 640, 704, 648, 968, 587, 699, 632, 591, 782, 510, 610, 524}
	eRn := []string{"ME", "NH", "VT", "MA", "RI", "CN", "NY", "NJ", "PA", "OH", "IN", "IL", "MI", "WI", "MN", "IA", "MO", "ND", "SD", "NE", "KS", "DE", "MD", "VA", "WV", "NC", "SC", "GA", "FL", "KY", "TN", "AL", "MS", "AR", "LA", "OK", "TX", "MT", "ID", "WY", "CO", "NM", "AZ", "UT", "NV", "WN", "OR", "CA"}
	evFitted := []float64{520.8173965873635, 549.9916010353209, 576.4309378974406, 482.1804859161812, 520.0988755930413, 434.104447111153, 329.6643807313973, 483.0191709398906, 496.2122520724045, 552.89424809106, 501.3963007307091, 459.26755971310706, 562.3921803591224, 564.0178419888504, 642.1253941906837, 612.8328092092445, 601.2046119489437, 590.4071792354257, 775.47775103501, 692.1464959511596, 699.9836899321745, 560.1100922844705, 411.37232440666594, 463.08754340856603, 503.0763494149816, 536.7345671637447, 590.287690578922, 620.5743167044114, 560.5722991918587, 474.78412158680385, 565.4680573936981, 579.4772338367183, 657.7195446104561, 609.8449068524513, 500.9849885798136, 719.4308470461744, 662.8752993246742, 641.4704339251605, 720.9238784563079, 729.2262752980691, 658.9125507874218, 626.2442701063292, 637.4282837052551, 544.57797281954, 699.5301943082306, 522.4960537859938, 665.1957433726243, 575.9285507809739}
	eResid := []float64{20.182603412636468, -25.991601035320855, -15.430937897440572, -68.18048591618123, -110.09887559304127, 22.89555288884702, 14.335619268602727, -16.019170939890614, -32.21225207240451, -54.894248091060035, 78.60369926929093, 11.732440286892938, -37.39218035912245, -56.017841988850364, -76.12539419068366, 22.167190790755512, 1.7953880510563067, 123.5928207645743, 89.52224896499001, -52.146495951159636, -50.98368993217446, -20.11009228447051, 52.62767559333406, 83.91245659143397, -43.07634941498162, 29.265432836255286, -13.287690578922025, 10.425683295588556, 13.427700808141253, 59.21587841319615, 5.531942606301868, -25.477233836718256, -80.71954461045607, 18.155093147548655, -13.984988579813603, -75.43084704617445, -22.875299324674188, 62.52956607483952, -72.92387845630788, 238.77372470193086, -71.91255078742176, 72.75572989367083, -5.428283705255126, 46.42202718045996, 82.46980569176935, -12.496053785993809, -55.19574337262429, -51.928550780973865}
	eStd_resid := []float64{0.32152880929372823, -0.4086818061632476, -0.2444510035746188, -1.0720516550262569, -1.6946063645714848, 0.3953645470140967, 0.24779127440840662, -0.25311079446791235, -0.49803520379933763, -0.8519649591836896, 1.213803447333082, 0.18776143950659802, -0.5834316632699668, -0.8677611013888814, -1.1763311039026982, 0.34166919031273957, 0.027675632971584586, 1.9306726707883226, 1.5062833365341226, -0.8567414812983033, -0.8079956786806023, -0.31521613061753384, 0.8444188065909397, 1.319664989897229, -0.6684402483091518, 0.4629688664173679, -0.2083597245282872, 0.16076273921754222, 0.2060921272018356, 0.9495834988380654, 0.08733737237568681, -0.40892873597448504, -1.3060407014345283, 0.2860013739760621, -0.22246343513682965, -1.191349839444032, -0.38961843340488694, 0.9676372026029857, -1.2111322707408247, 3.8027889481653117, -1.117411762506211, 1.1343737783804966, -0.08379756800838746, 0.7342803121158413, 1.37748875251569, -0.19638773403686438, -0.8562744710539871, -0.8162993862670482}
	evStdErrPred := []float64{17.90417605131339, 15.436481824195315, 16.896705473971906, 15.438070078070709, 9.98037309134902, 27.96205982291986, 28.05367756625469, 16.406256225112948, 11.372340983833181, 12.422676811148234, 11.010738305829234, 18.671354480949322, 13.745347693131732, 11.914453528085053, 11.211971290060657, 10.435898007813755, 10.468035360349436, 14.0162765191352, 25.325344920523193, 22.488934047064845, 16.972633759073545, 14.7760343812713, 19.091769889698163, 15.476746557622429, 12.379512857840762, 16.636987557723447, 14.860504871228946, 10.570610662099979, 8.989409335211507, 18.999861733240795, 16.25147381915012, 19.14753881027215, 20.376347704552977, 15.818842848646396, 17.64440138321412, 16.326315147578896, 26.614184497334136, 11.627816504898288, 23.834782098450038, 17.85347268649503, 12.729364113413695, 13.570828196441742, 10.916158728126854, 16.611281321468695, 24.502250610867364, 15.33608659949479, 12.308648095879507, 15.38466116724236}
	cv.Convey("The y, fitted, resid, std_resid and std_err_pred should match previous implementation", t, func() {
		cv.So(EpsSliceEqual(y, eY, 1e-10), cv.ShouldEqual, true)
		cv.So(StringSliceEqual(df.Rownames, eRn), cv.ShouldEqual, true)

		cv.So(EpsSliceEqual(vFitted, evFitted, 1e-10), cv.ShouldEqual, true)
		cv.So(EpsSliceEqual(resid, eResid, 1e-10), cv.ShouldEqual, true)
		cv.So(EpsSliceEqual(std_resid, eStd_resid, 1e-10), cv.ShouldEqual, true)
		cv.So(EpsSliceEqual(vStdErrPred, evStdErrPred, 1e-10), cv.ShouldEqual, true)

	})

}
