package lsq

import (
	"fmt"
	"testing"
	"time"

	cv "github.com/smartystreets/goconvey/convey"
)

// simple validation of stddev.go

func TestLsqOnBiggerData(t *testing.T) {

	load0 := time.Now()
	df, err := readData("bigger.dat")
	if err != nil {
		panic(err)
	}
	//fmt.Printf("df.Ncol is \n%v\n", df.Ncol)
	//fmt.Printf("df.Rows[0] is \n%v\n", df.Rows[0])

	loadElapsed := time.Since(load0)
	fmt.Printf("load Elapsed = %v\n", loadElapsed)

	nxvar := df.Ncol - 2 // has a rowname column, we omit it and the y
	nyvar := 1
	m := NewMillerLSQ(nxvar, nyvar)

	time1 := time.Now()

	// add dataframe's rows
	last := df.Ncol - 1
	for i := range df.Rows {
		m.Includ(1.0, df.Rows[i][1:last], df.Rows[i][last:], NAN_OMIT_ROW)
		// formula: g3 ~ .    same as: g3 ~ ad + bd + cd + dd + ed + g1 + g2
	}

	computeElap := time.Since(time1)
	fmt.Printf("compute Elapsed = %v\n", computeElap)

	/*
		lindep := make([]bool, m.Ncol)
		nSing := m.SingularCheck(&lindep)

		cv.Convey("The SingularCheck() call should find full rank.", t, func() {
			cv.So(nSing, cv.ShouldEqual, 0)
			for i := range lindep {
				cv.So(lindep[i], cv.ShouldEqual, false)
			}
		})
	*/

	vname := []string{"Constant"}
	vname = append(vname, df.Colnames[1:]...)
	vname = NormalizeNameLengths(vname)
	//fmt.Printf("vname is %v\n", vname)

	//y_name := df.Colnames[8]
	//fmt.Printf("y_name is %v\n", y_name)

	/*
		// now re-order
		list := []int64{2, 4, 5, 7}
		var startPosForListVarsInNewOrder int64 = 1

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

		cv.Convey("After m.Reorder([]int64{2, 4, 5, 7}, 1) we should have m.Vorder = {0,2,4,5,7,1,3,6}.", t, func() {
			cv.So(m.Vorder[0], cv.ShouldEqual, 0)
			cv.So(m.Vorder[1], cv.ShouldEqual, 2)
			cv.So(m.Vorder[2], cv.ShouldEqual, 4)
			cv.So(m.Vorder[3], cv.ShouldEqual, 5)
			cv.So(m.Vorder[4], cv.ShouldEqual, 7)
			cv.So(m.Vorder[5], cv.ShouldEqual, 1)
			cv.So(m.Vorder[6], cv.ShouldEqual, 3)
			cv.So(m.Vorder[7], cv.ShouldEqual, 6)
		})
	*/

	//     Calculate regression coefficients of Y against the first variable, which
	//     was just a constant = 1.0, and the next 4 predictors.

	m.Tolset(1e-12) // Calculate tolerances before calling
	// subroutine regcf.
	nreq := 8 // i.e. [Constant ad bd cd     dd ed g1 g2]
	wycol := 0
	// use Seq(nreq-1) to get regcf for all [Constant ad bd cd     dd ed g1 g2], since the constant is auto-added in.
	err, beta := m.Regcf(Seq(nreq-1), wycol)
	if err != nil {
		panic(err)
	}
	//fmt.Printf("at nreq-1 == %d, TestLsqOnBiggerData beta is %v\n", nreq-1, beta)

	knownGoodBeta := []float64{2.629355171810687, -0.03126581441662779, 0.011638188334615052, 0.012780715540548826, 8.820750124372001e-05, 0.010290090247564304, 0.976696861779601, -1.0000597616721187}
	cv.Convey("m.Regcf() should compute the correct parameter estimates or betas.", t, func() {
		cv.So(EpsSliceEqual(beta, knownGoodBeta, 1e-10), cv.ShouldEqual, true)
		if !EpsSliceEqual(beta, knownGoodBeta, 1e-10) {
			fmt.Printf("beta = %v\n\n vs knownGoodBeta=%v\n", beta, knownGoodBeta)
		}
	})

	m.SS(wycol) // Calculate residual sums of squares

	//     Calculate covariance matrix of the regression coefficients & their
	//     standard errors.

	//fmt.Printf("m.Rss[wycol] is %v\n", m.Rss[wycol])

	//	Variance := m.Rss[nreq] / float64(m.Nobs-nreq)
	Variance := m.Rss[wycol][nreq-1] / float64(m.Nobs-int64(nreq)) // why nreq-1 now, when nreq alone worked before?
	//fmt.Printf("Variance is %v\n", Variance)
	//fmt.Printf("m.Nobs is %v\n", m.Nobs)
	//fmt.Printf("nreq is %v\n", nreq)

	knownGoodRss := []float64{2.0334211954093957e+08, 1.9950106202615145e+08, 1.8416195500808036e+08, 1.7824795305768508e+08, 1.7810270992356238e+08, 1.78099909015074e+08, 1.7531515084464008e+08, 193.63577076132108}
	cv.Convey("m.SS() the residual sum-of-squares function should compute m.Rss that is correct", t, func() {
		cv.So(m.Nobs, cv.ShouldEqual, 50)
		cv.So(EpsSliceEqual(m.Rss[wycol], knownGoodRss, 1e-10), cv.ShouldEqual, true)
		cv.So(Variance, cv.ShouldAlmostEqual, 4.610375494317169, 1e-6)
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

	knownGoodSterrB := []float64{3.1214069698785494, 0.08327898815677234, 0.11348499151696355, 0.012088586772025463, 0.0037003724190309527, 0.008337974266046808, 0.030200518573496846, 0.00016217526544635535}
	knownGoodCovmatB := []float64{9.743181471606388, -0.02331183300712017, 0.001764208726048025, 0.009116416768602666, 0.000339217583374947, 0.0016044443577463972, -0.09042463746941756, 7.880897339339675e-05, 0.006935389868415828, 0.0005009344354861528, -0.00033596734136479337, 1.560986099979306e-06, -7.084914693713952e-05, -9.130655599324626e-05, 6.415002331800996e-07, 0.01287884329960529, 0.00012775769112102395, 3.240895778916354e-05, -6.505450882666637e-05, -0.0008246425071051237, -4.1098316909040395e-06, 0.000146133930144789, -5.723521101170216e-06, 9.996883651552994e-06, -9.443504615206631e-05, 4.0183055490145883e-07, 1.3692756039524985e-05, -2.9957302205232422e-06, -4.872236127971948e-06, -1.3376527923551208e-08, 6.95218148612588e-05, -1.714012158690983e-05, 5.9386858895031014e-09, 0.000912071322108128, -6.377130846002467e-07, 2.6300816722595816e-08}

	cv.Convey("m.Cov() should compute the correct Covariance of parameter estimate values", t, func() {
		cv.So(EpsSliceEqual(covmat, knownGoodCovmatB, 1e-10), cv.ShouldEqual, true)
		cv.So(EpsSliceEqual(sterr, knownGoodSterrB, 1e-10), cv.ShouldEqual, true)
		cv.So(Var2, cv.ShouldAlmostEqual, 4.610375494317169, 1e-10)
	})

	//     Calculate t-values
	Tstat := make([]float64, nreq)
	for k := 0; k < nreq; k++ {
		Tstat[k] = beta[k] / sterr[k]
	}

	//fmt.Printf("Tstat is %v\n", Tstat)

	TstatKnownGood := []float64{0.8423621774359632, -0.37543460972136244, 0.10255266515022292, 1.0572547297360713, 0.023837465869670393, 1.2341235315952865, 32.34040036109591, -6166.53691868271}
	cv.Convey("The t-values obtained should be correct; match previous implementation", t, func() {
		cv.So(EpsSliceEqual(Tstat, TstatKnownGood, 1e-10), cv.ShouldEqual, true)
	})

	//     Output regression table, residual sums of squares, and R-squared.
	rtab := RstyleRegressionTableString(nreq, vname, m.Vorder, beta, sterr, Tstat, m.Rss[wycol])
	//	fmt.Printf("Variable    Regn.coeff.   Std.error  t-value   Res.sum of sq.\n")
	//	for i := int64(0); i < nreq; i++ {
	//		fmt.Printf("%v  %12.4g  %11.4g  %7.2f  %14.6g\n", vname[m.Vorder[i]], beta[i], sterr[i], Tstat[i], m.Rss[i])
	//	}
	fmt.Printf(rtab)

	cv.Convey("The regression table obtained should be correct; match previous implementation", t, func() {
		cv.So(rtab, cv.ShouldEqual,
			`Variable      Regn.coeff.     Std.error         t-value       Res.sum of sq.
--------      -----------     ---------         -------       --------------
Constant   2.6293552e+00   3.1214070e+00       8.424e-01           2e+08
ad        -3.1265814e-02   8.3278988e-02      -3.754e-01           2e+08
bd         1.1638188e-02   1.1348499e-01       1.026e-01           2e+08
cd         1.2780716e-02   1.2088587e-02       1.057e+00           2e+08
dd         8.8207501e-05   3.7003724e-03       2.384e-02           2e+08
ed         1.0290090e-02   8.3379743e-03       1.234e+00           2e+08
g1         9.7669686e-01   3.0200519e-02       3.234e+01           2e+08
g2        -1.0000598e+00   1.6217527e-04      -6.167e+03           2e+02
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

}

// g3 ~ ad + bd + cd + dd + ed
func TestBiggerDataSmallerModelG3(t *testing.T) {

	load0 := time.Now()
	df, err := readData("bigger.dat")
	if err != nil {
		panic(err)
	}
	//fmt.Printf("df.Ncol is \n%v\n", df.Ncol)
	//fmt.Printf("df.Rows[0] is \n%v\n", df.Rows[0])

	loadElapsed := time.Since(load0)
	fmt.Printf("load Elapsed = %v\n", loadElapsed)

	nxvar := df.Ncol - 4 // has a rowname column, we omit it and the g1,g2,g3 as well
	nyvar := 1
	m := NewMillerLSQ(nxvar, nyvar)

	time1 := time.Now()

	// add dataframe's rows
	last := df.Ncol - 1
	for i := range df.Rows {
		m.Includ(1.0, df.Rows[i][1:last-2], df.Rows[i][last:], NAN_OMIT_ROW)
		// formula: g3 ~ .    same as: g3 ~ ad + bd + cd + dd + ed + g1 + g2
	}

	computeElap := time.Since(time1)
	fmt.Printf("compute Elapsed = %v\n", computeElap)

	vname := []string{"Constant"}
	vname = append(vname, df.Colnames[1:]...)
	vname = NormalizeNameLengths(vname)
	//fmt.Printf("vname is %v\n", vname)

	//y_name := df.Colnames[8]
	//fmt.Printf("y_name is %v\n", y_name)

	//     Calculate regression coefficients of Y against the first variable, which
	//     was just a constant = 1.0, and the next 4 predictors.

	m.Tolset(1e-12) // Calculate tolerances before calling
	// subroutine regcf.
	nreq := 6 // i.e. [Constant ad bd cd     dd ed g1 g2]
	wycol := 0

	// regress g3 on [ad bd cd dd ed]
	err, beta := m.Regcf(Seq(5), wycol)
	if err != nil {
		panic(err)
	}
	//fmt.Printf("at nreq-1 == %d, BiggerDataSmallerModel beta is %v\n", nreq-1, beta)

	knownGoodBeta := []float64{701.863273934701, 22.394916792915925, -180.25041095321492, 13.093753729871201, -0.6434176898397235, -0.20503547472106387}
	cv.Convey("m.Regcf() should compute the correct parameter estimates or betas.", t, func() {
		cv.So(EpsSliceEqual(beta, knownGoodBeta, 1e-10), cv.ShouldEqual, true)
	})

	m.SS(wycol) // Calculate residual sums of squares

	//     Calculate covariance matrix of the regression coefficients & their
	//     standard errors.

	//fmt.Printf("m.Rss is %v\n", m.Rss)

	//	Variance := m.Rss[nreq] / float64(m.Nobs-nreq)
	//Variance := m.Rss[wycol][nreq-1] / float64(m.Nobs-int64(nreq)) // why nreq-1 now, when nreq alone worked before?
	//fmt.Printf("Variance is %v\n", Variance)
	//fmt.Printf("m.Nobs is %v\n", m.Nobs)
	//fmt.Printf("nreq is %v\n", nreq)

	/*
		knownGoodRss := []float64{588366.4791666667, 468543.35883752245, 434888.65355693694, 401405.2109552021, 189049.96822312832, 177344.10625162232, 171100.06850570423, 150444.92571767746}
		cv.Convey("m.SS() the residual sum-of-squares function should compute m.Rss that is correct", t, func() {
			cv.So(m.Nobs, cv.ShouldEqual, 48)
			cv.So(nreq, cv.ShouldEqual, 5)
			cv.So(EpsSliceEqual(m.Rss, knownGoodRss, 1e-10), cv.ShouldEqual, true)
			cv.So(Variance, cv.ShouldAlmostEqual, 4124.281540735403, 1e-6)
		})
		cv.Convey("m.SS() the residual sum-of-squares function should compute m.Rss that is correct", t, func() {
			//	cv.So(Variance, cv.ShouldAlmostEqual, 0, 1e-6)
		})
	*/
	covmat := make([]float64, nreq*(nreq+1)/2)
	sterr := make([]float64, nreq)
	err, Var2 := m.Cov(nreq, covmat, sterr, wycol)
	if err != nil {
		panic(err)
	}
	//fmt.Printf("covmat is %v\n", covmat)
	//fmt.Printf("sterr is %v\n", sterr)
	//fmt.Printf("Var2 is %v\n", Var2)

	knownGoodSterrB := []float64{821.6233550104597, 77.90788488299032, 99.53252441444793, 10.765044311971389, 3.4625534191360896, 7.79444252154114}
	knownGoodCovmatB := []float64{675064.9374986439, -28720.13794384189, -67750.14244572366, -393.75032033043857, -117.39156461269752, -80.086303400403, 6069.638526941272, 459.24833281729155, -309.8529270832877, 1.2714670807429937, -63.59060299270026, 9906.723416312674, 90.64219190471397, 21.915355377989396, -71.68317429926249, 115.88617903870754, -5.276547787786842, 7.287708894484075, 11.989276180371023, -2.7139641037399898, 60.753334221608604}

	cv.Convey("m.Cov() should compute the correct Covariance of parameter estimate values", t, func() {
		cv.So(EpsSliceEqual(covmat, knownGoodCovmatB, 1e-7), cv.ShouldEqual, true)
		cv.So(EpsSliceEqual(sterr, knownGoodSterrB, 1e-10), cv.ShouldEqual, true)
		cv.So(Var2, cv.ShouldAlmostEqual, 4.047725204888046e+06, 1e-7)
	})

	//     Calculate t-values
	Tstat := make([]float64, nreq)
	for k := 0; k < nreq; k++ {
		Tstat[k] = beta[k] / sterr[k]
	}

	//fmt.Printf("Tstat is %v\n", Tstat)

	TstatKnownGood := []float64{0.8542396825194507, 0.28745379015937605, -1.8109699519190547, 1.2163213964024417, -0.1858217367229115, -0.026305341806603464}
	cv.Convey("The t-values obtained should be correct; match previous implementation", t, func() {
		cv.So(EpsSliceEqual(Tstat, TstatKnownGood, 1e-10), cv.ShouldEqual, true)
	})

	//     Output regression table, residual sums of squares, and R-squared.
	rtab := RstyleRegressionTableString(nreq, vname, m.Vorder, beta, sterr, Tstat, m.Rss[wycol])
	//	fmt.Printf("Variable    Regn.coeff.   Std.error  t-value   Res.sum of sq.\n")
	//	for i := int64(0); i < nreq; i++ {
	//		fmt.Printf("%v  %12.4g  %11.4g  %7.2f  %14.6g\n", vname[m.Vorder[i]], beta[i], sterr[i], Tstat[i], m.Rss[i])
	//	}
	fmt.Printf(rtab)

	cv.Convey("The regression table obtained should be correct; match previous implementation", t, func() {
		cv.So(rtab, cv.ShouldEqual,
			`Variable      Regn.coeff.     Std.error         t-value       Res.sum of sq.
--------      -----------     ---------         -------       --------------
Constant   7.0186327e+02   8.2162336e+02       8.542e-01           2e+08
ad         2.2394917e+01   7.7907885e+01       2.875e-01           2e+08
bd        -1.8025041e+02   9.9532524e+01      -1.811e+00           2e+08
cd         1.3093754e+01   1.0765044e+01       1.216e+00           2e+08
dd        -6.4341769e-01   3.4625534e+00      -1.858e-01           2e+08
ed        -2.0503547e-01   7.7944425e+00      -2.631e-02           2e+08
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

}
