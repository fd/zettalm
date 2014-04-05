package lsq

import (
	"math"
	"testing"

	cv "github.com/smartystreets/goconvey/convey"
)

// simple validation of stddev.go

var NAN float64 = math.NaN()

func TestSdTracker(t *testing.T) {

	x := []float64{4, 5, 6, 7}
	w := []float64{1, 2, 3, 4}

	eps := 1e-7

	expectedMean := []float64{4, 4.66666666666667, 5.33333333333333, 6}
	expectedSd := []float64{NAN, 0.577350269189626, 0.816496580927726, 1.05409255338946}

	ncol := 2

	// weighted is different from unweighted. The equivalent unweighted vector is
	// a=c(4, 5,5, 6,6,6, 7,7,7,7) # equivalent without weights
	cv.Convey("Given a series of weighted observations of a single variable, b=c(4, 5, 6, 7) and weights w=c(1,2,3,4)", t, func() {
		cv.Convey("When we compute a running weighted-mean and weighted-sd", func() {
			tracker := NewSdTracker(ncol)
			cv.Convey("We should see a weighted mean sequence of {4, 4.667, 5.33, 6}  and a weighted sd sequence of {NaN, 0.577, 0.816, 1.054}", func() {
				obsVector := make([]float64, ncol)

				for i := range x {

					for j := 0; j < 2; j++ {
						obsVector[j] = x[i] + float64(j)
					}

					tracker.AddObs(obsVector, w[i])
					mean := tracker.Mean()
					sd := tracker.Sd()

					//fmt.Printf("mean = %v   sd = %v\n", mean, sd)
					for j := 0; j < 2; j++ {
						cv.So(mean[j], cv.ShouldAlmostEqual, expectedMean[i]+float64(j), eps)

						if i == 0 {
							// special handling for NaN, as it won't be equal to NaN
							cv.So(math.IsNaN(sd[j]), cv.ShouldEqual, true)
						} else {
							cv.So(sd[j], cv.ShouldAlmostEqual, expectedSd[i], eps)
						}
					}
				}
			})
		})
	})

}

func TestMultiVariateSdTracker(t *testing.T) {

	df, err := readData("fuelcons.dat")
	if err != nil {
		panic(err)
	}
	nvar := df.Ncol - 2
	m := NewMillerLSQ(nvar-1, 2)

	last := df.Ncol - 1
	for i := range df.Rows {
		m.Includ(1.0, df.Rows[i][1:last-1], df.Rows[i][last-1:], NAN_OMIT_ROW)
	}

	eps := 1e-3

	XexpectedSd := []float64{4441.1087021, 0.9507698, 2384.9691183, 0.5736238, 3.4915072, 2124.2705778, 5.5470265, 111.8858156}
	XexpectedMean := []float64{4296.916667, 7.668333, 2362.083333, 4.241833, 5.565417, 2253.041667, 57.033333, 576.770833}

	YexpectedSd := []float64{5.547026549972453, 111.88581557530145}
	YexpectedMean := []float64{57.033333333333324, 576.7708333333331}

	cv.Convey("On the fuelcons.dat data with 1.0 weights for all", t, func() {
		cv.Convey("When we Includ() data", func() {
			cv.Convey("We should see, for X, the expected Mean and Sd when asked for", func() {
				sd := m.XStats.Sd()
				mean := m.XStats.Mean()
				//fmt.Printf("mean is %v\n", mean)
				//fmt.Printf("sd is %v\n", sd)
				for j := 0; j < m.Nxvar; j++ {
					cv.So(sd[j], cv.ShouldAlmostEqual, XexpectedSd[j], eps)
					cv.So(mean[j], cv.ShouldAlmostEqual, XexpectedMean[j], eps)
				}

			})

			cv.Convey("We should see, for Y, the expected Mean and Sd when asked for", func() {
				ysd := m.YStats.Sd()
				ymean := m.YStats.Mean()
				//fmt.Printf("ymean is %v\n", ymean)
				//fmt.Printf("ysd is %v\n", ysd)
				for j := 0; j < m.Nyvar; j++ {
					cv.So(ysd[j], cv.ShouldAlmostEqual, YexpectedSd[j], eps)
					cv.So(ymean[j], cv.ShouldAlmostEqual, YexpectedMean[j], eps)
				}

			})
		})
	})

}
