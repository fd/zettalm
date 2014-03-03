package lsq

import (
	"fmt"
	"math"
)

// facility for tracking running stddev and mean for each variable
// reference: http://en.wikipedia.org/wiki/Standard_deviation
/* R code for one covariate:

W = W0 = 0
A = A0 = 0
Q = Q0 = 0


for (i in 1:length(w)) {

 # save prev values before updating them
 W0 = W
 A0 = A
 Q0 = Q

 W = W0 + w[i] # W == s0, the sum of all weights seen.
 A = A0 + w[i]*(b[i]-A0)/W
 Q = Q0 + w[i]*(b[i]-A0)*(b[i]-A)
}

WMEAN = A
WVAR = Q/W
WSD  = sqrt( WVAR * W/(W -1 ))
*/

type SdTracker struct {
	// length of W, A, and Q
	// don't reuse Ncol since we embed inside MillerLSQ
	Nc int

	// W is the sum of all weights seen.
	W []float64

	// A is the weighted mean
	A []float64

	// Q is the weighted numerator for the variance (the Quadratic term)
	Q []float64
}

func NewSdTracker(ncol int) *SdTracker {
	s := &SdTracker{}
	s.ZeroTracker(ncol)
	return s
}

func (s *SdTracker) ZeroTracker(ncol int) {
	s.Nc = ncol
	s.W = make([]float64, ncol)
	s.A = make([]float64, ncol)
	s.Q = make([]float64, ncol)
}

// return the weighted mean vector for the rows to date.
// returns by reference, use DeepCopy on it if you
// want to preserve it
func (s *SdTracker) Mean() []float64 {
	return s.A
}

// return the weighted standard-deviation vector for rows to date.
func (s *SdTracker) Sd() []float64 {
	var wvar float64
	WSD := make([]float64, s.Nc)
	for i := range s.W {
		wvar = s.Q[i] / s.W[i]
		WSD[i] = math.Sqrt(wvar * s.W[i] / (s.W[i] - 1))
	}
	return WSD
}

func (s *SdTracker) AddObs(x []float64, weight float64) {

	if len(x) != s.Nc {
		panic(fmt.Sprintf("len(x) == %v but s.Nc == %v", len(x), s.Nc))
	}

	// note previous values, used in the update
	var A0i float64

	for i := range s.W {
		// W is the sum of all weights seen.
		s.W[i] += weight

		// need to save the old value for Q update as well as A update
		A0i = s.A[i]

		// A is the weighted mean
		s.A[i] = A0i + weight*(x[i]-A0i)/s.W[i]

		// update the quadratic term.
		// Q/W gives the weighted variance, when needed.
		s.Q[i] += weight * (x[i] - A0i) * (x[i] - s.A[i])
	}
}
