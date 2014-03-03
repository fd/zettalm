Project Zettalm
===============

build linear regression models on zettabytes of data - using Alan Miller's AS 274 online QR decomposition

* fits linear models to data *

This code does basic linear regression, otherwise know as least-squares fiting.

* but: handles really big data *

This particular solution to fitting linear models can handle infinite data,
and that's the point.

Again: there a tons of model fitting packages, but the point of interest here
is that this online QR decomposition algorithm can handle Zetta-Bytes of observations,
*actually unlimited* rows of data. It needs space only proportional to
O(p^2) for p variables. It only ever has to look at each row once. The original
Fortran90 did not handle more than one dependent variable, but the current
code does. Again, observations of x and y variables are not stored in memory.
Weighted observations are supported, including negative weights to delete
cases.

* is still quite fast *

You are certainly most likely to be I/O-bound. For CPU time, I benchmarked
the current code at 10x the performance of R's linear model solver.

Asymptodically, for p variables (both x and y variables), each additional
row of observations incurs O(p^2) time (see the Includ() routine).


* enhancements over the fortran code *

I've added the ability to handle multiple y-variables at once, as well
as online computation of means and standard deviations for all x and y variables.
Also I've fixed bugs that were present when using the downdating functionality
to remove rows from the set. Also I've added an extensive acceptance test-suite to
allow easy refactoring. Execute the tests, as usual in go, with "go test -v".

* online, but can be a moving window too *

Even though it only ever looks at any row once, rows are weighted, and
hence rows can be deleted, if you wish, from the data set considered
by adding the same row again but using a weight of -1 (instead of the default +1 weight)
on the second call to Includ().


Origins:

The structure and algorithm at the heart of this file is derived from
Alan Miller's AS 274 Algorithm and his Fortran90 source code for that publication.
AS 274 in turn builds upon W. Morven Gentleman's AG 75 algorithm.
http://www.jstor.org/stable/2347147

Additional information and the Fortran90 source is available at
the mirror, http://jblevins.org/mirror/amiller/lsq.f90.
The web page there indicates that Alan Miller's version of the code
has been placed in the public domain. Nonetheless, extensive credit is due to him
for this rigorous and numerically exacting algorithm. As a tribute, we
name the central struct the MillerLSQ struct. LSQ is short for least
squares, the method of regression implemented here.

