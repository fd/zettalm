Project Zettalm
===============

Go code to build linear regression models on zettabytes of data - using Alan Miller's AS 274 online QR decomposition

* fits linear models to data

This golang code does basic linear regression, otherwise know as least-squares fitting.

* but: handles really big data

This particular solution to fitting linear models can handle infinite data,
and that's the point.

There are lots of model fitting packages that can handle finite data; they handle data that fits in RAM, or fits on disk, or in the cloud. The point of interest here
is that this online QR decomposition algorithm can handle Zetta-Bytes of observations,
*actually unlimited* rows of data. It needs space only proportional to
O(p^2) for p variables. It only ever has to look at each row once. The original
Fortran90 did not handle more than one dependent variable, but the current
code does. Again, observations of x and y variables are not stored in memory.
Weighted observations are supported, including negative weights to delete
cases.

* is still quite fast

You are certainly most likely to be I/O-bound. For CPU time, I benchmarked
the current code at 10x the performance of R's linear model solver.

Asymptodically, for p variables (both x and y variables), each additional
row of observations incurs O(p^2) time at the time of addition (see the Includ() routine).
When you want to get the regression coefficients for the current model, 
it costs one O(p^2) time operation, a call to Regcf(), to extract and
return the betas. If you need standard errors (to compute p-values for the
betas), the Cov() routine will supply them, along with the parameter covariance matrix,
and needs O(p^3) time (as does its subroutine, Inv()). SingularCheck() 
and SS() are both O(p) time routines.

In summary, for a typical fit on N rows for p variables, where N > p,
the time complexity is O(N*p^2). If p < N, then cost is O(p^3) time due to the 
Cov() call. Overall, time complexity is O(N*p^2 + p^3), as is typical for
linear regression.  Note that the structure of the code means that you
can fit multiple dependent y-variables at once, and still only ever make one
pass through the data. Data need not fit in main memory, nor even ever be
stored all in one place. For big problems, the savings on space can make
a huge difference in what is viable to handle.


* enhancements over the fortran code

I've added the ability to handle multiple y-variables at once, as well
as online computation of means and standard deviations for all x and y variables.
Also I've fixed minor bugs that were present when using the downdating functionality
to remove rows from the set. 

Most importantly, I've added an extensive acceptance test-suite to
allow easy refactoring. Execute the tests, as usual in go, with "go test -v".

* online, but can be a moving window too

Even though it only ever looks at any row once, rows are weighted, and
hence rows can be deleted, if you wish, from the data set considered
by adding the same row again but using a weight of -1 (instead of the default +1 weight)
on the second call to Includ().


Origins:

The structure and algorithm at the heart of this code is derived from
Alan Miller's AS 274 Algorithm and his Fortran90 source code for that publication.
AS 274 in turn builds upon W. Morven Gentleman's AS 75 algorithm. Refer to:

http://www.jstor.org/stable/2347147

http://lib.stat.cmu.edu/apstat/274

http://lib.stat.cmu.edu/apstat/75

Additional information and the Fortran90 source is available at
http://jblevins.org/mirror/amiller and http://jblevins.org/mirror/amiller/lsq.f90
The web page http://jblevins.org/mirror/amiller/ indicates that Alan Miller's 
version of the code has been placed in the public domain. Extensive credit is due to
Dr. Miller for this rigorous and numerically exacting algorithm. As a tribute, we
name the central struct the MillerLSQ struct. LSQ is short for least
squares, the method of regression implemented here.

License: MIT

Copyright (c) 2014, Jason E. Aten, Ph.D. <j.e.aten@gmail.com>

