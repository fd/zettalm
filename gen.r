#!/bin/bash
exec R --vanilla -q --slave -e "source(file=pipe(\"tail -n +4 $0\"))" --args $@
#debug: exec R --vanilla --verbose -e "source(file=pipe(\"tail -n +6 $0\"))" --args $@
#
#gen.r : an R script to generate a large data set with known structure


N = 50

au = 5
as = 4
ad=rnorm(N, mean=au, sd=as)

bu = 6
bs = 3
bd=rnorm(N, mean=bu, sd=bs)

cu = 10
cs = 30
cd = rnorm(N, mean=cu, sd=cs)

du = 0
ds = 100
dd = rnorm(N, mean=du, sd=ds)

eu = ad + bd
es = 45
ed = eu + rnorm(N, mean=0, sd=es)

g1 = bd + rgamma(N, shape=100, rate=1)

g2 = ed + 50*rgamma(N, shape=1, rate=1) + rnorm(N, mean=0, sd=2000)

g3 = (g1 - g2) + rnorm(N, mean=0, sd=2)
w=is.infinite(g3)
g3[w]=0


df=data.frame(ad,bd,cd,dd,ed,g1,g2,g3)

p0=proc.time()
cat(file="bigger.dat", c("n", colnames(df)),"\n")
write.table(df, file="bigger.dat", quote=F, col.names=F, row.names=T, append=T)
elap = proc.time() - p0
print(paste("time to write out data: ", elap["elapsed"], " sec"))

p1=proc.time()
sm=summary(lm(g3 ~ ., data=df))
elap = proc.time() - p1
print(paste("time to compute regression: ", elap["elapsed"], " sec"))
            
sink(file="bigger.fit")
sm
sink()

sm
