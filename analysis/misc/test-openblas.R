# https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/mgcv-parallel.html

## illustration of multi-threading with gam...
require(mgcv)
set.seed(9)
dat <- gamSim(1,n=1e6,dist="poisson",scale=.1)
k <- 12;bs <- "cr";ctrl <- gam.control(nthreads=10, trace = TRUE)

times <- rep(NA_real_, 5)
names(times) <- c('gam', 'gam w 10 threads', 'bam', 'bam w 10 threads in control',
                  'bam w 10 threads in nthreads')

times[1] <- system.time(b1<-gam(y~s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k),
                                family=poisson,data=dat,method="REML"))[3]

times[2] <- system.time(b2<-gam(y~s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k),
                                family=poisson,data=dat,method="REML",control=ctrl))[3]

times[3] <- system.time(b3<-bam(y~s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k),
                                family=poisson,data=dat,method="fREML",discrete=TRUE))[3]

times[4] <- system.time(b4<-bam(y~s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k),
                                family=poisson,data=dat,method="fREML",discrete=TRUE,control=ctrl))[3]

times[5] <- system.time(b5<-bam(y~s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k),
                                family=poisson,data=dat,method="fREML",discrete=TRUE,nthreads=2))[3]

# repeat with a larger dataset
library(parallel)
set.seed(9)
dat_large <- gamSim(1,n=1e8,dist="poisson",scale=.1)

times_large <- rep(NA_real_, 3)
names(times_large) <- c('bam', 'bam with cluster', 'bam w 10 threads in nthreads')
times_large

cl <- makeCluster(5) # 5-core cluster

times_large[1] <- system.time(b3<-bam(y~s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k),
                                      family=poisson,data=dat_large,method="fREML",
                                      discrete=TRUE))[3]
times_large[2] <- system.time(b4<-bam(y~s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k),
                                      family=poisson,data=dat_large,method="fREML",
                                      cluster = cl))[3]
times_large[3] <- system.time(b5<-bam(y~s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k),
                                      family=poisson,data=dat_large,method="fREML",
                                      discrete=TRUE,nthreads=10))[3]
times_large

stopCluster(cl)
