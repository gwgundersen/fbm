# Process PBC data.
#
# This data was used as an example in T. R. Fleming and D. P. Harrington's
# book, Counting Processes and Survival Analysis.
#
# On page 188, they note two errors in the data:  the age of subject 253
# should be 54.4 years, and the protime value for subject 107 should be 10.7.
# These have not been fixed here, to retain comparability with analyses
# done with these errors present.

d.all <- read.table("data",head=T)  # All the data
d.rnd <- d.all[1:312,]              # Data for subjects in randomized trial

# Data actually used, with censoring indicated by negative time of death.
# Variables selected and their codings are:
#
#   1. drug      0=placebo 1=drug
#   2. age       in years
#   3. sex       0=male 1=female
#   4. ascites   0=no 1=yes
#   5. hepatom   0=no 1=yes
#   6. spiders   0=no 1=yes
#   7. edema     0=no 0.5=yes,but(maybe)treatable 1=yes,untreatable
#   8. bili      bilirubin in mg/dl
#   9. albumin   albumin in gm/dl
#  10. alkphos   alkaline phosphatase in U/liter
#  11. platelet  platelets per cubic ml / 1000
#  12. protime   prothrombin time in seconds
#  13. time      time to death in days, or -time to censoring in days
#  14. id        (just for possible reference)
#
# Missing values for platelet are replaced by the median value of 257.

d.use <- data.frame (list 
 (drug=2-d.rnd$drug, age=round(d.rnd$age/365.25,3), sex=d.rnd$sex, 
  ascites=d.rnd$ascites, hepatom=d.rnd$hepatom, spiders=d.rnd$spiders, 
  edema=d.rnd$edema,
  bili=d.rnd$bili, albumin=d.rnd$albumin, alkphos=d.rnd$alkphos, 
  platelet=d.rnd$platelet, protime=d.rnd$protime,
  time = d.rnd$fudays*(d.rnd$status==2) - d.rnd$fudays*(d.rnd$status!=2)),
  id = d.rnd$id)

d.use$platelet[is.na(d.use$platelet)] <- 257

write.table(d.use,"data.use",row.names=F,col.names=F)

# Plot log of hazard functions at a grid of levels of one covariate, for one 
# iteration.

hazard1 <- function(levels, run, iter, lo=0, xlab="time", ylab="log hazard",...)
{
  r <- paste(run,".log.use",sep="")
  cmd <- paste("net-eval",r,iter,"/ 0 15 30 /",levels[1],levels[2],levels[3],
               " > /tmp/svtmpv")
  # cat(cmd,"\n")
  system(cmd)
  v <- read.table("/tmp/svtmpv",head=F)$V3
  n <- levels[3]+1
  v <- matrix(v,n,31)

  plot(c(0,15.6),lo+c(min(v),max(v)),type="n",xlab=xlab,ylab=ylab,...)
  for (i in 1:n)
  { lines(seq(0,15,0.5),lo+v[i,])
    label <- as.character(round(lo+seq(levels[1],levels[2],length=n)[i],3))
    text(15.2,lo+v[i,31],label,adj=c(0,0.5))
  }
}

# Data sets with each variable changed, used to assess how important each 
# variable is.

d.vary.drug.0 <- d.use
d.vary.drug.0$drug[] <- 0

d.vary.drug.1 <- d.use
d.vary.drug.1$drug[] <- 1

d.vary.sex.0 <- d.use
d.vary.sex.0$sex[] <- 0

d.vary.sex.1 <- d.use
d.vary.sex.1$sex[] <- 1

d.vary.ascites.0 <- d.use
d.vary.ascites.0$ascites[] <- 0

d.vary.ascites.1 <- d.use
d.vary.ascites.1$ascites[] <- 1

d.vary.hepatom.0 <- d.use
d.vary.hepatom.0$hepatom[] <- 0

d.vary.hepatom.1 <- d.use
d.vary.hepatom.1$hepatom[] <- 1

d.vary.spiders.0 <- d.use
d.vary.spiders.0$spiders[] <- 0

d.vary.spiders.1 <- d.use
d.vary.spiders.1$spiders[] <- 1

d.vary.edema.0 <- d.use
d.vary.edema.0$edema[] <- 0

d.vary.edema.1 <- d.use
d.vary.edema.1$edema[] <- 1

d.vary.age.0 <- d.use
d.vary.age.0$age[] <- d.use$age-sqrt(var(d.use$age))/2

d.vary.age.1 <- d.use
d.vary.age.1$age[] <- d.use$age+sqrt(var(d.use$age))/2

d.vary.bili.0 <- d.use
d.vary.bili.0$bili[] <- exp(log(d.use$bili)-sqrt(var(log(d.use$bili)))/2)

d.vary.bili.1 <- d.use
d.vary.bili.1$bili[] <- exp(log(d.use$bili)+sqrt(var(log(d.use$bili)))/2)

d.vary.albumin.0 <- d.use
d.vary.albumin.0$albumin[] <- exp(log(d.use$albumin)-sqrt(var(log(d.use$albumin)))/2)

d.vary.albumin.1 <- d.use
d.vary.albumin.1$albumin[] <- exp(log(d.use$albumin)+sqrt(var(log(d.use$albumin)))/2)

d.vary.alkphos.0 <- d.use
d.vary.alkphos.0$alkphos[] <- exp(log(d.use$alkphos)-sqrt(var(log(d.use$alkphos)))/2)

d.vary.alkphos.1 <- d.use
d.vary.alkphos.1$alkphos[] <- exp(log(d.use$alkphos)+sqrt(var(log(d.use$alkphos)))/2)

d.vary.platelet.0 <- d.use
d.vary.platelet.0$platelet[] <- exp(log(d.use$platelet)-sqrt(var(log(d.use$platelet)))/2)

d.vary.platelet.1 <- d.use
d.vary.platelet.1$platelet[] <- exp(log(d.use$platelet)+sqrt(var(log(d.use$platelet)))/2)

d.vary.protime.0 <- d.use
d.vary.protime.0$protime[] <- exp(log(d.use$protime)-sqrt(var(log(d.use$protime)))/2)

d.vary.protime.1 <- d.use
d.vary.protime.1$protime[] <- exp(log(d.use$protime)+sqrt(var(log(d.use$protime)))/2)

# Function for producing a plot of the effect of a covariate.

cov.effect <- function (cov, run, iters, max=14,
  xlab=paste("median survival with low value for",cov), 
  ylab=paste("median survival with high value for",cov), 
  pch=20, ...)
{
  r <- paste(run,".log.use",sep="")

  p0 <- c()
  for (i in iters)
  { p0 <- c(p0,cov.pred(cov,r,i,0))
  }

  p1 <- c()
  for (i in iters)
  { p1 <- c(p1,cov.pred(cov,r,i,1))
  }

  mini <- min(p0,p1)

  p0[p0>=max] <- max
  p1[p1>=max] <- max
  p0 <- p0 + (p0==max)*runif(length(p0))
  p1 <- p1 + (p1==max)*runif(length(p1))

  plot(c(mini,max+1),c(mini,max+1),type="n",xlab=xlab,ylab=ylab,...)
  points(p0,p1,pch=pch,...)
  lines(c(mini,max),c(mini,max),...)
  lines(c(max,max),c(mini,max),...)
  lines(c(mini,max),c(max,max),...)
}

cov.pred <- function (cov, run, iter, n)
{
  d <- get(paste("d.vary.",cov,".",n,sep=""))
  write.table(d,"/tmp/svtmpd",row.names=F,col.names=F)
  system(paste("net-pred Db",run,iter," / /tmp/svtmpd > /tmp/svtmpp"))
  scan("/tmp/svtmpp",quiet=T)
}

cov.effect1 <- function (cov, run, iters, max=14,
  xlab=paste("median survival with low value for",cov), 
  ylab=paste("median survival with high value for",cov), 
  pch=20, ...)
{
  r <- paste(run,".log.use",sep="")

  p0 <- c()
  for (i in iters)
  { p0 <- c(p0,cov.pred1(cov,r,i,0))
  }

  p1 <- c()
  for (i in iters)
  { p1 <- c(p1,cov.pred1(cov,r,i,1))
  }

  mini <- min(p0,p1)

  p0[p0>=max] <- max
  p1[p1>=max] <- max
  p0 <- p0 + (p0==max)*runif(length(p0))
  p1 <- p1 + (p1==max)*runif(length(p1))

  plot(c(mini,max+1),c(mini,max+1),type="n",xlab=xlab,ylab=ylab,...)
  points(p0,p1,pch=pch,...)
  lines(c(mini,max),c(mini,max),...)
  lines(c(max,max),c(mini,max),...)
  lines(c(mini,max),c(max,max),...)
}

cov.pred1 <- function (cov, run, iter, n)
{
  d <- get(paste("d.vary.",cov,".",n,sep=""))
  write.table(d[,cov],"/tmp/svtmpd",row.names=F,col.names=F,quote=F)
  cmd <- paste("net-pred Db",run,iter," / /tmp/svtmpd > /tmp/svtmpp")
  # cat(cmd,"\n")
  system(cmd)
  scan("/tmp/svtmpp",quiet=T)
}
