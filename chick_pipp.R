library(sjedrp)
library(readxl)
library(sjevor)
library(sjedist)
library(sjedmin)
library(sjedrp)
library(parallel)



read_pts = function(fieldname="DN1", type="Blue") {
  datafile = file.path("data", "journal.pone.0008992.s004.XLS")
  dat = read_excel(datafile, sheet=fieldname)
  names = dat[1,]
  matching_column = which(names(dat) == type)
  realx = truncate(dat[-c(1),matching_column])
  realy = truncate(dat[-c(1),matching_column+1])
  pts = cbind(realx, realy)
}


try_field <- function(prefix, delta, sigma, kappa) {

  theta = c(delta, sigma, kappa)
  x.lut <- seq(from=1, to=50, by=0.5)
  y.lut <- hpar(x.lut, theta)
  par$x.lut <- x.lut
  par$y.lut <- y.lut
  
  filename = sprintf('%s_%.2f_%.2f_%.2f.pdf', prefix, delta, sigma, kappa)

  sim = make.univsim("chick 1", par)
  pdf(filename, width=11, height=8)
  par(mfcol=c(2,3), oma=c(0,0,1,0))

  ## Plot dataset and LUT in column one.
  plot(par$pts.1, asp=1, pch=19, col=expt.col)
  plot(x.lut, y.lut, type='l', xlab='distance (um)', ylab='h()')

  if (sim$okay) {
    all.plots = sim$dist.arr$get()
    plot.spat.array(all.plots$g1)
    plot.spat.array(all.plots$l1)
    plot.spat.array(all.plots$f1)
    plot.spat.ri(all.plots$ri)
    msg <- "P-values"
  } else {
    msg <- "NO CONVERGENCE"
  }
  mtext(outer=TRUE, side=3, paste(filename, msg))
  dev.off()
}


hpar <- function(d,theta) {
  ## Choice of h() suggested by Peter.
  delta<-theta[1]
  sigma<-theta[2]
  kappa<-theta[3]
  res <- (0*(d<delta)) + (d>=delta)*(1-exp(-((d-delta)/sigma)^kappa))
  if (any (is.nan(res)))
    res[ which(is.nan(res))] <- 0

  res
}

sjespatdists.univ <- function (pts1, w, note, plot=F, param=NULL) {
  ## GENERAL ANALYSIS ROUTINE, UNIVARIATE VERSION.

  if( length(w) != 4)
    stop(paste("w (", paste(w, collapse=' '), ") should be of length 4."))
  else {
    xmin=w[1]; xmax=w[2]; ymin=w[3]; ymax=w[4]; 
  }
  
  ht <- ymax - ymin
  wid <- xmax - xmin


  steps   <- param$steps
  datapoly <- spoints( c(xmin,ymin, xmin, ymax,   xmax, ymax,  xmax,ymin))
  ## Check to see that all points are within the polygon.
  outside <- pip(pts1, datapoly, out=T, bound=TRUE)
  if (dim(outside)[1]>0) {
    print("some points outside polygon\n")
    browser()
  }

  
  null.xylist <-  list(x=NULL, y=NULL)

  g1 <- null.xylist
  if (!is.null(param$distribs$g1))
    g1 <- list(x=steps, y=Ghat(pts1,steps))
  
  f1 <- null.xylist
  if (!is.null(param$distribs$f1))
    f1 <- list(x=steps, y=Fhat(pts1, gridpts(datapoly, dim(pts1)[1]),steps))
  
  l1 <- null.xylist
  if (!is.null(param$distribs$l1))
    l1 <- list(x=steps, y= sqrt(khat(pts1, datapoly, steps)/pi))

  ## Voronoi areas.
  vd1 <- null.xylist
  if (!is.null(param$distribs$vd1))
    vd1 <- sje.vorarea(pts1, w, param$vd1.breaks, need.v=TRUE)

  if (!is.null(param$distribs$ri)) {
    ri <- calc.ri(pts1, w)
  } else{
    ri <- NULL
  }


  ## Before returning results, remove Voronoi plot, don't think we
  ## neeed to return that.
  vd1$v <- NULL

  res <- list(
    note = note,
    pts1=pts1,
    w = w,
    g1=g1,
    f1=f1,
    l1=l1,
    vd1=vd1,
    ri=ri,
    param=param)
  class(res) <- "sjespatdistsuniv"

  res
}


make.univsim <- function(field, allpar) {
  ## Make a univariate simulation
  with(allpar, {
    n1 <- nrow(pts.1)
    dist.real <- sjespatdists.univ(pts.1, w, "note", param=upar)
    dist.arr <- new.dist.arr(dist.real, nreps)
  
    rej <- matrix(NA, nrow=nreps, ncol=2)
    okay <- TRUE; t1.sim <- NA
    for (i in 1:nreps) {
      print(i)
      sim <- pipp.lookup(w=w, n=n1, pts=NULL, h=y.lut, d=x.lut,
                         nsweeps=10, verbose=FALSE)
      
      if (sim$okay) {
        rej[i,] <- 0 #TODO
        simpts <- cbind(sim$x, sim$y)
        t1.sim <- simpts
        dist.sim <- sjespatdists.univ(t1.sim, w, "note", param=upar)
        dist.arr$set.row(dist.sim, i+1)
      } else {
        okay = FALSE
        break #quit the loop, no point in continuing
      }
    }
    
    ##psfile <- paste("eg_", field, '.ps', sep='')
    ##postscript(file=psfile)
    
    ##par(mfrow=c(2,3), oma=c(1,0,0,0))
    ##dist.arr$plot()
    
    ## label <-"sim params here"
    ## mtext(label, side=1, outer=T)
      
    ##dev.off()
    
    ##image.name <- sprintf("last_%s.Rda", field)
    ##save.image(image.name)
    
    res <- list(okay=okay, dist.arr=dist.arr, t1.sim=t1.sim, allpar=allpar,
                rej=rej)

  })

}


truncate = function(d) {
  v = unlist(d)
  ok = !is.na(v)
  as.numeric(v[ok])
}


######################################################################
## End of functions
######################################################################


nreps = 99

type = "blue"
field = "DN1"

## full field list.
## [1] "Summary"       "DN1"           "DN2"           "DN3"          
##  [5] "DN4"           "DN5"           "DN6"           "DN7"          
##  [9] "DT1"           "DT2"           "DT3"           "DT4"          
## [13] "DT5"           "DT6"           "DT7"           "VN1"          
## [17] "VN2"           "VN3"           "VN4"           "VN5"          
## [21] "VN6"           "VN7"           "VT1"           "VT2"          
## [25] "VT3"           "VT4"           "VT5"           "VT6"          
## [29] "VT7"           "E18"           "P0"            "P6"           
## [33] "P. pubescens"  "P. domesticus" "C. livia"     


types = c("Red", "Green",  "Blue", "Violet", "Double")
expt.col = "blue"
prefix = sprintf("%s_%s", field, type)


## apply(par$pts.1, 2, range)
par <- list()
par$pts.1 <- read_pts()
par$w <- c(0, 430, 0, 325)
par$upar <- list(steps=seq(from=0, to=30, by=0.5),
                 distribs=list(l1=1, f1=1, g1=1,ri=1))
par$nreps <- nreps

## try_field(prefix, 8, 3, 1)
## try_field(prefix, 20, 3, 1)


deltas = seq(from=5, to=10, by=2.5)
sigma = seq(from=1, to=10, by=2)
kappa = seq(from=1, to=6, by=3)

grid = as.matrix(expand.grid(deltas, sigma, kappa))

ngrid = nrow(grid)
run.one = function(i) { try_field(prefix, grid[i,1], grid[i,2], grid[i,3]) }

options(mc.cores=2)
getOption("mc.cores",1)
mclapply(1:ngrid, run.one)

