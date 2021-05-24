## Hex lattice code.
## Copied from hex.model.fns.R
## [2021-05-24 Mon]
## sp: spacing between adjacent pts in regular lattice
## w: window, vector of length 4 (x1,x2,y1,y2)
## n: number of cells (ON : n1, OFF : n2)
## sd: noise (FRACTION OF SPACING) 
## v: displacement vector (fraction of spacing) (bivariate case)
## theta: rotation (bivariate case) (radians)

## FUNCTIONS:
## hex_1 (for univariate parameter fits)
## hex_2_v (for bivariate parameter fits, the Troy model)
## hex_2_theta (for bivariate parameter fits, the Ringach model)
 
## exclusion zones act between cells of the same type rather than different types

spacing <- function(w, n) {
  ## given the window (w) and number of cells (n), calculates the spacing 
  area <- (w[2]-w[1])*(w[4]-w[3])
  density <- 2/sqrt(3) ## density of a hexagonal grid with sp=1
  sqrt(density * area / n)
}

hex <- function(w, sp) {  
  ## makes hexagonal lattice larger than window (w) with spacing (sp)
  ## might have problems with this if w contains negative values
  w1 <- max(w)
  A <- matrix(c(1,sqrt(3),1,-sqrt(3)),2,2)
  H <- matrix(NA, (15*round(w1/sp)+1)^2, 2)
  k <- 1
  for (i in (-5*round(w1/sp)):(10*round(w1/sp))) { 
    for (j in (-5*round(w1/sp)):(10*round(w1/sp))) {
      H[k,] <- 1/2 * sp * A %*% c(i,j)
      k <- k+1
    }
  }
  return(H)
}


hex_1 <- function(w, n, sd, ex) {
  ## univariate fits (for univ.plot)
 
  sp <- spacing(w, n)
  L <- hex(w, sp) 

  rstart <- runif(2, 0, sp) # random start point
  r <- t(t(L) + rstart)

  ## ignore points more than 150um from the window
  r<-r[r[,1]<(w[2]+150),]; r<-r[r[,1]>(w[1]-150),]
  r<-r[r[,2]<(w[4]+150),]; r<-r[r[,2]>(w[3]-150),] 

  randorder <- sample(1:nrow(r)) 
  
  ## ADD NOISE WHILST IMPOSING EXCLUSION ZONES
  rej <- NULL

  for (i in randorder) { ## for each point in turn...
    accept <- FALSE
    j <- 0
    while((accept == FALSE) && (j <= 100)) {
      ## keep generating noise vectors until one is accepted
      noise <- rnorm(2, 0, sp*sd)
      rtrial <- r[i,1:2] + noise
      ## calculate distance from trial point to other points
      disti <- apply(r[-c(i,rej),1:2], 1, 
                       function(x) sqrt(sum((x - rtrial)^2)))
      samei <- min(disti)
      accept <- (samei > ex)
      j <- j + 1
    }
    if (j > 100) {  ## give up
      rej <- c(rej,i)   
    } else {  ## accept
      r[i,1:2] <- r[i,1:2] + noise
    }
  }

  if(is.null(rej)==FALSE)  r <- r[-rej,]  ## remove rejected points
  
  inc <- inout(r, poly_win(w)) ## only return points within window
  r <- r[inc,]
 
  sim.d <- list(x=r[,1], y=r[,2], n=length(r[,1]), w=w, spacing=sp)
  return(sim.d)  
}


hex_2_v <- function(w, n1, n2, sd, v, ex) {
  ## Bivariate (the Troy model with exclusion zones)
  ## density of mosaics the same
  n <- round(mean(c(n1, n2)))        
  sp <- spacing(w, n)
  L <- hex(w, sp) 

  rstart  <- c(runif(1, 0, sp), runif(1, 0, sqrt(3)*sp))  ## random start
  r1 <- t(t(L) + rstart) 
  r2 <- t(t(L) + rstart + v*sp) ## 2nd mosaic, displaced by v

    ## ignore points more than 150um from the window
  r1<-r1[r1[,1]<(w[2]+150),]; r1<-r1[r1[,1]>(w[1]-150),]
  r1<-r1[r1[,2]<(w[4]+150),]; r1<-r1[r1[,2]>(w[3]-150),] 
  r2<-r2[r2[,1]<(w[2]+150),]; r2<-r2[r2[,1]>(w[1]-150),]
  r2<-r2[r2[,2]<(w[4]+150),]; r2<-r2[r2[,2]>(w[3]-150),] 

  ## NOISE + EXCLUSION ZONES FOR 1st LATTICE
  randorder <- sample(1:nrow(r1)) 
  rej <- NULL

  for (i in randorder) { ## for each point in turn...
    accept <- FALSE
    j <- 0
    while((accept == FALSE) && (j <= 100)) {
      ## keep generating noise vectors until one is accepted
      noise <- rnorm(2, 0, sp*sd)
      rtrial <- r1[i,1:2] + noise
      ## calculate distance from trial point to other points
      disti <- apply(r1[-c(i,rej),1:2], 1, 
                       function(x) sqrt(sum((x - rtrial)^2)))
      samei <- min(disti)
      accept <- (samei > ex)
      j <- j + 1
    }
    if (j > 100) {  ## give up
      rej <- c(rej,i)   
    } else {  ## accept
      r1[i,1:2] <- r1[i,1:2] + noise
    }
  }

  if(is.null(rej)==FALSE)  r1 <- r1[-rej,]  ## remove rejected points


  ## NOISE + EXCLUSION ZONES FOR 2st LATTICE
  randorder <- sample(1:nrow(r2)) 
  rej <- NULL

  for (i in randorder) { ## for each point in turn...
    accept <- FALSE
    j <- 0
    while((accept == FALSE) && (j <= 100)) {
      ## keep generating noise vectors until one is accepted
      noise <- rnorm(2, 0, sp*sd)
      rtrial <- r2[i,1:2] + noise
      ## calculate distance from trial point to other points
      disti <- apply(r2[-c(i,rej),1:2], 1, 
                       function(x) sqrt(sum((x - rtrial)^2)))
      samei <- min(disti)
      accept <- (samei > ex)
      j <- j + 1
    }
    if (j > 100) {  ## give up
      rej <- c(rej,i)   
    } else {  ## accept
      r2[i,1:2] <- r2[i,1:2] + noise
    }
  }

  if(is.null(rej)==FALSE)  r2 <- r2[-rej,]  ## remove rejected points
  
  inc <- inout(r1, poly_win(w)) ## only return points within window
  r1 <- r1[inc,]
  inc <- inout(r2, poly_win(w)) ## only return points within window
  r2 <- r2[inc,]

  d <- list(x=c(r1[,1], r2[,1]), y=c(r1[,2], r2[,2]), on=r1, of=r2, 
            n1=length(r1[,1]), n2=length(r2[,1]), 
            spacing=sp, w=w)
  return(d)
}


hex_2_theta <- function(w, n1, n2, sd1, sd2, theta, ex1, ex2) {
  ## The Ringach model with exclusion zones (bivariate)
  sp1 <- spacing(w, n1)  ## spacing for ON mosaic, with n1 points
  sp2 <- spacing(w, n2) ## spacing for OFF mosaic, with n2 points

  ## parameter from paik and ringach 2011
  alpha <- (sp2 - sp1)/sp1  
  ## (1+alpha)sp1 = sp2  

  L <- hex(w, sp1) 

  ## rotation matrix (from paik and ringach 2011)
  R_theta <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
  
  ## random start point for 1st lattice
  rstart1 <- runif(2, 0, sp1)
  r1 <- t(t(L) + rstart1)

  ## random start point for 2nd lattice and
  ## change spacing from sp to sp1 and rotation
  rstart2  <- runif(2, 0, sp2)
  r2 <- t(t(L) + rstart2)
  r2 <- (1+alpha) * t(apply(r1, 1, function(l) R_theta %*% l))

  ## period of interference patterns
  S <- (1+alpha) / (sqrt(alpha^2 + 2*(1-cos(theta))*(1+alpha)))

  ## ignore points more than 150um from the window
  r1<-r1[r1[,1]<(w[2]+150),]; r1<-r1[r1[,1]>(w[1]-150),]
  r1<-r1[r1[,2]<(w[4]+150),]; r1<-r1[r1[,2]>(w[3]-150),] 
  r2<-r2[r2[,1]<(w[2]+150),]; r2<-r2[r2[,1]>(w[1]-150),]
  r2<-r2[r2[,2]<(w[4]+150),]; r2<-r2[r2[,2]>(w[3]-150),] 

  ## NOISE + EXCLUSION ZONES FOR 1st LATTICE
  randorder <- sample(1:nrow(r1)) 
  rej <- NULL

  for (i in randorder) { ## for each point in turn...
    accept <- FALSE
    j <- 0
    while((accept == FALSE) && (j <= 100)) {
      ## keep generating noise vectors until one is accepted
      noise <- rnorm(2, 0, sp1*sd1)
      rtrial <- r1[i,1:2] + noise
      ## calculate distance from trial point to other points
      disti <- apply(r1[-c(i,rej),1:2], 1, 
                       function(x) sqrt(sum((x - rtrial)^2)))
      samei <- min(disti)
      accept <- (samei > ex1)
      j <- j + 1
    }
    if (j > 100) {  ## give up
      rej <- c(rej,i)   
    } else {  ## accept
      r1[i,1:2] <- r1[i,1:2] + noise
    }
  }

  if(is.null(rej)==FALSE)  r1 <- r1[-rej,]  ## remove rejected points

  ## NOISE + EXCLUSION ZONES FOR 2st LATTICE
  randorder <- sample(1:nrow(r2)) 
  rej <- NULL

  for (i in randorder) { ## for each point in turn...
    accept <- FALSE
    j <- 0
    while((accept == FALSE) && (j <= 100)) {
      ## keep generating noise vectors until one is accepted
      noise <- rnorm(2, 0, sp2*sd2)
      rtrial <- r2[i,1:2] + noise
      ## calculate distance from trial point to other points
      disti <- apply(r2[-c(i,rej),1:2], 1, 
                       function(x) sqrt(sum((x - rtrial)^2)))
      samei <- min(disti)
      accept <- (samei > ex2)
      j <- j + 1
    }
    if (j > 100) {  ## give up
      rej <- c(rej,i)   
    } else {  ## accept
      r2[i,1:2] <- r2[i,1:2] + noise
    }
  }

  if(is.null(rej)==FALSE)  r2 <- r2[-rej,]  ## remove rejected points
  
  inc <- inout(r1, poly_win(w)) ## only return points within window
  r1 <- r1[inc,]
  inc <- inout(r2, poly_win(w))
  r2 <- r2[inc,]

  d <- list(x=c(r1[,1], r2[,1]), y=c(r1[,2], r2[,2]), on=r1, of=r2, 
            n1=length(r1[,1]), n2=length(r2[,1]), 
            spacing=c(sp1, sp2), w=w, S=S)
  return(d)
}




