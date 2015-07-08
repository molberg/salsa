library(FITSio)

## A routine to look-up a variable in a FITS header
lookup <- function(hdr, keyword) {
    i <- which(hdr$key == keyword)
    if (length(i) == 0) {
        value <- NA
    } else {
        value <- hdr$value[i[1]]
    }
    value
}

## The serial numbers of our spectra
ids <- seq(1909,1937)

MHz <- 1.0e6              # define a MHz in Hz
c.kms <- 2.99792458e5     # the speed of light in km/s

## Create an empty data frame, one line per spectrum.
head <- data.frame(id=ids, obs.utc=as.POSIXlt(Sys.time()),
                   target=NA, LII=0, BII=0, f0=0, v.lsr=0, dt=0)
## Empty matrices for data and frequency vectors, one column per spectrum.
data <- matrix(0, nrow=256, ncol=nrow(head))
freq <- matrix(0, nrow=256, ncol=nrow(head))

## Loop over FITS files
for (i in seq(nrow(head))) {
    id <- head$id[i]
    fits <- paste("spectrum_", id, ".fits", sep="")
    f <- readFITS(file=paste("FITS", fits, sep="/"), hdu = 1, phdu = 1)
    even <- seq(2, length(f$hdr), by=2)
    odd <- even-1
    hdr <- data.frame(key=f$hdr[odd], value=f$hdr[even], stringsAsFactors=FALSE)
    freq[,i] <- (f$axDat$crval[1]+f$axDat$cdelt[1]*(seq(f$axDat$len[1])-f$axDat$crpix[1]))/MHz
    bzero <- as.double(lookup(hdr, "BZERO"))
    if (is.na(bzero)) bzero <- 0.0
    bscale <- as.double(lookup(hdr, "BSCALE"))
    if (is.na(bscale)) bscale <- 1.0
    d <- f$imDat*bscale + bzero
    data[,i] <- d
    target <- lookup(hdr, "OBJECT")
    epoch <- as.double(lookup(hdr, "EPOCH"))
    if (is.na(epoch)) epoch <- as.double(lookup(hdr, "EQUINOX"))
    onx <- as.double(lookup(hdr, "CRVAL2"))
    ony <- as.double(lookup(hdr, "CRVAL3"))
    dt <- as.double(lookup(hdr, "OBSTIME"))
    f0 <- as.double(lookup(hdr, "CRVAL1"))/MHz
    fr <- as.double(lookup(hdr, "RESTFREQ"))/MHz
    if (is.na(fr)) fr <- 1420.40575177
    df <- as.double(lookup(hdr, "CDELT1"))
    vs <- as.double(lookup(hdr, "VLSR"))
    if (is.na(vs)) vs <- as.double(lookup(hdr, "VELO-LSR"))
    tsys <- as.double(lookup(hdr, "TSYS"))
    nch <- length(d)
    dv <- -df/f0*c.kms/MHz
    v <- dv*(seq(f$axDat$len[1])-f$axDat$crpix[1])-vs

    date <- lookup(hdr, "DATE-OBS")
    tstamp <- as.POSIXlt(sub("T", " ", date))
    ## fill in data frame
    head$obs.utc[i] <- tstamp
    head$target[i] <- target
    head$LII[i] <- onx
    head$BII[i] <- ony
    head$f0[i] <- f0
    head$v.lsr[i] <- vs
    head$dt[i] <- dt
    plot(v, d, , type='s', xlab="velocity [km/s]", ylab="",
         lty=1, xlim=c(-250,250), ylim=c(-20, 100))
}
grid()

## order by galactic longitude
o <- order(head$LII)
sd <- list(head=head[o,], freq=freq[,o], data=data[,o])

nch <- nrow(sd$data)
ns <- nrow(sd$head)
f0 <- matrix(sd$head$f0, nrow=nch, ncol=ns, byrow=TRUE)
vs <- matrix(sd$head$v.lsr, nrow=nch, ncol=ns, byrow=TRUE)
dv <- (f0-freq)/f0*c.kms
vel <- dv-vs

## plot all of the spectra in one go, offset each spectrum by 10 K in
## order to get a stacked picture
dT <- matrix(seq(ns)*10, nrow=nch, ncol=ns, byrow=TRUE)
matplot(vel, sd$data+dT, type='s', xlab="velocity [km/s]", ylab="", lty=1)
grid()

