library(Spectra)

## binary - vendor-specific formats
## open formats: mzML, mzXML

## proteowizard msconvert: https://proteowizard.sourceforge.io/
## ThermoRawFileParser: https://github.com/compomics/ThermoRawFileParser

spd <- DataFrame(msLevel = c(1L, 2L),
                 rtime = c(5, 10))

spd$mz <- list(c(100, 110, 140),
               c(200, 210, 300, 400))

spd$intensity <- list(c(101, 711, 50),
                      c(10, 78, 837, 3931))

spd

sp <- Spectra(spd)

sp

sp[1]

sp[2]


sp[c(1, 1, 2)]

spectraData(sp)

spectraVariables(sp)

peaksData(sp)

peaksData(sp)[[1]]
peaksData(sp)[[2]]


sp <- Spectra(f)
sp

## Exercise:
## - how many spectra? (see length())

sp
length(sp)

## - how many MS1 and MS2 spectra? (see sp$msLevel or msLevel(sp))

sp$msLevel
table(msLevel(sp))


## - how many of each precursor charges? (see precursorCharge(sp))

sp$precursorCharge
table(precursorCharge(sp))

rtime(sp)

peaksData(sp[100])[[1]]
mz(sp[1:2])
intensity(sp[1:2])

## MS3TMT11

library(msdata)
(f2 <- proteomics(full.names = TRUE)[3])

sp2 <- Spectra(f2)

## - How many spectra are there in that file?

length(sp2)

## - How many MS levels, and how many spectra per MS level?

table(msLevel(sp2))

## - What is the index of the MS2 spectrum with the highest
##   precursor intensity?

max(precursorIntensity(sp2), na.rm = TRUE)

which.max(precursorIntensity(sp2))

scanIndex(sp2)[which.max(precursorIntensity(sp2))]

## - Plot one spectrum of each level.

which(msLevel(sp2) == 1)

i <- 355

plot(peaksData(sp2[i])[[1]],
     type = "l")

plotSpectra(sp2[355])

plotSpectra(sp2[355:358])

spectraData(sp2)[i, ]

which(msLevel(sp2) == 2)

i <- 204

plot(peaksData(sp2[i])[[1]],
     type = "h")

spectraData(sp2)[i, ]


which(msLevel(sp2) == 3)

i <- 206

plot(peaksData(sp2[i])[[1]],
     type = "h")


sp


print(object.size(sp), units = "Mb")


fsciex <- dir(system.file("sciex", package = "msdata"),
              full.names = TRUE)

sciex <- Spectra(fsciex)
sciex

sciex_mem <- setBackend(sciex, MsBackendDataFrame())

print(object.size(sciex), units = "Mb")
print(object.size(sciex_mem), units = "Mb")


## Visualisation exercise

sp

spectraVariables(sp)

## The *chromatogram* be created by extracting the totIonCurrent and
## rtime variables for all MS1 spectra. Annotate the spectrum of
## interest.
##
## The vertical line identifies one scan in particular at retention
## time 1800.68 seconds (the 2807th scan).

plot(rtime(sp)[msLevel(sp) == 1],
     tic(sp)[msLevel(sp) == 1],
     type = "l")

sp1 <- filterMsLevel(sp, 1)

plot(rtime(sp1), tic(sp1), type = "l")
abline(v = rtime(sp)[2807], col = "red")


scanIndex(sp)[2807]
rtime(sp)[2807]



which(scanIndex(sp1) == 2807)

abline(v = rtime(sp1)[976], col = "blue")

which(rtime(sp) > 1896 & rtime(sp) < 1905)

filterRt(sp, c(1895, 1905))

filterRt(sp, c(1895, 1905), 2)

sp2807 <- filterPrecursorScan(sp, 2807)

## Plot the MS1 spectrum of interest and highlight all the peaks that
## will be selected for MS2 analysis. (hint: plotSpectra())
