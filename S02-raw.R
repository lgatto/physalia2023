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

sp |>
    filterMsLevel(1) |>
    spectraData() |>
    data.frame() |>
    ggplot(aes(x = rtime,
               y = totIonCurrent)) +
    geom_line()

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

plotSpectra(sp2807[1])

plotSpectra(sp2807[1], xlim = c(400, 1000))
abline(v = precursorMz(sp2807)[-1], col = "red")


plotSpectra(sp2807[1], xlim = c(400, 1000))

## Zoom in mz values 521.1 and 522.5 to reveal the isotopic envelope
## of that peak.

pdf("spectrum.pdf")
plotSpectra(sp2807[1], xlim = c(521.1, 522.5))
abline(v = precursorMz(sp2807)[-1], col = "red")
dev.off()

## The plotSpectra() function is used to plot all 10 MS2 spectra.

letters[3]
letters[-1]

sp2807[-1]

plotSpectra(sp2807[-1])


ms_2 <- sp2807[-1]

plotSpectra(ms_2[7], xlim = c(126, 132))

plotSpectra(ms_2[7], xlim = c(126, 132),
            labels = letters)


z <- ms_2[7]


plotSpectra(z, xlim = c(126, 132))

mzLabel <- function(z) {
    pd <- peaksData(z)[[1]]
    mylabels <- format(pd[, "mz"], digits = 4)
    mylabels[pd[, "intensity"] < 1e5] <- ""
    mylabels
}

plotSpectra(ms_2[7], xlim = c(126, 132),
            labels = mzLabel,
            labelCol = "red")

## Filter MS2 level spectra and find any 2 MS2 spectra that have
## matching precursor peaks based on the precursor m/z values.

sp2 <- filterMsLevel(sp, 2L)

anyDuplicated(precursorMz(sp2))

sp2i <- sp2[which(precursorMz(sp2) == precursorMz(sp2)[37])]

sp2i


## plotSpectraOverlay(x, y) ## Spectra of length 2
plotSpectraOverlay(sp2i, col = c("red", "steelblue"),
                   xlim = c(150, 600))

plotSpectraOverlay(sp2[c(32, 100)], col = c("red", "steelblue"),
                   xlim = c(150, 600))


## plotSpectraMirror(x) ## two Spectra of lengths 1
plotSpectraMirror(sp2i[1], sp2i[2])

plotSpectraMirror(sp2i[1], sp2[100],
                  xlim = c(100, 600))

BiocManager::install("RforMassSpectrometry/SpectraVis")

library(SpectraVis)

plotlySpectra(sp2i[1])

browseSpectra(sp)


plotSpectra(sp[2807], xlim = c(521.2, 522.5))
grid()

peaksData(sp[2807])[[1]] |>
    as_tibble() |>
    filter(mz > 521.2, mz < 521.4)


## profile mode

plotSpectra(sp[2807], xlim = c(521.2, 521.4))


## peak picking / centroiding

## centroided
sp[2807] |>
    pickPeaks() |>
    plotSpectra(xlim = c(521.2, 522.5))

peaksData(pickPeaks(sp[2807])) |>
    as_tibble() |>
    filter(mz > 521.2, mz < 521.4)

library(magrittr)

pickPeaks(sp[2807]) %>%
    filterIntensity(1e7) %>%
    plotSpectra(xlim = c(521.2, 522.5))


sp2 <- Spectra(f2)

msLevel(sp2)
centroided(sp2)

table(msLevel(sp2), centroided(sp2))

table(msLevel(sp), centroided(sp))

## A note on efficiency

table(msLevel(sp))
filterMsLevel(sp, 2)

length(sp)

length(filterMsLevel(sp, 2))

sp[msLevel(sp) == 2]

fsciex <- dir(system.file("sciex", package = "msdata"),
              full.names = TRUE)

sciex <- Spectra(fsciex)
sciex

range(intensity(sciex[1]))

sciex2 <- filterIntensity(sciex, intensity = c(50000, Inf))


range(intensity(sciex2[1]))

plotSpectra(sciex[1])

x11()

plotSpectra(sciex2[1])

## BiocParallel
