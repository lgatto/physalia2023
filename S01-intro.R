## Setup

## Bioconductor - http://bioconductor.org/

install.packages("BiocManager")
BiocManager::install("Spectra")

BiocManager::install("rpx")



BiocManager::install("RforMassSpectrometry/SpectraVis")

## Mass spectrometry - separation

## 0. liquid chromatography
## 1. source - ioinise
## 2. analyser - separate ions based m/z
## 3. detector - measurements

## Getting data

library(rpx)

px <- PXDataset("PXD000001")

px

pxfiles(px)

f <- pxget(px, "PXD000001_mztab.txt")

f <- pxget(px, "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML")

pxurl(px)
pxref(px)
pxtax(px)

rpxCache()

library(msdata)
msdata::proteomics()

msdata::quant()

msdata::ident()
