library(rpx)
library(pheatmap)
library(tidyverse)
library(rpx)
library(Spectra)
library(msdata)
library(PSMatch)

idf <- msdata::ident(full.names = TRUE)

basename(idf)

id <- PSM(idf)
id


dim(id)

names(id)

id$spectrumID


px <- PXDataset("PXD000001")
f <- pxget(px, "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML")

sp <- Spectra(f)
sp
spectraVariables(sp)

sp$spectrumId

## scan
## match/PSM
## peptide

## Verify that this table contains 5802 matches for 5343 scans and
## 4938 peptides sequences.

as.data.frame(id)

nrow(id)

length(unique(id$spectrumID))

length(unique(id$sequence))

table(table(id$spectrumID))

table(id$isDecoy)

(i <- which(id$spectrumID == "controllerType=0 controllerNumber=1 scan=1774"))

data.frame(id[i, ])[, c(1:5, 15, 30)]

id2 <- reducePSMs(id, id$spectrumID)


(j <- which(id2$spectrumID == "controllerType=0 controllerNumber=1 scan=1774"))

data.frame(id[i, ])[, c(1, 15)]

id2[j, "DatabaseAccess"]

id_tbl <- as_tibble(id)

## Remove decoy hits

id_tbl <- id_tbl |>
    filter(!isDecoy)

## Keep first rank matches

id_tbl <- id_tbl |>
    filter(rank == 1)

## Remove shared peptides. Start by identifying scans that match
## different proteins. For example scan 4884 matches proteins
## XXX_ECA3406 and ECA3415. Scan 4099 match XXX_ECA4416_1,
## XXX_ECA4416_2 and XXX_ECA4416_3. Then remove the scans that match
## any of these proteins.

mltm <- id_tbl |>
    group_by(spectrumID) |>
    mutate(nProts = length(unique(DatabaseAccess))) |>
    select(spectrumID, nProts) |>
    filter(nProts > 1)

id_tbl <- id_tbl |>
    filter(!spectrumID %in% mltm$spectrumID)

id_filtered <- filterPSMs(id)

describeProteins(id_filtered)

describeProteins(filterPsmDecoy(id))

describePeptides(id_filtered)

describePeptides(filterPsmDecoy(id))

id_filtered2 <-
    filterPsmDecoy(id) |>
    filterPsmRank()

describeProteins(id_filtered2)

describePeptides(id_filtered2)

am <- makeAdjacencyMatrix(id_filtered2)
dim(am)

cc <- ConnectedComponents(am)

connectedComponents(cc, 1)

connectedComponents(cc, 527)

connectedComponents(cc, 38)

connectedComponents(cc, 920)

i <- which(nrows(cc) > 2 & ncols(cc) > 2)

dims(cc)[i, ]

cx <- connectedComponents(cc, 1082)

plotAdjacencyMatrix(cx)

## Starting from the _non-filtered_ identification data (that still
## contains the decoy hits), compare the distribution of raw
## identification scores of the decoy and non-decoy hits. Interpret
## the figure.

as_tibble(id) |>
    ggplot(aes(x = MS.GF.RawScore,
               colour = isDecoy)) +
    geom_density()

id_filtered$MS.GF.QValue

id_filtered

sp

names(id_filtered)

spectraVariables(sp)

table(table(id_filtered$spectrumID))

which(table(id_filtered$spectrumID) == 4)

id_filtered[id_filtered$spectrumID == "controllerType=0 controllerNumber=1 scan=5490", ] |>
    as.data.frame() |>
    DT::datatable()

id_filtered <- reducePSMs(id_filtered, id_filtered$spectrumID)

sp <- joinSpectraData(sp, id_filtered,
                      by.x = "spectrumId",
                      by.y = "spectrumID")

spectraVariables(sp)

## Verify that the identification data has been added to the correct
## spectra: no identifications for MS1, only to (some) MS2 scans.

all(is.na(filterMsLevel(sp, 1)$sequence))

sp2 <- filterMsLevel(sp, 2)

table(!is.na(sp2$sequence))

sp2id <- sp2[!is.na(sp2$sequence)]

sp2id


sp <- countIdentifications(sp)

spectraVariables(sp)

table(msLevel(sp), sp$countIdentifications)

sp |>
    filterMsLevel(1) |>
    spectraData() |>
    data.frame() |>
    ggplot(aes(x = rtime,
               y = totIonCurrent)) +
    geom_line()

sp |>
    filterMsLevel(1) |>
    spectraData() |>
    as_tibble() |>
    ggplot(aes(x = rtime,
               y = totIonCurrent)) +
    geom_line(alpha = 0.25) +
    geom_point(aes(colour =
                       ifelse(countIdentifications == 0,
                              NA,
                              countIdentifications)),
               size = 2,
               alpha = 0.5) +
    labs(colour = "Nb ids")


plotSpectra(sp2[234])

## Find an MS2 scan that has an identication score > 100

i <- which(sp$MS.GF.RawScore > 100)[1]

plotSpectra(sp[i])

sp[i]$sequence

nchar(sp[i]$sequence)

calculateFragments(sp[i]$sequence)

sp[i]

addFragments(sp[i])

mz(sp[i])[[1]]


plotSpectra(sp[i], labels = addFragments,
            labelPos = 3, labelCol = "red",
            main = sp[i]$sequence)


## Create a new Spectra object containing the MS2 spectra with
## sequences "SQILQQAGTSVLSQANQVPQTVLSLLR" and
## "TKGLNVMQNLLTAHPDVQAVFAQNDEMALGALR".

k <- which(sp$sequence %in% c("SQILQQAGTSVLSQANQVPQTVLSLLR", "TKGLNVMQNLLTAHPDVQAVFAQNDEMALGALR"))

sp_k <- sp[k]

mat <- compareSpectra(sp_k)

colnames(mat) <- rownames(mat) <- strtrim(sp_k$sequence, 2)

pheatmap(mat)

## Compare the spectra with the plotting function seen previously.

sp_k

plotSpectra(sp_k)


plotSpectra(filterIntensity(sp_k, 1e4),
            labels = addFragments, labelCol = "red")


plotSpectraMirror(sp_k[1], sp_k[2], labels = addFragments,
                  labelCol = "red", labelPos = 3)

plotSpectraMirror(sp_k[3], sp_k[4], labels = addFragments,
                  labelCol = "red", labelPos = 3)

plotSpectraMirror(sp_k[1], sp_k[4], labels = addFragments,
                  labelCol = "red", labelPos = 3)

k2 <- sample(length(sp2id), 10)

mat2 <- compareSpectra(sp2id[k2])

sp2id[k2]$sequence

spectraVariables(sp2id)

summary(sp2id$experimentalMassToCharge - sp2id$calculatedMassToCharge)

## Download the 3 first mzML and mzID files from the _PXD022816_
## project from Morgenstern, David, Rotem Barzilay, and Yishai
## Levin. 2021. “RawBeans: A Simple, Vendor-Independent, Raw-Data
## Quality-Control Tool.” Journal of Proteome
## Research. (https://doi.org/10.1021/acs.jproteome.0c00956.).

px2 <- PXDataset("PXD022816")

pxfiles(px2)

fmz <- pxget(px2, grep("mzML$", pxfiles(px2), value = TRUE)[1:3])
fmz

pxget(px2, grep("mzML$", pxfiles(px2), value = TRUE)[4])

fid <- pxget(px2, grep("mzID.gz$", pxfiles(px2), value = TRUE)[1:3])

## Generate a Spectra object and a table of filtered PSMs.

sp <- Spectra(fmz)

sp

table(basename(sp$dataOrigin))

table(msLevel(sp), centroided(sp))

filterMsLevel(sp, 1) |>
    spectraData() |>
    as_tibble() |>
    ggplot(aes(x = rtime,
               y = totIonCurrent,
               colour = basename(dataOrigin))) +
    geom_line()

id <- PSM(fid)

id_filtered <- filterPSMs(id)

table(basename(id_filtered$spectrumFile))

## Visualise the total ion chromatograms and check the quality of the
## identification data by comparing the density of the decoy and
## target PSMs id scores for each file.


as_tibble(id) |>
    ggplot(aes(x = MetaMorpheus.score,
               colour = isDecoy)) +
    geom_density() +
    facet_wrap(~ spectrumFile)

table(id$isDecoy)

summary(id$PSM.level.q.value)

table(id$rank)
