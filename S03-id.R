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
