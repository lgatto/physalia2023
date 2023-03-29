library(PSMatch)
library(limma)
library(msdata)
library(tidyverse)
library(QFeatures)

## quantitative assay: features x samples (5000 by 16)
## feature meta-data: rowData
## sample meta-data: colData

library(SummarizedExperiment)

data(se_na2)


assay(se_na2)

rowData(se_na2)

colData(se_na2)


data(feat1)
feat1

colData(feat1)

colData(feat1)$X <- c("X1", "X2")
colData(feat1)

feat1[[1]]

feat1[, , 1]

names(feat1)

feat1[["psms"]]


assay(feat1[["psms"]])

rowData(feat1[["psms"]])

feat1 <- aggregateFeatures(feat1,
                  i = "psms",
                  fcol = "Sequence",
                  name = "peptides",
                  fun = colMeans)

feat1

assay(feat1[["peptides"]])

rowData(feat1[["peptides"]])

## Exercise: aggregate peptides into proteins

feat1 <- aggregateFeatures(feat1,
                  i = "peptides",
                  fcol = "Protein",
                  name = "proteins",
                  fun = colMedians)

feat1

feat1["ProtA", , ]


filterFeatures(feat1, ~ pval < 0.05)

filterFeatures(feat1, ~ pval < 0.05, keep = TRUE)

filterFeatures(feat1, ~ location == "Mitochondrion")

## Ex: filter features dans do not localise to the MT

filterFeatures(feat1, ~ location != "Mitochondrion")


data(hlpsms)

hl <- readQFeatures(hlpsms,
                    ecol = 1:10,
                    name = "psms")

hl

assay(hl[[1]])
rowData(hl[[1]])
colData(hl)

se <- readSummarizedExperiment(hlpsms,
                               ecol = 1:10)

se

QFeatures(list(psms = se))

## Load data as a SummarizedExperiment

f <- msdata::quant(full.names = TRUE)

names(read.delim(f))

grep("Intensity\\.", names(read.delim(f)), value = TRUE)

i <- grep("Intensity\\.", names(read.delim(f)))


cptac_se <- readSummarizedExperiment(f, ecol = i,
                                     sep = "\t",
                                     fnames = "Sequence")
cptac_se

colnames(cptac_se)

colData(cptac_se)

## Ex: populate the colData

colnames(cptac_se) <- sub("^I.+\\.", "", colnames(cptac_se))

colData(cptac_se)$condition <- rep(c("A", "B"), each = 3)
colData(cptac_se)$id <- rep(7:9, 2)

names(rowData(cptac_se))

keep_var <-  c("Sequence", "Proteins", "Leading.razor.protein",
               "PEP", "Score", "Potential.contaminant", "Reverse")


rowData(cptac_se) <- rowData(cptac_se)[, keep_var]

## Missing data

anyNA(assay(cptac_se))

cptac_se <- zeroIsNA(cptac_se)

table(is.na(assay(cptac_se)))

nNA(cptac_se)

barplot(nNA(cptac_se)$nNAcols$nNA)

table(nNA(cptac_se)$nNArows$nNA)

cptac_se <- filterNA(cptac_se, pNA = 4/6)



data(se_na2)
se_na2

MsCoreUtils::imputeMethods()

assay(impute(se_na2, method = "knn"))

assay(impute(se_na2, method = "zero"))

table(rowData(se_na2)$randna)

impute(se_na2, method = "mixed",
       randna = rowData(se_na2)$randna,
       mar = "knn", mnar = "zero")

install.packages("imputeLCMD")

## Compare multiple imputation methods:
##  - knn, zero, MinDet, bpca
##  - plot density of values

plot(density(na.omit(assay(se_na2))))
lines(density(assay(impute(se_na2, method = "knn"))), col = "red")
lines(density(assay(impute(se_na2, method = "zero"))), col = "blue")
lines(density(assay(impute(se_na2, method = "MinDet"))), col = "green")
lines(density(assay(impute(se_na2, method = "bpca"))), col = "orange")

table(rowData(cptac_se)$Reverse)

table(rowData(cptac_se)$Potential.contaminant)

## Compare the score (and PEP) distributions for FWD and Reverse peptides

rowData(cptac_se) |>
    as_tibble() |>
    ggplot(aes(x = Score,
               colour = Reverse)) +
    geom_density()



rowData(cptac_se) |>
    as_tibble() |>
    ggplot(aes(x = PEP,
               colour = Reverse)) +
    geom_density()



summary(rowData(cptac_se)$PEP)

prot <- rowData(cptac_se)$Proteins
names(prot) <- rowData(cptac_se)$Sequence

am <- PSMatch::makeAdjacencyMatrix(prot, split = ";")

cc <- ConnectedComponents(am)

which(ncols(cc) > 2 & nrows(cc) > 2)

connectedComponents(cc, 24)

## QFeatures

cptac <- QFeatures(list(peptides = cptac_se))
cptac

colData(cptac) <- colData(cptac_se)

## Using the filterFeatures() function, filter out the reverse and
## contaminant hits, and also retain those that have a posterior error
## probability smaller than 0.05.

cptac <- cptac |>
    filterFeatures(~ Reverse != "+") |>
    filterFeatures(~ Potential.contaminant != "+") |>
    filterFeatures(~ PEP < 0.05)


plot(density(na.omit(assay(cptac[[1]]))))

cptac <- logTransform(cptac,
                      i = "peptides",
                      name = "log_peptides")

plot(density(na.omit(assay(cptac[[2]]))))

## normalize(cptac, i, name, method = "center.median")
cptac <- normalize(cptac,
                   i = "log_peptides",
                   name  = "lognorm_peptides",
                   method = "center.median")

cptac

par(mfrow = c(1, 3))

limma::plotDensities(assay(cptac[[1]]))
limma::plotDensities(assay(cptac[[2]]))
limma::plotDensities(assay(cptac[[3]]))

## Use median aggregation to aggregation peptides into protein
## values. Hint: use Leading.razor.protein variable.
##
## Careful: NA

cptac <- aggregateFeatures(cptac,
                           "lognorm_peptides",
                           name = "proteins_med",
                           fcol = "Leading.razor.protein",
                           fun = colMedians,
                           na.rm = TRUE)

## MsCoreUtils::robustSummary

cptac

library(factoextra)
library(patchwork)

colData(cptac)$condition <- rep(c("A", "B"), each = 3)

pca_pep <- cptac[["lognorm_peptides"]] |>
    filterNA() |>
    assay() |>
    t() |>
    prcomp(scale = TRUE, center = TRUE) |>
    fviz_pca_ind(habillage = cptac$condition)

pca_prot <- cptac[["proteins_med"]] |>
    filterNA() |>
    assay() |>
    t() |>
    prcomp(scale = TRUE, center = TRUE) |>
    fviz_pca_ind(habillage = cptac$condition)

pca_pep + pca_prot


plot(cptac)

normalize(cptac, "log_peptides",
          name = "logquantiles_peptides",
          method = "quantiles") |>
    aggregateFeatures(
        "logquantiles_peptides",
        name = "proteins_med2",
        fcol = "Leading.razor.protein",
        fun = colMedians,
        na.rm = TRUE) |>
    plot()

cptac["P02787ups|TRFE_HUMAN_UPS", ,
      c("lognorm_peptides", "proteins_med")] |>
    longFormat(colvars = "condition") |>
    as_tibble() |>
    ggplot(aes(
        x = colname,
        y = value,
        colour = rowname,
        shape = condition)) +
    geom_point(size = 3) +
    geom_line(aes(group = rowname)) +
    facet_grid(~ assay) +
    ggtitle("P02787ups|TRFE_HUMAN_UPS")

library(limma)

prots <- getWithColData(cptac, "proteins_med")

prots

colData(prots)

design <- model.matrix(~ prots$condition)
fit <- lmFit(assay(prots), design)
fit <- eBayes(fit)

res <- topTable(fit, coef = "prots$conditionB",
         number = Inf) |>
    rownames_to_column("protein") |>
    as_tibble() |>
    mutate(TP = grepl("UPS", protein))

k <- filter(res, is.na(t)) |>
    pull(1)

assay(prots[k, ])

res |>
    ggplot(aes(x = logFC,
               y = -log10(adj.P.Val),
               colour = TP)) +
    geom_point() +
    geom_vline(xintercept = c(-1, 1)) +
    geom_hline(yintercept = -log10(0.05)) +
    scale_color_manual(values = c("black","red"))
