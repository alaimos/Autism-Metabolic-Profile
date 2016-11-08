####################################################################################################
# Load required files and data
####################################################################################################
rm(list = ls())
source("src/setup.R")
source("src/analyze.R")
load("../data/raw_data.RData")
####################################################################################################
# Some constants
####################################################################################################
p.values.thresh <- 0.05 # p-value threshold
####################################################################################################
# Filters and prepare raw data
####################################################################################################
common.ids     <- intersect(rownames(clinici), rownames(concentrazioni))
clinici        <- clinici[common.ids,]
concentrazioni <- concentrazioni[common.ids, ]
concentrazioni <- t(data.matrix(concentrazioni))
concentrazioni[concentrazioni == 0] <- NA
log.correction <- function (c) {
    c[is.na(c)] <- 0
    return (log2(c+1))
}
concentrazioni <- log.correction(concentrazioni)
clinici$sample.id <- as.vector(clinici$sample.id)
tmp <- character(length(clinici$eta))
tmp[clinici$eta <= 5] <- "M5"
tmp[clinici$eta > 5]  <- "P5"
tmp[tmp == ""] <- NA
clinici$strat.eta <- factor(tmp)
tmp <- paste0(clinici$tipo, clinici$strat.eta)
tmp[is.na(clinici$strat.eta)] <- NA
tmp <- gsub("patologico", "P", tmp, fixed = TRUE)
tmp <- gsub("controllo sano", "S", tmp, fixed = TRUE)
clinici$strat.eta.tipo <- factor(tmp)
clean.factor <- function (f, na.value = NA) {
    tmp <- gsub(" ", "", tolower(as.vector(f)), fixed = TRUE)
    tmp[tmp == "si"] <- "S"
    tmp[tmp == "no"] <- "N"
    tmp[is.na(tmp)] <- na.value
    return (factor(tmp))
}
clinici$presenza.disturbi.gi <- clean.factor(clinici$presenza.disturbi.gi, "C")
clinici$regressione <- clean.factor(clinici$regressione, "C")
clinici$rm          <- clean.factor(clinici$rm, "C")
tmp <- paste0(clinici$regressione, clinici$strat.eta)
tmp[is.na(clinici$strat.eta)] <- NA
clinici$strat.eta.regr <- factor(tmp)
tmp <- as.vector(clinici$CSS)
tmp1 <- character(length(tmp))
tmp1[tmp <= 3] <- "NS"
tmp1[tmp > 3 & tmp <= 5] <- "AS"
tmp1[tmp > 5]  <- "A"
clinici$severity <- factor(tmp1)
tmp <- as.vector(clinici$qit)
tmp1 <- character(length(tmp))
tmp1[tmp > 70] <- "NR"
tmp1[tmp >= 55 & tmp <= 70] <- "L"
tmp1[tmp >= 40 & tmp < 55] <- "M"
tmp1[tmp >= 25 & tmp < 40] <- "G"
tmp1[tmp <  25]  <- "GG"
tmp1[tmp1 == ""] <- "NR"
clinici$rm.cat <- factor(tmp1)
rm(common.ids, tmp,tmp1)
####################################################################################################
# Find altered metabolites in ASD vs TD
####################################################################################################
analyze(
    clinici = clinici, 
    concentrazioni = concentrazioni, 
    model.formula = ~0+tipo, 
    levels.vector = c("C", "P"),
    levels.variable = NULL, 
    contrasts = "P-C", 
    out.de.file = "../results/altered-metabolites/ASD_vs_TD.txt"
)

####################################################################################################
# Find altered metabolites in ASD vs TD
####################################################################################################
analyze(
    clinici = clinici, 
    concentrazioni = concentrazioni, 
    model.formula = ~0+tipo, 
    levels.vector = c("C", "P"),
    levels.variable = NULL, 
    contrasts = "P-C", 
    out.de.file = "../results/altered-metabolites/ASD_vs_TD.txt"
)
####################################################################################################
# Regressive Autism and Age
####################################################################################################
tmp <- paste0(clinici$strat.eta,clinici$regressione)
tmp[tmp == "NAC"] <- NA
clinici$regr.eta <- factor(tmp)
analyze(
    clinici = clinici, 
    concentrazioni = concentrazioni, 
    model.formula = ~0+regr.eta, 
    levels.variable = "regr.eta",
    contrasts = c("M5S-M5C","M5N-M5C","P5S-P5C","P5N-P5C","M5S-M5N","P5S-P5N"), 
    out.de.file = "../results/altered-metabolites/RegrAutism_and_Age.txt"
)
####################################################################################################
# Intellectual disability and Age
####################################################################################################
tmp <- paste0(clinici$strat.eta,clinici$rm)
tmp[tmp == "NAC"] <- NA
clinici$rm.eta <- factor(tmp)
analyze(
    clinici = clinici, 
    concentrazioni = concentrazioni, 
    model.formula = ~0+rm.eta, 
    levels.variable = "rm.eta",
    contrasts = c("M5S-M5C","M5N-M5C","P5S-P5C","P5N-P5C","M5S-M5N","P5S-P5N"), 
    out.de.file = "../results/altered-metabolites/id_and_age.txt"
)

####################################################################################################
# Severity Score and Age
####################################################################################################
tmp <- paste0(clinici$strat.eta,clinici$severity)
clinici$sev.eta <- factor(tmp)
analyze(
    clinici = clinici, 
    concentrazioni = concentrazioni, 
    model.formula = ~0+sev.eta, 
    levels.variable = "sev.eta",
    contrasts = c("M5A-M5NS","M5AS-M5NS","P5A-P5NS","M5AS-M5NS"), 
    out.de.file = "../results/altered-metabolites/severity_and_age.txt"
)
