rm(list = ls())

####################################################################################################
# Install required packages
####################################################################################################

source("src/setup.R")
.setup("pwr")

####################################################################################################
# Load data
####################################################################################################

load("../data/raw_data.RData")

td  <- concentrazioni[as.vector(clinici$sample.id[clinici$tipo == "controllo sano"]),]
asd <- concentrazioni[as.vector(clinici$sample.id[clinici$tipo == "patologico"]),]

rm(clinici, concentrazioni)

load("../data/classifier.data.RData")

most.significant.metabolites <- colnames(data[,-c(1,2)])

rm(data)

td  <- td[,most.significant.metabolites]
asd <- asd[,most.significant.metabolites]

####################################################################################################
# Run power analysis
####################################################################################################
n1      <- nrow(td)
n2      <- nrow(asd)
sig.lvl <- 0.05

r <- list(
    pwr.t2n.test(n1=n1,n2=n2,sig.level=sig.lvl,
                 d=mean(abs(colMeans(td)-colMeans(asd))/unlist(Map(max, apply(td, 2, sd), apply(asd, 2, sd))))),
    pwr.t2n.test(n1=n1,n2=n2,sig.level=sig.lvl,
                 d=mean(abs(colMeans(td)-colMeans(asd))/unlist(Map(mean, apply(td, 2, sd), apply(asd, 2, sd))))),
    pwr.t2n.test(n1=n1,n2=n2,sig.level=sig.lvl,
                 d=mean(abs(colMeans(td)-colMeans(asd))/unlist(Map(min, apply(td, 2, sd), apply(asd, 2, sd)))))
)

print(r[[which.min(sapply(r, function (x) (x$power)))]])


