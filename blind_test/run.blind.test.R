rm(list = ls())

source("src/setup.R")
source("src/algorithms.R")


####################################################################################################
# Runs a blind test on a single classification algorithm
####################################################################################################

run.single.test <- function (classifier, dataset, n.repeat = 1000, 
                             test.set.fraction = 1/3, cases.fraction = 1/2, ...) {
    set.seed(1234) # Set random seed for reproducibility
    cn    <- colnames(data)
    dataset$classi <- factor(as.vector(dataset$classi))
    ## Change as needed
    control.name <- levels(dataset$classi)[1]
    case.name    <- levels(dataset$classi)[2]
    ## END
    cat("Case Name:",case.name,"Control Name:",control.name,"\n")
    
    dataset.controls <- dataset[dataset$classi == control.name,]
    dataset.cases    <- dataset[dataset$classi == case.name,]
    test.set.size        <- round(nrow(dataset) * test.set.fraction)
    train.set.size       <- nrow(dataset) - test.set.size
    number.test.cases    <- round(test.set.size * cases.fraction)
    number.test.controls <- test.set.size - number.test.cases
    train.sets <- vector("list", n.repeat)
    test.sets  <- vector("list", n.repeat)
    tables     <- vector("list", n.repeat)
    accuracy   <- numeric(n.repeat)
    TPR        <- numeric(n.repeat)
    TNR        <- numeric(n.repeat)
    FPR        <- numeric(n.repeat)
    FNR        <- numeric(n.repeat)
    PPV        <- numeric(n.repeat)
    NPV        <- numeric(n.repeat)
    LRP        <- numeric(n.repeat)
    LRM        <- numeric(n.repeat)
    DOR        <- numeric(n.repeat)
    
    for (i in 1:n.repeat) {
        repeat {
            test.cases    <- sample(1:nrow(dataset.cases), number.test.cases)
            test.controls <- sample(1:nrow(dataset.controls), number.test.controls)
            test.set  <- rbind(dataset.cases[test.cases,], dataset.controls[test.controls,])
            test.set$classi  <- factor(as.vector(test.set$classi))
            train.set <- rbind(dataset.cases[-test.cases,], dataset.controls[-test.controls,])
            train.set$classi <- factor(as.vector(train.set$classi))
            prediction <- classifier(train.set, test.set, ...)
            if (!all(is.null(prediction))) {
                cat(".")
                break
            }
        }
        if (!(control.name %in% levels(prediction))) {
            levels(prediction) <- c(levels(prediction), control.name)
        }
        if (!(case.name %in% levels(prediction))) {
            levels(prediction) <- c(levels(prediction), case.name)
        }
        test.set$prediction <- prediction
        train.sets[[i]] <- train.set
        test.sets[[i]]  <- test.set
        tables[[i]]     <- table(test.set$classi, test.set$prediction, dnn = c("real", "prediction"))
        tables[[i]]     <- tables[[i]][c(control.name, case.name),c(control.name, case.name)]
        
        accuracy[i] <- sum(diag(tables[[i]])) / sum(tables[[i]])
        TPR[i]      <- tables[[i]][2,2] / (tables[[i]][2,2]+tables[[i]][2,1])
        TNR[i]      <- tables[[i]][1,1] / (tables[[i]][1,1]+tables[[i]][1,2])
        FPR[i]      <- tables[[i]][1,2] / (tables[[i]][1,1]+tables[[i]][1,2])
        FNR[i]      <- tables[[i]][2,1] / (tables[[i]][2,2]+tables[[i]][2,1])
        PPV[i]      <- tables[[i]][2,2] / (tables[[i]][2,2]+tables[[i]][1,2])
        NPV[i]      <- tables[[i]][1,1] / (tables[[i]][1,1]+tables[[i]][2,1])
        LRP[i]      <- TPR[i]/FPR[i]
        LRM[i]      <- FNR[i]/TNR[i]
        DOR[i]      <- LRP[i]/LRM[i]
    }
    cat("\n")
    
    print.cint <- function (vals) {
        vals <- vals[!is.na(vals) & !is.infinite(vals)]
        m    <- format(round(mean(vals), digits = 4), scientific = FALSE)
        tmp  <- tryCatch(t.test(vals), error=function(e) {return (list(conf.int=c(NaN, NaN))); })
        cint <- format(round(tmp$conf.int, digits = 4), scientific = FALSE)
        cint <- paste0("[", cint[1],"; ", cint[2], "]")
        return (paste0(m,"  ",cint))
    }
    
    mean.table <- Reduce("+", tables)/n.repeat
    mean.table.fraction <- mean.table
    for (j in 1:nrow(mean.table.fraction)) {
        mean.table.fraction[j,] <- mean.table.fraction[j,] / sum(mean.table.fraction[j,]) 
    }
    return (list(
        summary=list(accuracy=print.cint(accuracy),
                     TPR=print.cint(TPR),
                     TNR=print.cint(TNR),
                     FPR=print.cint(FPR),
                     FNR=print.cint(FNR),
                     PPV=print.cint(PPV),
                     NPV=print.cint(NPV),
                     DOR=print.cint(DOR)
        ),
        table=mean.table,
        table.fraction=mean.table.fraction,
        single.run=list(train.sets=train.sets,
                        test.sets=test.sets,
                        tables=tables,
                        accuracy=accuracy,
                        TPR=TPR,
                        TNR=TNR,
                        FPR=FPR,
                        FNR=FNR,
                        LRP=LRP,
                        LRM=LRM,
                        DOR=DOR
        )
    ))
}


####################################################################################################
# Some definitions
####################################################################################################

nr                <- 1000      # number of repeats
test.set.fraction <- 0.2407407 # Percentage computed from the original training/test set
cases.fraction    <- 0.4358974 # Percentage computed from the original test set

####################################################################################################
# Classifies patients youger than 5 y.o
####################################################################################################

load("../data/classifier.data.RData")
data <- data[data$eta < 5 & !is.na(data$eta), -c(2)]

r  <- run.single.test(std.ctree, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r1 <- run.single.test(std.randfor, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r2 <- run.single.test(std.svm, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r3 <- run.single.test(std.lm, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r4 <- run.single.test(std.rpart, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r5 <- run.single.test(std.svm.radial, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r6 <- run.single.test(std.lrm, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r7 <- run.single.test(std.naiveBayes, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)

finale <- rbind(r$summary, r1$summary, r2$summary, 
                r5$summary, r7$summary, r3$summary, 
                r6$summary, r4$summary)
rownames(finale) <- c("C-Tee", "RF", "SVM Linear", "SVM Radial Basis", "Naive Bayes", "LM", "LRM", "PART")

# Save results
write.table(finale, file = "../results/blind-test/results-lt-5yo.txt", sep = "\t", append = FALSE, quote = FALSE, dec = ".", row.names = TRUE, col.names = TRUE)
save(r,r1,r2,r3,r4,r5,r6,r7, file = "../results/blind-test/results-lt-5yo.RData", ascii = FALSE, compress = TRUE)

####################################################################################################
# Classifies patients older than 5 y.o
####################################################################################################

load("../data/classifier.data.RData")
data <- data[data$eta > 5 & !is.na(data$eta), -c(2)]

r  <- run.single.test(std.ctree, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r1 <- run.single.test(std.randfor, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r2 <- run.single.test(std.svm, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r3 <- run.single.test(std.lm, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r4 <- run.single.test(std.rpart, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r5 <- run.single.test(std.svm.radial, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r6 <- run.single.test(std.lrm, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r7 <- run.single.test(std.naiveBayes, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)

finale <- rbind(r$summary, r1$summary, r2$summary, 
                r5$summary, r7$summary, r3$summary, 
                r6$summary, r4$summary)
rownames(finale) <- c("C-Tee", "RF", "SVM Linear", "SVM Radial Basis", "Naive Bayes", "LM", "LRM", "PART")

# Save results
write.table(finale, file = "../results/blind-test/results-gt-5yo.txt", sep = "\t", append = FALSE, quote = FALSE, dec = ".", row.names = TRUE, col.names = TRUE)
save(r,r1,r2,r3,r4,r5,r6,r7, file = "../results/blind-test/results-gt-5yo.RData", ascii = FALSE, compress = TRUE)

####################################################################################################
# Classifies all patients
####################################################################################################

load("../data/classifier.data.RData")
data <- data[, -c(2)]

r  <- run.single.test(std.ctree, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r1 <- run.single.test(std.randfor, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r2 <- run.single.test(std.svm, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r3 <- run.single.test(std.lm, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r4 <- run.single.test(std.rpart, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r5 <- run.single.test(std.svm.radial, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r6 <- run.single.test(std.lrm, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)
r7 <- run.single.test(std.naiveBayes, data, n.repeat = nr, test.set.fraction=test.set.fraction, cases.fraction = cases.fraction)

finale <- rbind(r$summary, r1$summary, r2$summary, 
                r5$summary, r7$summary, r3$summary, 
                r6$summary, r4$summary)
rownames(finale) <- c("C-Tee", "RF", "SVM Linear", "SVM Radial Basis", "Naive Bayes", "LM", "LRM", "PART")

# Save results
write.table(finale, file = "../results/blind-test/results-all-patients.txt", sep = "\t", append = FALSE, quote = FALSE, dec = ".", row.names = TRUE, col.names = TRUE)
save(r,r1,r2,r3,r4,r5,r6,r7, file = "../results/blind-test/results-all-patients.RData", ascii = FALSE, compress = TRUE)
