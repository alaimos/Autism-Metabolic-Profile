####################################################################################################
# Required libraries
####################################################################################################

.setup(c(
    "e1071",
    "party",
    "rpart",
    "randomForest",
    "bnlearn",
    "ROCR"
))

####################################################################################################
# Utilities
####################################################################################################

#Change class names into numeric values
change.classes <- function (train) {
    classes <- levels(train$classi)
    map <- 0:(length(classes) - 1)
    names(map) <- classes
    train$classi <- unname(map[train$classi])
    return (list(new.train=train, map=map))
}

#Revert class numeric values into class names
revert.classes <- function (cv, map) {
    nc <- character(length(cv))
    for (i in 1:length(map)) {
        nc[cv == unname(map[i])] <- names(map)[i]
    }
    return (factor(nc))
}

#Estimate best cutoff on the basis of a sensitivity/specificity curve
estimate.best.cutoff <- function (model, train) {
    pred <- prediction(unname(model$fitted.values), 
                       train[names(model$fitted.values),"classi"])
    perf1 <- performance(pred, "sens")
    perf2 <- performance(pred, "spec")
    intersection <- which.min(abs(unlist(perf1@y.values) - unlist(perf2@y.values)))
    return (unlist(perf1@x.values)[intersection])
}

#Predict classes by employing a cutoff value
cutoff.predict <- function (model, test, cutoff, ...) {
    preds <- predict(model, test, ...)
    preds[preds <= cutoff] <- 0
    preds[preds >  cutoff] <- 1
    return (preds)
}

####################################################################################################
# Implementation of standard classifiers functions. These functions class a classifier algorithm on
# a training set, and return predictions made on a test set.
####################################################################################################

## C-tree
std.ctree <- function (train, test) {
    ct <- ctree(classi ~ ., data = train)
    return (predict(ct, test))
}

## PART
std.rpart <- function (train, test) {
    ct <- rpart(classi ~ ., data = train, method = "class")
    return (predict(ct, test, type = "class"))
}

## Random Forest
std.randfor <- function (train, test) {
    ct <- randomForest(classi ~ ., data = train)
    return (predict(ct, test))
}

## SVM Linear Kernel
std.svm <- function (train, test, return.model = FALSE) {
    svm.model.l <<- svm(classi ~ ., data = train, kernel="linear")
    if (return.model) return (svm.model.l)
    return (predict(svm.model.l, test))
}

## SVM Radial-basis Kernel
std.svm.radial <- function (train, test) {
    svm.model.l <- svm(classi ~ ., data = train, kernel="radial")
    return (predict(svm.model.l, test))
}

## Linear Regression Model
std.lm <- function (train, test) {
    tmp <- change.classes(train)
    model <- lm(classi~., data = tmp$new.train)
    preds <- cutoff.predict(model, test, estimate.best.cutoff(model, tmp$new.train))
    return (revert.classes(preds, tmp$map))
}

## Logistic Regression Model
std.lrm <- function (train, test) {
    tmp <- change.classes(train)
    model <- glm(classi~., data = tmp$new.train, family = binomial(link = "logit"))
    preds <- cutoff.predict(model, test, estimate.best.cutoff(model, tmp$new.train), type="response")
    return (revert.classes(preds, tmp$map))
}

## Naive-Bayes Model
std.naiveBayes <- function (train, test) {
    model <- naiveBayes(classi ~ ., train)
    preds <- predict(model, test)
    return (preds)
}