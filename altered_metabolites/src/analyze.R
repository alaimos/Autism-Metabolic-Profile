####################################################################################################
# Required libraries
####################################################################################################

.setup(c(
    "limma"
))

####################################################################################################
# Helper function used to get altered metabolites
####################################################################################################

analyze <- function (clinici, concentrazioni, model.formula, levels.variable, 
                     contrasts, out.de.file, levels.vector = NULL, p.values.thresh = 0.05, 
                     adj="BH") {
    
    design <- model.matrix(model.formula, data = clinici)
    if (all(!is.null(levels.vector)) && length(levels.vector) > 0) {
        colnames(design) <- levels.vector
    } else {
        colnames(design) <- levels(clinici[, levels.variable])
    }
    if (nrow(design) != ncol(concentrazioni)) {
        concentrazioni <- concentrazioni[,rownames(design)]
    }
    first.fit <- eBayes(lmFit(concentrazioni, design))
    cntrs <- makeContrasts(contrasts = contrasts, levels=design)
    secon.fit <- eBayes(contrasts.fit(first.fit, cntrs))
    
    if (length(contrasts) == 1) {
        all.de <- topTable(secon.fit, number = nrow(concentrazioni), adjust.method = adj,
                           p.value = p.values.thresh)
        write.table(all.de, file=out.de.file, row.names=TRUE, col.names=TRUE, quote = FALSE, sep="\t", dec = ",")
    } else {
        all.de <- vector("list", length(contrasts))
        names(all.de) <- contrasts
        fz <- file(out.de.file, "w")
        for (i in 1:length(contrasts)) {
            cat(contrasts[i], "\n", file = fz, append = TRUE, sep = "")
            all.de[[i]] <- topTable(secon.fit, number = nrow(concentrazioni), adjust.method = adj,
                                    p.value = p.values.thresh, coef = i)
            suppressWarnings(write.table(all.de[[i]], file=fz, append = TRUE, row.names=TRUE, col.names=TRUE, 
                                         quote = FALSE, sep="\t", dec = ","))
        }
        close(fz)
    }
    return (all.de)
}

get.all.de.metab <- function (all.de) {
    if (is.data.frame(all.de)) {
        return (rownames(all.de))
    } else {
        return (unique(Reduce(union, Map(rownames, tmp))))
    }
}