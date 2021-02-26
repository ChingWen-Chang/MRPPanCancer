mafSurvival1<-function (maf, genes = NULL, samples = NULL, clinicalData = NULL, 
          time = "Time", Status = "Status", groupNames = c("Mutant", 
                                                           "WT"), showConfInt = TRUE, addInfo = TRUE, col = c("maroon", 
                                                                                                              "royalblue"), isTCGA = FALSE, textSize = 12, fn = NULL, 
          width = 6, height = 6) 
{
  if (is.null(genes) & is.null(samples)) {
    stop("Either provide Gene names or Sample names to group by.")
  }
  if (!is.null(genes) & !is.null(samples)) {
    stop("Either provide Gene names or Sample names to group by. Not both!")
  }
  if (is.null(clinicalData)) {
    message("Looking for clinical data in annoatation slot of MAF..")
    clinicalData = getClinicalData(x = maf)
    clinicalData = data.table::setDT(clinicalData)
  }
  else {
    clinicalData = data.table::setDT(clinicalData)
  }
  if (!"Tumor_Sample_Barcode" %in% colnames(clinicalData)) {
    print(colnames(clinicalData))
    stop("Column Tumo_Sample_Barcode not found in clinical data. Check column names and rename it to Tumo_Sample_Barcode if necessary.")
  }
  if (isTCGA) {
    clinicalData$Tumor_Sample_Barcode = substr(x = clinicalData$Tumor_Sample_Barcode, 
                                               start = 1, stop = 12)
  }
  if (length(colnames(clinicalData)[colnames(clinicalData) %in% 
                                    time]) == 0) {
    print(colnames(clinicalData))
    stop(paste0(time, " not found in clinicalData. Use argument time to povide column name containing time to event."))
  }
  else {
    colnames(clinicalData)[colnames(clinicalData) %in% time] = "Time"
  }
  if (length(colnames(clinicalData)[colnames(clinicalData) %in% 
                                    Status]) == 0) {
    print(colnames(clinicalData))
    stop(paste0(Status, " not found in clinicalData. Use argument Status to povide column name containing events (Dead or Alive)."))
  }
  else {
    colnames(clinicalData)[colnames(clinicalData) %in% Status] = "Status"
  }
  if (!is.null(genes)) {
    genesTSB = genesToBarcodes(maf = maf, genes = genes, 
                               justNames = TRUE)
    genesTSB = genesTSB[sapply(genesTSB, FUN = function(x) length(x) != 
                                 0)]
    message("Number of mutated samples for given genes: ")
    print(sapply(genesTSB, FUN = length))
    genesMissing = genes[!genes %in% names(genesTSB)]
    if (length(genesMissing) > 0) {
      genes = genes[!genes %in% genesMissing]
      genesMissing = paste(genesMissing, collapse = ", ")
      message(paste0("genes ", genesMissing, " does not seeem to be mutated. Removing them."))
    }
    if (length(genes) == 0) {
      stop("None of the given genes are mutated!")
    }
    else {
      genes = paste(genes, collapse = ", ")
    }
    genesTSB = unique(as.character(unlist(genesTSB)))
  }
  else {
    genesTSB = samples
  }
  data.table::setDT(clinicalData)
  clinicalData$Time = suppressWarnings(as.numeric(as.character(clinicalData$Time)))
  clinicalData$Status = suppressWarnings(as.integer(clinicalData$Status))
  clinicalData$Group = ifelse(test = clinicalData$Tumor_Sample_Barcode %in% 
                                genesTSB, yes = groupNames[1], no = groupNames[2])
  clin.mut.dat = clinicalData[, .(medianTime = median(Time, 
                                                      na.rm = TRUE), N = .N), Group][order(Group)]
  message("Median survival..")
  print(clin.mut.dat)
  clinicalData$Time = ifelse(test = is.infinite(clinicalData$Time), 
                             yes = 0, no = clinicalData$Time)
  if (nrow(clinicalData[!is.na(Time)][!is.na(Status)]) < nrow(clinicalData)) {
    message(paste0("Removed ", nrow(clinicalData) - nrow(clinicalData[!is.na(Time)][!is.na(Status)]), 
                   " samples with NA's"))
    clinicalData = clinicalData[!is.na(Time)][!is.na(Status)]
  }
  surv.km = survival::survfit(formula = survival::Surv(time = Time, 
                                                       event = Status) ~ Group, data = clinicalData, conf.type = "log-log")
  res = summary(surv.km)
  surv.diff = survival::survdiff(formula = survival::Surv(time = Time, 
                                                          event = Status) ~ Group, data = clinicalData)
  surv.diff.pval = signif(1 - pchisq(surv.diff$chisq, length(surv.diff$n) - 
                                       1), digits = 3)
  surv.cox = survival::coxph(formula = survival::Surv(time = Time, 
                                                      event = Status) ~ Group, data = clinicalData)
  hr = signif(1/exp(stats::coef(surv.cox)), digits = 3)
  surv.dat = data.table::data.table(Group = res$strata, Time = res$time, 
                                    survProb = res$surv, survUp = res$upper, survLower = res$lower)
  surv.dat$Group = gsub(pattern = "Group=", replacement = "", 
                        x = surv.dat$Group)
  par(mar = c(4, 4, 2, 2))
  x_lims = pretty(surv.km$time)
  y_lims = seq(0, 1, 0.2)
  plot(NA, xlim = c(min(x_lims), max(x_lims)), ylim = c(0, 
                                                        1), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
  abline(h = y_lims, v = x_lims, lty = 2, col = grDevices::adjustcolor(col = "gray70", 
                                                                       alpha.f = 0.5), lwd = 0.75)
  points(surv.dat[Group %in% "Mutant", Time], y = surv.dat[Group %in% 
                                                             "Mutant", survProb], pch = 8, col = col[1])
  points(surv.dat[Group %in% "WT", Time], y = surv.dat[Group %in% 
                                                         "WT", survProb], pch = 8, col = col[2])
  lines(surv.km[1], col = col[1], lty = 1, lwd = 2, conf.int = FALSE)
  lines(surv.km[2], col = col[2], lty = 1, lwd = 2, conf.int = FALSE)
  axis(side = 1, at = x_lims)
  axis(side = 2, at = y_lims, las = 2)
  mtext(text = "Survival probability", side = 2, line = 2.5)
  mtext(text = "Time", side = 1, line = 2)
  legend(x = "topright", legend = c("Mutant", "WT"), col = col, 
         bty = "n", lwd = 2, pch = 8)
  title(main = NA, sub = paste0("P-value: ", surv.diff.pval, 
                                "; ", "HR: ", hr), cex.main = 1, font.main = 4, col.main = "black", 
        cex.sub = 1, font.sub = 3, col.sub = ifelse(test = surv.diff.pval < 
                                                      0.05, yes = "red", no = "black"), line = 2.5, adj = 0)
  title(main = paste0(groupNames[1], " v/s ", groupNames[2]), 
        adj = 0, font.main = 4)
  return(clinicalData)
}
<bytecode: 0x7ff0f6eda800>
  <environment: namespace:maftools>