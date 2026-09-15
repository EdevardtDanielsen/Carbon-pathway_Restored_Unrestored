.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT  <- if (length(.file)) dirname(normalizePath(.file[1])) else getwd()
OUTDIR <- file.path(ROOT, "sem_output")
dir.create(OUTDIR, showWarnings = FALSE)

suppressPackageStartupMessages({library(nlme); library(MuMIn); library(piecewiseSEM)})
source(file.path(ROOT, "00_load_sem_data.R"))
OUT <- file.path(OUTDIR, "tableS4_dsep_claims.txt"); sink(OUT)
lmc <- lmeControl(opt="optim", msMaxIter=500, maxIter=500)

dat <- load_SEM_data(); dat$MOB_II <- dat$P18w7c; dat$Bac <- dat$Pa15
dat$MOB_I <- rowMeans(dat[, c("P16w7c","P16w5c")], na.rm=FALSE)
macro <- dat[dat$Dominance=="macrophyte", ]
PLFA_TISSUE <- c("P16w7c","P16w5c","P18w7c","P18_3n3","P18_2w6","Pi15","Pa15","MOB_I","MOB_II","Bac","d13C_Macro")
macro[macro$Mesocosm=="TH2", PLFA_TISSUE] <- NA; macro$Mesocosm <- droplevels(macro$Mesocosm)
phyto <- dat[dat$Dominance=="phytoplankton", ]; phyto$Mesocosm <- droplevels(phyto$Mesocosm)
z_col <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
phyto$Oxidized_N <- rowMeans(cbind(z_col(phyto$NO3_N), z_col(phyto$NO2_N)), na.rm=FALSE)

fit_lme <- function(formula, data) tryCatch(lme(formula, random=~1|Mesocosm, data=data,
             na.action=na.omit, method="REML", control=lmc), error=function(e) NULL)
dsep_p <- function(x, y, parents, data) {
  preds <- unique(c(x, parents, "Day"))
  sub <- na.omit(data[, unique(c(y, preds, "Mesocosm"))]); sub$Mesocosm <- droplevels(sub$Mesocosm)
  if (nrow(sub) < length(preds) + 3) return(NA_real_)
  f <- as.formula(paste(y, "~", paste(preds, collapse="+")))
  m <- tryCatch(lme(f, random=~1|Mesocosm, data=sub, na.action=na.omit, method="REML", control=lmc),
                error=function(e) tryCatch(lm(f, data=sub), error=function(e2) NULL))
  if (is.null(m)) return(NA_real_)
  tt <- if (inherits(m,"lme")) summary(m)$tTable else coef(summary(m))
  if (!(x %in% rownames(tt))) return(NA_real_)
  tt[x, if ("p-value" %in% colnames(tt)) "p-value" else 4]
}
m1_dsep <- list(
  c("Temp","d13C_DIC","TN"), c("TN","d13C_DOC","d13C_DIC,Temp"),
  c("Temp","d13C_CO2_Conc","d13C_DIC,TN"), c("d13C_DOC","d13C_CO2_Conc","d13C_DIC,TN,Temp"),
  c("Temp","d13C_CH4_Conc","d13C_DOC"), c("TN","d13C_CH4_Conc","d13C_DOC"),
  c("d13C_DIC","d13C_CH4_Conc","d13C_DOC"), c("d13C_CO2_Conc","d13C_CH4_Conc","d13C_DOC,d13C_DIC,TN"),
  c("Temp","d13C_Macro","d13C_DIC,d13C_CO2_Conc"), c("TN","d13C_Macro","d13C_DIC,d13C_CO2_Conc"),
  c("d13C_DOC","d13C_Macro","d13C_DIC,d13C_CO2_Conc,Temp"),
  c("d13C_CH4_Conc","d13C_Macro","d13C_DIC,d13C_CO2_Conc,d13C_DOC"),
  c("Temp","MOB_I","d13C_DOC"), c("TN","MOB_I","d13C_DOC"), c("d13C_DIC","MOB_I","d13C_DOC"),
  c("d13C_CH4_Conc","MOB_I","d13C_DOC"), c("d13C_CO2_Conc","MOB_I","d13C_DIC,d13C_DOC,TN"),
  c("d13C_Macro","MOB_I","d13C_DIC,d13C_DOC,d13C_CO2_Conc"),
  c("Temp","MOB_II","d13C_DIC"), c("TN","MOB_II","d13C_DIC"), c("d13C_DOC","MOB_II","d13C_DIC,Temp"),
  c("d13C_CH4_Conc","MOB_II","d13C_DIC,d13C_DOC"), c("d13C_CO2_Conc","MOB_II","d13C_DIC,TN"),
  c("d13C_Macro","MOB_II","d13C_DIC,d13C_CO2_Conc"), c("MOB_I","MOB_II","d13C_DIC,d13C_DOC"),
  c("Temp","Bac","MOB_I,MOB_II,d13C_DOC"), c("TN","Bac","MOB_I,MOB_II,d13C_DOC"),
  c("d13C_DIC","Bac","MOB_I,MOB_II,d13C_DOC"), c("d13C_CH4_Conc","Bac","MOB_I,MOB_II,d13C_DOC"),
  c("d13C_CO2_Conc","Bac","MOB_I,MOB_II,d13C_DOC,d13C_DIC,TN"),
  c("d13C_Macro","Bac","MOB_I,MOB_II,d13C_DOC,d13C_DIC,d13C_CO2_Conc"))
m2_dsep <- c(m1_dsep, list(
  c("Temp","d13C_Zoo","MOB_II,P18_2w6"), c("TN","d13C_Zoo","MOB_II,P18_2w6"),
  c("d13C_DIC","d13C_Zoo","MOB_II,P18_2w6"), c("d13C_DOC","d13C_Zoo","MOB_II,P18_2w6,Temp"),
  c("d13C_CH4_Conc","d13C_Zoo","MOB_II,P18_2w6"), c("d13C_CO2_Conc","d13C_Zoo","MOB_II,P18_2w6"),
  c("d13C_Macro","d13C_Zoo","MOB_II,P18_2w6,d13C_DIC,d13C_CO2_Conc"),
  c("MOB_I","d13C_Zoo","MOB_II,P18_2w6"), c("Bac","d13C_Zoo","MOB_II,P18_2w6")))

cat("=========== MACROPHYTE-DOMINATED (MACRO-zoo) d-sep basis (hand-rolled, m2_dsep) ===========\n")
macro_rows <- lapply(m2_dsep, function(tr){
  parents <- if (nchar(tr[3])) strsplit(tr[3],",")[[1]] else character(0)
  p <- dsep_p(tr[1], tr[2], parents, macro)
  data.frame(Claim=sprintf("%s _|_ %s | %s", tr[2], tr[1], tr[3]), P.Value=round(p,4), stringsAsFactors=FALSE)
})
macro_ds <- do.call(rbind, macro_rows)
print(macro_ds, row.names=FALSE)
pv <- macro_ds$P.Value[!is.na(macro_ds$P.Value) & macro_ds$P.Value>0]
cat(sprintf("\n  n_claims=%d  fails(<0.05)=%d  Fisher C=%.2f df=%d p=%.4f\n",
    length(pv), sum(pv<0.05), -2*sum(log(pv)), 2*length(pv), 1-pchisq(-2*sum(log(pv)),2*length(pv))))

cat("\n=========== PHYTOPLANKTON-DOMINATED (PHYTO-zoo) d-sep basis (piecewiseSEM::dSep) ===========\n")
PHYTO_zoo <- list(
  P1=d13C_DIC~PO4_P+Day, PDOC=d13C_DOC~d13C_DIC+Day, P2=d13C_CO2_Conc~d13C_DIC+MOB_I+Day,
  P3=d13C_CH4_Conc~Chl_a+Oxidized_N+d13C_Zoo+Day, P4=Pa15~MOB_I+Chl_a+Day,
  P5=MOB_II~d13C_DIC+Day, P6=P18_2w6~d13C_DIC+Chl_a+Day, P7=d13C_Zoo~P18_2w6+Pa15+Day)
pm <- Filter(function(m) inherits(m,"lme"), lapply(PHYTO_zoo, fit_lme, data=phyto))
mdl <- do.call(psem, unname(pm))
ds <- dSep(mdl, .progressBar=FALSE)
print(ds[, c("Independ.Claim","P.Value")], row.names=FALSE)
f0 <- fisherC(mdl)
cat(sprintf("\n  n_claims=%d  fails(<0.05)=%d  Fisher C=%.2f df=%d p=%.4f\n",
    nrow(ds), sum(ds$P.Value<0.05,na.rm=TRUE), f0$Fisher.C, f0$df, f0$P.Value))
cat("\nDONE\n"); sink()
