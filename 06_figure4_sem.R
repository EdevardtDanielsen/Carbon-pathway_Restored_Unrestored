.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT  <- if (length(.file)) dirname(normalizePath(.file[1])) else getwd()
OUTDIR <- file.path(ROOT, "sem_output")
dir.create(OUTDIR, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(nlme); library(MuMIn); library(MCMCglmm); library(coda); library(piecewiseSEM)
})

source(file.path(ROOT, "00_load_sem_data.R"))

set.seed(20260707L)
B     <- 1000L
NITT  <- 33000L; BURNIN <- 3000L; THIN <- 10L
lmc   <- lmeControl(opt="optim", msMaxIter=500, maxIter=500)

dat        <- load_SEM_data()
dat$MOB_II <- dat$P18w7c
dat$Bac    <- dat$Pa15
dat$MOB_I  <- rowMeans(dat[, c("P16w7c","P16w5c")], na.rm=FALSE)

macro <- dat[dat$Dominance=="macrophyte", ]
PLFA_TISSUE <- c("P16w7c","P16w5c","P18w7c","P18_3n3","P18_2w6","Pi15","Pa15",
                 "MOB_I","MOB_II","Bac","d13C_Macro")
macro[macro$Mesocosm=="TH2", PLFA_TISSUE] <- NA
macro$Mesocosm <- droplevels(macro$Mesocosm)

phyto <- dat[dat$Dominance=="phytoplankton", ]
phyto$Mesocosm <- droplevels(phyto$Mesocosm)
z_col <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
phyto$Oxidized_N <- rowMeans(cbind(z_col(phyto$NO3_N), z_col(phyto$NO2_N)), na.rm=FALSE)

TS   <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
LOG  <- file.path(OUTDIR, paste0("figure4_sem_", TS, ".txt"))
CSV  <- file.path(OUTDIR, paste0("figure4_sem_coefs_", TS, ".csv"))
FCSV <- file.path(OUTDIR, paste0("figure4_sem_fisher_", TS, ".csv"))
SCSV <- file.path(OUTDIR, paste0("figure4_sem_stdcoef_", TS, ".csv"))
sink(LOG, split=TRUE)
cat("Figure 4 SEM   (freq+bootstrap+Bayes+Fisher; both panels)\n TS:", TS, "\n")
cat("MACRO rows=", nrow(macro), " mesocosms=", paste(levels(macro$Mesocosm),collapse=","),
    "   PHYTO rows=", nrow(phyto), "\n\n")

fit_lme <- function(formula, data) {
  tryCatch(lme(formula, random=~1|Mesocosm, data=data, na.action=na.omit,
               method="REML", control=lmc),
           error=function(e) NULL)
}

resample_within_mesocosm <- function(dt) {
  bs  <- split(seq_len(nrow(dt)), dt$Mesocosm)
  idx <- unlist(lapply(bs, function(i) if (length(i)) sample(i, length(i), replace=TRUE) else integer(0)))
  out <- dt[idx,,drop=FALSE]; out$Mesocosm <- droplevels(factor(out$Mesocosm)); out
}

boot_ci <- function(formula, data, B) {
  cc <- na.omit(data[, c(all.vars(formula), "Mesocosm")]); cc$Mesocosm <- droplevels(cc$Mesocosm)
  fit0 <- tryCatch(lme(formula, random=~1|Mesocosm, data=cc, method="REML", control=lmc), error=function(e) NULL)
  if (is.null(fit0)) return(NULL)
  trm <- rownames(summary(fit0)$tTable)
  M <- matrix(NA_real_, B, length(trm), dimnames=list(NULL, trm))
  for (b in seq_len(B)) {
    rs <- resample_within_mesocosm(cc)
    fb <- tryCatch(suppressWarnings(lme(formula, random=~1|Mesocosm, data=rs, method="REML", control=lmc)),
                   error=function(e) NULL)
    if (!is.null(fb)) { tt <- summary(fb)$tTable; cm <- intersect(trm, rownames(tt)); M[b, cm] <- tt[cm, "Value"] }
  }
  data.frame(term=trm,
             q025=apply(M,2,quantile,0.025,na.rm=TRUE),
             q975=apply(M,2,quantile,0.975,na.rm=TRUE),
             boot_n=apply(M,2,function(x) sum(!is.na(x))), row.names=NULL)
}

prior_PE <- list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000)))
bayes_fit <- function(formula, data) {
  cc <- na.omit(data[, c(all.vars(formula), "Mesocosm")]); cc$Mesocosm <- droplevels(cc$Mesocosm)
  if (nrow(cc) < length(all.vars(formula))+2) return(NULL)
  tryCatch(suppressWarnings(MCMCglmm(fixed=formula, random=~Mesocosm, data=cc, prior=prior_PE,
              nitt=NITT, burnin=BURNIN, thin=THIN, verbose=FALSE, pr=FALSE)),
           error=function(e) NULL)
}

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
fisher_C_hand <- function(ps) { v <- ps[!is.na(ps) & ps > 0]
  if (!length(v)) return(NULL)
  list(Fisher.C=-2*sum(log(v)), df=2*length(v), P.Value=1-pchisq(-2*sum(log(v)),2*length(v)),
       k=length(v), n_fails=sum(v<0.05)) }

run_variant <- function(vlabel, forms, data, dsep=NULL) {
  cat("\n", paste(rep("=",74),collapse=""), "\n ", vlabel, "\n",
      paste(rep("=",74),collapse=""), "\n", sep="")
  rows <- list(); models <- list()
  for (k in names(forms)) {
    f <- forms[[k]]
    m <- fit_lme(f, data); models[[k]] <- m
    if (is.null(m) || !inherits(m,"lme")) { cat(sprintf("  %-4s [failed]\n", k)); next }
    r2 <- tryCatch(r.squaredGLMM(m), error=function(e) matrix(c(NA,NA),1,2))
    tt <- summary(m)$tTable
    bc <- boot_ci(f, data, B)
    by <- bayes_fit(f, data); bs <- if (!is.null(by)) summary(by)$solutions else NULL
    cat(sprintf("\n  %s : %s\n", k, deparse(f, width.cutoff=200)))
    cat(sprintf("     R2m=%.3f R2c=%.3f  n=%d\n", r2[1,1], r2[1,2], nobs(m)))
    for (i in 2:nrow(tt)) {
      nm <- rownames(tt)[i]
      ci <- if (!is.null(bc) && nm %in% bc$term) bc[bc$term==nm, ] else NULL
      pm <- if (!is.null(bs) && nm %in% rownames(bs)) bs[nm, ] else NULL
      cat(sprintf("       %-16s b=%11.4f p=%.4f", nm, tt[i,"Value"], tt[i,"p-value"]))
      if (!is.null(ci)) cat(sprintf("  boot95=[%.4g, %.4g]", ci$q025, ci$q975))
      if (!is.null(pm)) cat(sprintf("  pMCMC=%.4f CrI=[%.4g, %.4g]", pm["pMCMC"], pm["l-95% CI"], pm["u-95% CI"]))
      cat("\n")
      rows[[length(rows)+1]] <- data.frame(
        variant=vlabel, eq=k, response=all.vars(f)[1], term=nm,
        estimate=tt[i,"Value"], se=tt[i,"Std.Error"], p=tt[i,"p-value"],
        R2m=r2[1,1], R2c=r2[1,2], n=nobs(m),
        boot_q025=if(!is.null(ci)) ci$q025 else NA_real_,
        boot_q975=if(!is.null(ci)) ci$q975 else NA_real_,
        pMCMC=if(!is.null(pm)) unname(pm["pMCMC"]) else NA_real_,
        crI_l=if(!is.null(pm)) unname(pm["l-95% CI"]) else NA_real_,
        crI_u=if(!is.null(pm)) unname(pm["u-95% CI"]) else NA_real_,
        row.names=NULL, stringsAsFactors=FALSE)
    }
  }
  ok <- Filter(function(m) inherits(m,"lme"), models)
  fc <- NULL; std <- NULL; n_claims <- NA
  if (!is.null(dsep)) {
    ps <- sapply(dsep, function(tr) {
      parents <- if (nchar(tr[3])) strsplit(tr[3], ",")[[1]] else character(0)
      p <- dsep_p(tr[1], tr[2], parents, data)
      if (!is.na(p) && p < 0.05) cat(sprintf("    FAIL: %s _|_ %s | %s  p=%.4f\n", tr[1], tr[2], tr[3], p))
      p })
    fc <- fisher_C_hand(ps); n_claims <- if (!is.null(fc)) fc$k else NA
    std <- tryCatch({ m <- do.call(psem, unname(ok)); coefs(m, standardize="scale") }, error=function(e) NULL)
  } else {
    mdl <- tryCatch(do.call(psem, unname(ok)), error=function(e){cat("  psem err:",conditionMessage(e),"\n"); NULL})
    if (!is.null(mdl)) {
      f0 <- tryCatch(fisherC(mdl), error=function(e) NULL)
      ds <- tryCatch(dSep(mdl, .progressBar=FALSE), error=function(e) NULL)
      std <- tryCatch(coefs(mdl, standardize="scale"), error=function(e) NULL)
      nf <- if (!is.null(ds)) sum(ds$P.Value < 0.05, na.rm=TRUE) else NA
      if (!is.null(ds)) { n_claims <- nrow(ds)
        fails <- ds[!is.na(ds$P.Value) & ds$P.Value < 0.05, , drop=FALSE]
        if (nrow(fails)) { cat("  d-sep fails:\n"); print(fails[,c("Independ.Claim","P.Value")], row.names=FALSE) } }
      if (!is.null(f0)) fc <- list(Fisher.C=f0$Fisher.C, df=f0$df, P.Value=f0$P.Value, k=n_claims, n_fails=nf)
    }
  }
  if (!is.null(fc))
    cat(sprintf("\n  Fisher C = %.2f  df = %d  p = %.4f   (%d/%d d-sep fails)\n",
                fc$Fisher.C, fc$df, fc$P.Value, fc$n_fails, fc$k))
  list(rows=do.call(rbind, rows), fisher=fc, std=std, label=vlabel)
}

MACRO_no_zoo <- list(
  M1 = d13C_DIC      ~ TN + Day,
  M2 = d13C_DOC      ~ d13C_DIC + Temp + Day,
  M3 = d13C_CO2_Conc ~ d13C_DIC + TN + Day,
  M4 = d13C_CH4_Conc ~ d13C_DOC + Day,
  M5 = d13C_Macro    ~ d13C_DIC + d13C_CO2_Conc + Day,
  M6 = MOB_I         ~ d13C_DOC + Day,
  M7 = MOB_II        ~ d13C_DIC + Day,
  M8 = Bac           ~ MOB_I + MOB_II + d13C_DOC + Day)
MACRO_zoo <- c(MACRO_no_zoo, list(
  M9 = d13C_Zoo      ~ MOB_II + P18_2w6 + Day))

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

PHYTO_no_zoo <- list(
  P1   = d13C_DIC      ~ PO4_P + Day,
  PDOC = d13C_DOC      ~ d13C_DIC + Day,
  P2   = d13C_CO2_Conc ~ d13C_DIC + MOB_I + Day,
  P3   = d13C_CH4_Conc ~ Chl_a + Oxidized_N + Day,
  P4   = Pa15          ~ MOB_I + Chl_a + Day,
  P5   = MOB_II        ~ d13C_DIC + Day,
  P6   = P18_2w6       ~ d13C_DIC + Chl_a + Day)
PHYTO_zoo <- list(
  P1   = d13C_DIC      ~ PO4_P + Day,
  PDOC = d13C_DOC      ~ d13C_DIC + Day,
  P2   = d13C_CO2_Conc ~ d13C_DIC + MOB_I + Day,
  P3   = d13C_CH4_Conc ~ Chl_a + Oxidized_N + d13C_Zoo + Day,
  P4   = Pa15          ~ MOB_I + Chl_a + Day,
  P5   = MOB_II        ~ d13C_DIC + Day,
  P6   = P18_2w6       ~ d13C_DIC + Chl_a + Day,
  P7   = d13C_Zoo      ~ P18_2w6 + Pa15 + Day)

R1 <- run_variant("MACRO-no-zoo (Macrophyte-dominated, var-specific TH2)", MACRO_no_zoo, macro, dsep=m1_dsep)
R2 <- run_variant("MACRO-zoo (Macrophyte-dominated, var-specific TH2)",    MACRO_zoo,    macro, dsep=m2_dsep)
R3 <- run_variant("PHYTO-no-zoo (Phytoplankton-dominated, +DIC->DOC)",      PHYTO_no_zoo, phyto)
R4 <- run_variant("PHYTO-zoo (Phytoplankton-dominated, +DIC->DOC +Zoo->CH4) = FIGURE", PHYTO_zoo, phyto)

coef_rows <- do.call(rbind, list(R1$rows,R2$rows,R3$rows,R4$rows))
write.csv(coef_rows, CSV, row.names=FALSE)

fisher_rows <- do.call(rbind, lapply(list(R1,R2,R3,R4), function(r) {
  if (is.null(r$fisher)) return(NULL)
  data.frame(variant=r$label, Fisher_C=r$fisher$Fisher.C, df=r$fisher$df, p=r$fisher$P.Value,
             n_claims=r$fisher$k, n_fails=r$fisher$n_fails, row.names=NULL)
}))
write.csv(fisher_rows, FCSV, row.names=FALSE)

std_rows <- do.call(rbind, lapply(list(R1,R2,R3,R4), function(r) {
  if (is.null(r$std)) return(NULL)
  keep <- intersect(c("Response","Predictor","Std.Estimate","Std.Error","P.Value"), names(r$std))
  cbind(variant=r$label, r$std[, keep])
}))
write.csv(std_rows, SCSV, row.names=FALSE)

cat("\n", paste(rep("=",74),collapse=""), "\n TH3 day-1 CO2 sensitivity (MACRO M3)\n",
    paste(rep("=",74),collapse=""), "\n", sep="")
m3_full <- fit_lme(d13C_CO2_Conc ~ d13C_DIC + TN + Day, macro)
macro_s <- macro; macro_s$d13C_CO2_Conc[macro_s$Mesocosm=="TH3" & macro_s$Day==1] <- NA
m3_sens <- fit_lme(d13C_CO2_Conc ~ d13C_DIC + TN + Day, macro_s)
for (tag in c("full","sens")) {
  m <- if (tag=="full") m3_full else m3_sens
  r2 <- tryCatch(r.squaredGLMM(m), error=function(e) matrix(c(NA,NA),1,2)); tt <- summary(m)$tTable
  cat(sprintf("\n  M3 [%s]  R2m=%.3f  n=%d\n", tag, r2[1,1], nobs(m)))
  for (i in 2:nrow(tt)) cat(sprintf("     %-12s b=%11.4f p=%.4f\n", rownames(tt)[i], tt[i,"Value"], tt[i,"p-value"]))
}

cat("\n", paste(rep("=",74),collapse=""), "\n SUMMARY (four variants)\n",
    paste(rep("=",74),collapse=""), "\n", sep="")
print(fisher_rows, row.names=FALSE)
cat("\n COEFS CSV : ", CSV, "\n FISHER CSV: ", FCSV, "\n STDCOEF   : ", SCSV, "\n LOG       : ", LOG, "\n")
sink()
