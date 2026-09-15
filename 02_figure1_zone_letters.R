.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT  <- if (length(.file)) dirname(normalizePath(.file[1])) else getwd()
DATA  <- file.path(ROOT, "data")

suppressWarnings(suppressMessages({
  ok <- require(lme4) & require(lmerTest) & require(emmeans)
}))
if(!ok){ cat("MISSING_PACKAGES\n"); quit(status=1) }

compact_letters <- function(means, pmat, alpha=0.05){
  g   <- names(sort(means))
  k   <- length(g)
  conn <- function(a,b) pmat[a,b] >= alpha
  cols <- list(rep(TRUE, k)); names(cols[[1]]) <- g
  for(i in 1:(k-1)) for(j in (i+1):k){
    if(!conn(g[i], g[j])){
      newcols <- list()
      for(col in cols){
        if(col[g[i]] && col[g[j]]){
          c1 <- col; c1[g[j]] <- FALSE
          c2 <- col; c2[g[i]] <- FALSE
          newcols <- c(newcols, list(c1), list(c2))
        } else newcols <- c(newcols, list(col))
      }
      cols <- newcols
    }
  }
  keep <- rep(TRUE, length(cols))
  for(a in seq_along(cols)) for(b in seq_along(cols)) if(a!=b && keep[a] && keep[b]){
    if(all(which(cols[[a]]) %in% which(cols[[b]]))) keep[a] <- FALSE
  }
  cols <- cols[keep]
  lab <- setNames(rep("", k), g)
  for(ci in seq_along(cols)) for(nm in g) if(cols[[ci]][nm]) lab[nm] <- paste0(lab[nm], letters[ci])
  lab[names(means)]
}

path <- file.path(DATA, "field_survey.csv")
d <- read.csv(path, fileEncoding="UTF-8", check.names=FALSE)
names(d) <- trimws(names(d))
fch4 <- grep("FCH", names(d), value=TRUE)[1]
fco2 <- grep("FCO", names(d), value=TRUE)[1]
conv <- 86400*12.011/1e6
d$FCH4 <- as.numeric(d[[fch4]])*conv
d$FCO2 <- as.numeric(d[[fco2]])*conv
d$Position <- as.integer(d[["Position"]])
d$Month <- factor(trimws(as.character(d[["Month"]])))
zmap <- function(p) ifelse(p %in% 1:4,"R1",ifelse(p %in% 9:12,"R2",ifelse(p %in% 5:8,"U1","U2")))
d$Zone <- factor(zmap(d$Position), levels=c("R1","R2","U1","U2"))
d$Position <- factor(d$Position)
d <- d[!is.na(d$Position),]

emit_letters <- function(m, label){
  em   <- emmeans(m, ~Zone)
  emm  <- summary(em)
  means <- setNames(emm$emmean, as.character(emm$Zone))
  pr   <- summary(pairs(em), adjust="tukey")
  zl   <- levels(m@frame$Zone)
  pmat <- matrix(1, length(zl), length(zl), dimnames=list(zl, zl))
  for(i in seq_len(nrow(pr))){
    ab <- strsplit(as.character(pr$contrast[i]), " - ")[[1]]
    a <- trimws(gsub("[()]","",ab[1])); b <- trimws(gsub("[()]","",ab[2]))
    pmat[a,b] <- pmat[b,a] <- pr$p.value[i]
  }
  lab <- compact_letters(means, pmat)
  cat("\n--- ", label, " : emmeans + Tukey compact letters ---\n", sep="")
  out <- data.frame(Zone=names(means), emmean=round(unname(means),4),
                    letters=unname(lab[names(means)]))
  print(out, row.names=FALSE)
  cat("Pairwise (Tukey-adjusted p):\n")
  print(pr)
  invisible(lab)
}

dch <- d[!is.na(d$FCH4) & d$FCH4>0,]
cat("=== CH4  log(FCH4) ~ Zone + (1|Position)+(1|Month)   n =",nrow(dch),
    "(dropped",sum(!is.na(d$FCH4) & d$FCH4<=0),"non-positive)\n")
m1 <- lmer(log(FCH4) ~ Zone + (1|Position)+(1|Month), data=dch,
           control=lmerControl(optimizer="bobyqa"))
cat("Type-III test of Zone:\n"); print(anova(m1))
cat("singular:", isSingular(m1), "\n")
emit_letters(m1, "CH4 (letters on log scale, back-transform monotone)")

dco <- d[!is.na(d$FCO2),]
cat("\n=== CO2  FCO2 ~ Zone + (1|Position)+(1|Month)   n =",nrow(dco),"\n")
m2 <- lmer(FCO2 ~ Zone + (1|Position)+(1|Month), data=dco,
           control=lmerControl(optimizer="bobyqa"))
cat("Type-III test of Zone:\n"); print(anova(m2))
cat("singular:", isSingular(m2), "\n")
emit_letters(m2, "CO2 (raw scale)")
