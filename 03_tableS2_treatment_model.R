.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT  <- if (length(.file)) dirname(normalizePath(.file[1])) else getwd()
DATA  <- file.path(ROOT, "data")

suppressWarnings(suppressMessages({
  ok <- require(lme4) & require(lmerTest)
}))
if(!ok){ cat("MISSING_PACKAGES\n"); quit(status=1) }

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
restored <- c(1,2,3,4,9,10,11,12)
d$State <- factor(ifelse(d$Position %in% restored,"Restored","Unrestored"),
                      levels=c("Unrestored","Restored"))
zmap <- function(p) ifelse(p %in% 1:4,"R1",ifelse(p %in% 9:12,"R2",ifelse(p %in% 5:8,"U1","U2")))
d$Zone <- factor(zmap(d$Position)); d$Position <- factor(d$Position)
d <- d[!is.na(d$Position),]

cat("rows:", nrow(d), " | positions:", nlevels(d$Position),
    " | zones:", nlevels(d$Zone), " | months:", nlevels(d$Month), "\n")
cat("medians (mg C m-2 d-1): CH4 R/U =",
    round(median(d$FCH4[d$State=="Restored"],na.rm=T),1),"/",
    round(median(d$FCH4[d$State=="Unrestored"],na.rm=T),1),
    " CO2 R/U =",
    round(median(d$FCO2[d$State=="Restored"],na.rm=T),1),"/",
    round(median(d$FCO2[d$State=="Unrestored"],na.rm=T),1),"\n\n")

dch <- d[!is.na(d$FCH4) & d$FCH4>0,]
cat("=== CH4  log(FCH4) ~ State + (1|Zone)+(1|Position)+(1|Month)   n =",nrow(dch),
    "(dropped",sum(!is.na(d$FCH4) & d$FCH4<=0),"non-positive)\n")
m1 <- lmer(log(FCH4) ~ State + (1|Zone)+(1|Position)+(1|Month), data=dch,
           control=lmerControl(optimizer="bobyqa"))
print(round(summary(m1)$coefficients,4)); cat("Variance components:\n"); print(VarCorr(m1))
cat("singular:", isSingular(m1), "\n\n")

dco <- d[!is.na(d$FCO2),]
cat("=== CO2  FCO2 ~ State + (1|Zone)+(1|Position)+(1|Month)   n =",nrow(dco),"\n")
m2 <- lmer(FCO2 ~ State + (1|Zone)+(1|Position)+(1|Month), data=dco,
           control=lmerControl(optimizer="bobyqa"))
print(round(summary(m2)$coefficients,4)); cat("Variance components:\n"); print(VarCorr(m2))
cat("singular:", isSingular(m2), "\n\n")

cat("=== SENSITIVITY: zone means (n = 2 restored, 2 unrestored) ===\n")
zm <- aggregate(cbind(FCH4,FCO2)~Zone+State, d, mean, na.rm=TRUE)
print(zm)
cat("CH4 Welch t on zone means: p =",
    round(t.test(FCH4~State, zm)$p.value,3),
    " | CO2: p =", round(t.test(FCO2~State, zm)$p.value,3), "\n")
