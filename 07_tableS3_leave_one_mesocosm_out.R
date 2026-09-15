.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT  <- if (length(.file)) dirname(normalizePath(.file[1])) else getwd()
OUTDIR <- file.path(ROOT, "sem_output")
dir.create(OUTDIR, showWarnings = FALSE)

suppressPackageStartupMessages({library(nlme);library(MuMIn);library(piecewiseSEM)})
source(file.path(ROOT, "00_load_sem_data.R"))
OUT <- file.path(OUTDIR, "tableS3_leave_one_mesocosm_out.txt"); sink(OUT)
lmc<-lmeControl(opt="optim",msMaxIter=500)
dat<-load_SEM_data();dat$MOB_II<-dat$P18w7c;dat$Bac<-dat$Pa15
dat$MOB_I<-rowMeans(dat[,c("P16w7c","P16w5c")],na.rm=FALSE)
PT<-c("P16w7c","P16w5c","P18w7c","P18_3n3","P18_2w6","Pi15","Pa15","MOB_I","MOB_II","Bac","d13C_Macro")
z<-function(x)(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
sl<-function(f,d) tryCatch(lme(f,random=~1|Mesocosm,data=d,na.action=na.omit,method="REML",control=lmc),error=function(e)NULL)
key<-function(m,resp,pred){ if(is.null(m))return('fit-fail')
  tt<-tryCatch(summary(m)$tTable,error=function(e)NULL); if(is.null(tt)||!(pred%in%rownames(tt)))return('NA')
  sprintf('%+.3f p=%.3f (n=%d)',tt[pred,'Value'],tt[pred,'p-value'],nobs(m))}

cat('LEAVE-ONE-MESOCOSM-OUT SENSITIVITY (REML refit; FINAL model)\n\n')

mk_macro<-function(drop=NULL){
  m<-dat[dat$Dominance=='macrophyte',]; m[m$Mesocosm=='TH2',PT]<-NA
  if(!is.null(drop)) m<-m[m$Mesocosm!=drop,]; m$Mesocosm<-droplevels(m$Mesocosm); m}
cat('=== MACRO (macrophyte-dominated) — pool/gas paths across mesocosm drops ===\n')
cat(sprintf('%-8s %-22s %-22s %-22s %-22s\n','drop','M1 TN->DIC','M2 DIC->DOC','M3 DIC->CO2','M4 DOC->CH4'))
for(drp in c('none','TH1','TH2','TH3')){
  d<-mk_macro(if(drp=='none')NULL else drp)
  m1<-sl(d13C_DIC~TN+Day,d); m2<-sl(d13C_DOC~d13C_DIC+Temp+Day,d)
  m3<-sl(d13C_CO2_Conc~d13C_DIC+TN+Day,d); m4<-sl(d13C_CH4_Conc~d13C_DOC+Day,d)
  cat(sprintf('%-8s %-22s %-22s %-22s %-22s\n',drp,
    key(m1,'DIC','TN'),key(m2,'DOC','d13C_DIC'),key(m3,'CO2','d13C_DIC'),key(m4,'CH4','d13C_DOC')))
}
cat('\n(PLFA/tissue MACRO equations have only 2 mesocosms to begin with; dropping one leaves a single\n mesocosm and the random intercept is unidentifiable, so LOBO is reported for pool/gas only on the macrophyte-dominated side.)\n')

mk_phyto<-function(drop=NULL){
  p<-dat[dat$Dominance=='phytoplankton',]
  if(!is.null(drop)) p<-p[p$Mesocosm!=drop,]; p$Mesocosm<-droplevels(p$Mesocosm)
  p$Oxidized_N<-rowMeans(cbind(z(p$NO3_N),z(p$NO2_N)),na.rm=FALSE); p}
cat('\n=== PHYTO (phytoplankton-dominated) — key paths + Fisher C across mesocosm drops ===\n')
PZ<-list(P1=d13C_DIC~PO4_P+Day,PDOC=d13C_DOC~d13C_DIC+Day,P2=d13C_CO2_Conc~d13C_DIC+MOB_I+Day,
  P3=d13C_CH4_Conc~Chl_a+Oxidized_N+d13C_Zoo+Day,P4=Pa15~MOB_I+Chl_a+Day,P5=MOB_II~d13C_DIC+Day,
  P6=P18_2w6~d13C_DIC+Chl_a+Day,P7=d13C_Zoo~P18_2w6+Pa15+Day)
for(drp in c('none','URT1','URT2','URT3')){
  p<-mk_phyto(if(drp=='none')NULL else drp)
  ms<-lapply(PZ,sl,d=p)
  fc<-tryCatch({mdl<-do.call(psem,unname(Filter(Negate(is.null),ms))); fisherC(mdl)},error=function(e)NULL)
  cat(sprintf('\n-- drop %s (n=%d) --\n',drp,nrow(p)))
  cat('  P1 PO4->DIC   ',key(ms$P1,'DIC','PO4_P'),'\n')
  cat('  P3 OxN->CH4   ',key(ms$P3,'CH4','Oxidized_N'),'\n')
  cat('  P3 Chl_a->CH4 ',key(ms$P3,'CH4','Chl_a'),'\n')
  cat('  P3 Zoo->CH4   ',key(ms$P3,'CH4','d13C_Zoo'),'\n')
  cat('  P2 DIC->CO2   ',key(ms$P2,'CO2','d13C_DIC'),'\n')
  cat('  PDOC DIC->DOC ',key(ms$PDOC,'DOC','d13C_DIC'),'\n')
  if(!is.null(fc)) cat(sprintf('  Fisher C=%.2f df=%d p=%.4f\n',fc$Fisher.C,fc$df,fc$P.Value))
}
cat('\nDONE\n'); sink()
