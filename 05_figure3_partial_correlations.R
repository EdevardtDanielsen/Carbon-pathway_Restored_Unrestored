.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT  <- if (length(.file)) dirname(normalizePath(.file[1])) else getwd()
DATA  <- file.path(ROOT, "data")
FIGS  <- file.path(ROOT, "figures")

library(GGally)
library(ggplot2)
library(cowplot)
library(svglite)
library(ragg)

DATA_CSV <- file.path(DATA, "mesocosm_experiment.csv")
OUTDIR   <- FIGS

cat("=======================================================================\n")
cat("Figure 3: (a) macrophyte over (b) phytoplankton, vertical stack\n")
cat("=======================================================================\n\n")

FAM        <- "Arial"
STRIP_PT   <- 11.0
R_SIZE     <- 2.9
STAR_SIZE  <- 3.1
PT_SIZE    <- 0.40
LINE_LW    <- 0.35
HIST_LW    <- 0.10
COL_POS    <- "#CC3311"
COL_NEG    <- "#0072B2"
LABEL_PT   <- 13.0

dat <- read.csv(DATA_CSV, fileEncoding = "latin1")

rename_col <- function(d, p, n) { m <- grep(p, names(d)); if (length(m) == 1) names(d)[m] <- n; d }
dat <- rename_col(dat, "^.*13C.DIC$", "d13C_DIC")
dat <- rename_col(dat, "^.*13C.DOC$", "d13C_DOC")
dat <- rename_col(dat, "13C.CH4.Conc", "d13C_CH4")
dat <- rename_col(dat, "13C.CO2.Conc", "d13C_CO2")
dat <- rename_col(dat, "^FCH4", "FCH4")
dat <- rename_col(dat, "^FCO2", "FCO2")
dat <- rename_col(dat, "13C.*Macro", "d13C_Macro")
dat <- rename_col(dat, "Zooplankton", "d13C_Zoo")
dat <- rename_col(dat, "16.1.*7c.*PLFA", "PLFA_16_1w7c")
dat <- rename_col(dat, "16.1.*5c.*PLFA", "PLFA_16_1w5c")
dat <- rename_col(dat, "18.1.*7c.*PLFA", "PLFA_18_1w7c")
dat <- rename_col(dat, "18.3.*3.*PLFA", "PLFA_18_3n3")
dat <- rename_col(dat, "18.2.*6.*PLFA", "PLFA_18_2w6")
dat <- rename_col(dat, "i15.0.*PLFA", "PLFA_i15_0")
dat <- rename_col(dat, "a15.0.*PLFA", "PLFA_a15_0")

day_col <- grep("^Day$|^Dag$|sampling", names(dat), ignore.case = TRUE, value = TRUE)
if (length(day_col) == 0) day_col <- names(dat)[1]
dat$Day <- as.numeric(dat[[day_col]])
meso_col <- grep("^Mesocosm$", names(dat), ignore.case = TRUE, value = TRUE)
if (length(meso_col) == 0) meso_col <- names(dat)[4]
dat$Mesocosm <- dat[[meso_col]]

LABELS <- c(
  "d13C_DIC" = "DIC", "d13C_DOC" = "DOC", "d13C_Macro" = "Macro",
  "d13C_CH4" = "CH4", "d13C_CO2" = "CO2", "d13C_Zoo" = "Zoo",
  "PLFA_18_3n3" = "18:3n3", "PLFA_18_2w6" = "18:2w6",
  "PLFA_16_1w7c" = "16:1w7c", "PLFA_16_1w5c" = "16:1w5c",
  "PLFA_18_1w7c" = "18:1w7c", "PLFA_i15_0" = "i15:0", "PLFA_a15_0" = "a15:0"
)

build_macro_resid <- function() {
  dat_all  <- dat[dat$Dominance == "macrophyte" & dat$Label == "yes", ]
  dat_plfa <- dat_all[dat_all$Mesocosm != "TH2", ]

  sem_vars <- c("d13C_DIC","d13C_DOC","d13C_Macro","d13C_CH4","d13C_CO2","d13C_Zoo",
                "PLFA_18_3n3","PLFA_18_2w6","PLFA_16_1w7c","PLFA_16_1w5c",
                "PLFA_18_1w7c","PLFA_i15_0","PLFA_a15_0")
  sem_vars <- sem_vars[sem_vars %in% names(dat_all)]
  plfa_vars <- sem_vars[grepl("^PLFA_", sem_vars)]

  rd <- data.frame(row.names = seq_len(nrow(dat_all)))
  for (v in sem_vars) {
    is_plfa <- v %in% plfa_vars
    src <- if (is_plfa) dat_plfa else dat_all
    x_src <- as.numeric(src[[v]]); d_src <- as.numeric(src$Day)
    ok <- complete.cases(x_src, d_src)
    r_all <- rep(NA_real_, nrow(dat_all))
    if (sum(ok) >= 4) {
      fit <- lm(x_src[ok] ~ d_src[ok])
      if (is_plfa) {
        idx <- which(dat_all$Mesocosm != "TH2")
        r_plfa <- rep(NA_real_, nrow(src)); r_plfa[ok] <- residuals(fit)
        r_all[idx] <- r_plfa
      } else r_all[ok] <- residuals(fit)
    }
    rd[[v]] <- r_all
  }
  names(rd) <- LABELS[names(rd)]
  cat(sprintf("(a) macrophyte: n_all=%d, n_plfa=%d, vars=%d\n",
              nrow(dat_all), nrow(dat_plfa), length(sem_vars)))
  list(resid = rd, df_fun = function(n) max(n - 3, 1))
}

build_phyto_resid <- function() {
  dat_p <- dat[dat$Dominance == "phytoplankton" & dat$Label == "yes", ]

  sem_vars <- c("d13C_DIC","d13C_DOC","d13C_CH4","d13C_CO2","d13C_Zoo",
                "PLFA_18_3n3","PLFA_18_2w6","PLFA_16_1w7c","PLFA_16_1w5c",
                "PLFA_18_1w7c","PLFA_i15_0","PLFA_a15_0")
  sem_vars <- sem_vars[sem_vars %in% names(dat_p)]

  MIN_OBS <- 10
  rd <- data.frame(row.names = seq_len(nrow(dat_p))); kept <- c()
  for (v in sem_vars) {
    x <- as.numeric(dat_p[[v]]); d <- as.numeric(dat_p$Day)
    ok <- complete.cases(x, d)
    if (sum(ok) < MIN_OBS) { cat("  (b) omitted", v, ": n =", sum(ok), "< ", MIN_OBS, "\n"); next }
    r_all <- rep(NA_real_, nrow(dat_p))
    fit <- lm(x[ok] ~ d[ok]); r_all[ok] <- residuals(fit)
    rd[[v]] <- r_all; kept <- c(kept, v)
  }
  names(rd) <- LABELS[names(rd)]
  cat(sprintf("(b) phytoplankton: n=%d, vars=%d (Day-continuous detrend)\n",
              nrow(dat_p), length(kept)))
  list(resid = rd, df_fun = function(n) max(n - 3, 1))
}

make_upper <- function(df_fun, s = 1, y_off = 0) {
  function(data, mapping, ...) {
    x <- eval_data_col(data, mapping$x); y <- eval_data_col(data, mapping$y)
    ok <- complete.cases(x, y); n_ok <- sum(ok)
    if (n_ok < 5) {
      return(ggplot() + theme_void() +
        annotate("text", x = .5, y = .5 + y_off, label = "n/a", size = R_SIZE * s, colour = "grey60"))
    }
    r_val <- cor(x[ok], y[ok])
    df_corr <- df_fun(n_ok)
    t_val <- r_val * sqrt(df_corr / (1 - r_val^2))
    p_val <- 2 * pt(-abs(t_val), df = df_corr)
    st <- ifelse(p_val < .001, "***", ifelse(p_val < .01, "**", ifelse(p_val < .05, "*", "")))
    cl <- ifelse(r_val > 0, COL_POS, COL_NEG)
    r_y <- ifelse(st == "", 0.50, 0.70) + y_off
    r_lab <- sub("^(-?)0[.]", "\\1.", sprintf("%.2f", r_val))
    ggplot() + theme_void() +
      annotate("text", x = .5, y = r_y, label = r_lab,
               size = R_SIZE * s, colour = cl, family = FAM) +
      annotate("text", x = .5, y = 0.20 + y_off, label = st,
               size = STAR_SIZE * s, colour = cl, family = FAM) +
      xlim(0, 1) + ylim(0, 1)
  }
}

make_lower <- function(line_col, s = 1) {
  function(data, mapping, ...) {
    ggplot(data, mapping) +
      geom_point(size = PT_SIZE * s, alpha = 0.7, colour = "grey30") +
      geom_smooth(method = "lm", se = FALSE, colour = line_col,
                  linewidth = LINE_LW * s, na.rm = TRUE) +
      theme_minimal(base_size = STRIP_PT * s, base_family = FAM)
  }
}

make_diag <- function(hist_fill, s = 1) {
  function(data, mapping, ...) {
    ggplot(data, mapping) +
      geom_histogram(bins = 7, fill = hist_fill, colour = "white", linewidth = HIST_LW * s) +
      theme_minimal(base_size = STRIP_PT * s, base_family = FAM)
  }
}

matrix_theme <- function(s = 1) {
  theme_bw(base_size = STRIP_PT, base_family = FAM) +
    theme(
      text             = element_text(family = FAM),
      strip.text.x     = element_text(size = STRIP_PT, angle = 90,
                                       hjust = 0.5, vjust = 0.5,
                                       margin = margin(2, 2, 2, 2, "pt")),
      strip.text.y     = element_text(size = STRIP_PT, angle = 0,
                                       hjust = 0.5, vjust = 0.5,
                                       margin = margin(2, 4, 2, 4, "pt")),
      strip.background = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.title       = element_blank(),
      panel.grid       = element_blank(),
      panel.spacing    = unit(1.0, "pt"),
      plot.margin      = margin(t = 11, r = 2, b = 2, l = 2, unit = "pt")
    )
}

to_plotmath <- function(lab) {
  if (lab == "CH4") return("'CH'[4]*'(aq)'")
  if (lab == "CO2") return("'CO'[2]*'(aq)'")
  m <- regmatches(lab, regexec("^([0-9]+:[0-9]+)[wn]([0-9]+c?)$", lab))[[1]]
  if (length(m) == 3) return(paste0("'", m[2], "ω", m[3], "'"))
  paste0("'", lab, "'")
}

build_matrix <- function(parts, line_col, hist_fill, s = 1, y_off = 0) {
  labs <- vapply(names(parts$resid), to_plotmath, character(1), USE.NAMES = FALSE)
  ggpairs(parts$resid,
          columnLabels = labs, labeller = "label_parsed",
          upper = list(continuous = make_upper(parts$df_fun, s, y_off)),
          lower = list(continuous = make_lower(line_col, s)),
          diag  = list(continuous = make_diag(hist_fill, s))) +
    matrix_theme(s)
}

cat("\nBuilding panels...\n")
macro <- build_macro_resid()
phyto <- build_phyto_resid()

n_a <- ncol(macro$resid); n_b <- ncol(phyto$resid)
SCALE_B <- n_a / n_b
cat(sprintf("Panels: (a) %d vars, (b) %d vars -> panel (b) glyph scale = %.3f\n",
            n_a, n_b, SCALE_B))

B_YOFF <- -0.05
p_a <- build_matrix(macro, line_col = "#0072B2", hist_fill = "#0072B2", s = 1)
p_b <- build_matrix(phyto, line_col = "#009E73", hist_fill = "#009E73", s = SCALE_B, y_off = B_YOFF)

cat("Converting to gtables and combining (vertical stack, square + symmetric)...\n")
g_a <- GGally::ggmatrix_gtable(p_a); g_a$respect <- TRUE
g_b <- GGally::ggmatrix_gtable(p_b); g_b$respect <- TRUE

combined <- plot_grid(
  g_a, g_b, ncol = 1, rel_heights = c(1, 1),
  labels = c("(a)", "(b)"),
  label_size = LABEL_PT, label_fontfamily = FAM, label_fontface = "plain",
  label_x = 0.012, label_y = 0.988, hjust = 0, vjust = 1
)
combined <- combined + theme(plot.margin = margin(t = 8, r = 3, b = 3, l = 6, unit = "pt"))

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
W_MM <- 122; H_MM <- 236

ggsave(file.path(OUTDIR, "Figure_3.pdf"), combined,
       width = W_MM, height = H_MM, units = "mm", device = cairo_pdf)
ggsave(file.path(OUTDIR, "Figure_3.svg"), combined,
       width = W_MM, height = H_MM, units = "mm", device = svglite::svglite, bg = "white")
ggsave(file.path(OUTDIR, "Figure_3.png"), combined,
       width = W_MM, height = H_MM, units = "mm", dpi = 600,
       device = ragg::agg_png, bg = "white")

cat("\nSaved to:", OUTDIR, "\n")
cat("  Figure_3.pdf / .svg / .png  (", W_MM, "x", H_MM, "mm )\n")
cat("=======================================================================\n")
cat("DONE\n")
cat("=======================================================================\n")
