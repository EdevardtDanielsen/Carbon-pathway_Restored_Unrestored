.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grep("^--file=", .args)])
ROOT  <- if (length(.file)) dirname(normalizePath(.file[1])) else getwd()
DATA  <- file.path(ROOT, "data")

SEM_COLUMNS <- list(
  c("^Day$",                    "Day"),
  c("^Mesocosm$",               "Mesocosm"),
  c("^Dominance$",              "Dominance"),
  c("C-CH4 Conc$",              "d13C_CH4_Conc"),
  c("C-CO2 Conc$",              "d13C_CO2_Conc"),
  c("C-DIC$",                   "d13C_DIC"),
  c("C-DOC$",                   "d13C_DOC"),
  c("Zooplankton$",             "d13C_Zoo"),
  c("Macrophytes$",             "d13C_Macro"),
  c("16:1.*7c.*PLFA$",          "P16w7c"),
  c("16:1.*5c.*PLFA$",          "P16w5c"),
  c("18:1.*7c.*PLFA$",          "P18w7c"),
  c("18:3n3 PLFA$",             "P18_3n3"),
  c("18:2.*6 PLFA$",            "P18_2w6"),
  c("i15:0 PLFA$",              "Pi15"),
  c("a15:0 PLFA$",              "Pa15"),
  c("^NO₃",                "NO3_N"),
  c("^NO₂",                "NO2_N"),
  c("^PO₄",                "PO4_P"),
  c("^TN ",                     "TN"),
  c("^Temp",                    "Temp"),
  c("^Chl",                     "Chl_a")
)

NUMERIC_COLUMNS <- c("Day", "NO3_N", "NO2_N", "PO4_P", "TN", "Temp", "Chl_a",
                     "d13C_CH4_Conc", "d13C_CO2_Conc", "d13C_DIC", "d13C_DOC",
                     "d13C_Zoo", "d13C_Macro",
                     "P16w7c", "P16w5c", "P18w7c", "P18_3n3", "P18_2w6", "Pi15", "Pa15")

load_SEM_data <- function() {
  CSV_PATH <- file.path(DATA, "mesocosm_experiment.csv")
  dat <- read.csv(CSV_PATH, fileEncoding = "UTF-8", check.names = FALSE)

  for (spec in SEM_COLUMNS) {
    hits <- grep(spec[1], names(dat), perl = TRUE)
    if (length(hits) != 1L)
      stop(sprintf("'%s' matched %d columns with pattern /%s/ in %s",
                   spec[2], length(hits), spec[1], basename(CSV_PATH)))
    names(dat)[hits] <- spec[2]
  }

  if (!"Label" %in% names(dat)) stop("column 'Label' is missing from ", basename(CSV_PATH))
  dat <- dat[dat$Label == "yes", ]

  for (v in NUMERIC_COLUMNS) dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))

  dat$Mesocosm <- factor(dat$Mesocosm)
  dat
}

if (sys.nframe() == 0) {
  d <- load_SEM_data()
  cat("Load check: n =", nrow(d), "  cols =", ncol(d), "\n")
  for (v in c("Temp", "TN", "NO3_N", "NO2_N", "PO4_P", "Chl_a",
              "d13C_CH4_Conc", "d13C_DIC", "Pa15"))
    cat(sprintf("  %-16s range=[%.4g, %.4g]  n=%d\n", v,
                min(d[[v]], na.rm = TRUE), max(d[[v]], na.rm = TRUE), sum(!is.na(d[[v]]))))
}
