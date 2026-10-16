## ---------------------------------------------------------------------------
## Extract the inter-subunit interface residue set directly from the
## EXISTING FoldX binding_ddg output (NOT re-running FoldX, per instruction).
## binding_ddg.py tested every residue within 10A of the partner chain, so
## the unique positions appearing in Binding_DDG_{AD,DA}.csv mutation codes
## (format [WT_AA][chain][position][MUT_AA], e.g. "FE196D" = Phe196 on
## chain E -> Asp) ARE the interface residue set by the Xu repo's own
## documented methodology.
## ---------------------------------------------------------------------------
library(dplyr)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "structural"
XU_DIR <- "Xu_KCNA5_Mutations-main/binding_ddg/results"

parse_ddg <- function(path, chain_letter) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  ## Mutation code: [WT_AA][chain][position][MUT_AA]
  m <- regmatches(df$Mutation, regexec(paste0("^([A-Z])", chain_letter, "([0-9]+)([A-Z])$"), df$Mutation))
  df$wt_aa <- sapply(m, `[`, 2)
  df$position <- as.integer(sapply(m, `[`, 3))
  df$mut_aa <- sapply(m, `[`, 4)
  df$chain <- chain_letter
  df
}

ad <- parse_ddg(file.path(XU_DIR, "Binding_DDG_AD.csv"), "E")
da <- parse_ddg(file.path(XU_DIR, "Binding_DDG_DA.csv"), "F")
log_progress("Binding_DDG_AD.csv (chain E): ", nrow(ad), " rows, ", length(unique(ad$position)), " unique positions")
log_progress("Binding_DDG_DA.csv (chain F): ", nrow(da), " rows, ", length(unique(da$position)), " unique positions")

all_ddg <- bind_rows(ad, da)

## Per-residue summary: max |DDG| across all substitutions tested at that
## position (used later for the summary figure's point-size encoding)
residue_ddg <- all_ddg %>%
  group_by(position) %>%
  summarise(n_substitutions = n(), max_abs_ddg = max(abs(DDG), na.rm = TRUE),
            mean_ddg = mean(DDG, na.rm = TRUE), .groups = "drop") %>%
  arrange(position)

interface_residues <- sort(unique(all_ddg$position))
log_progress("Total interface residue positions (union of both chains' FoldX runs): ", length(interface_residues))
log_progress("Interface residues within T1 domain (", T1_DOMAIN_START, "-", T1_DOMAIN_END, "): ",
             sum(interface_residues >= T1_DOMAIN_START & interface_residues <= T1_DOMAIN_END), " of ",
             T1_DOMAIN_END - T1_DOMAIN_START + 1, " T1 positions: ",
             paste(interface_residues[interface_residues >= T1_DOMAIN_START & interface_residues <= T1_DOMAIN_END], collapse = ", "))

save_table(residue_ddg, "interface_residue_ddg_summary", ANALYSIS)
save_table(data.frame(position = interface_residues), "interface_residue_positions", ANALYSIS)

saveRDS(list(interface_residues = interface_residues, residue_ddg = residue_ddg, all_ddg = all_ddg),
        file.path(PROCESSED_DIR, "structural_interface_residues.rds"))
log_progress("Interface residue extraction complete.")
