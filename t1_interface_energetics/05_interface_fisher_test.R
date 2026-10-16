## ---------------------------------------------------------------------------
## Analysis 2: interface enrichment Fisher's exact test.
## Primary test: are T1-domain COSMIC mutations enriched at interface
## positions (from the existing FoldX binding_ddg analysis) vs non-interface
## T1 positions?
## Secondary test: restricted to RASA<0.25 (partially buried) positions, to
## confirm the enrichment isn't just "interface residues are bigger/more
## mutable targets" -- depends on 06_rasa_conservation.R's RASA output, so
## that script must be run first despite the file numbering.
## ---------------------------------------------------------------------------
library(dplyr)
source("R/constants.R")
source("R/progress.R")
source("R/plot_theme.R")

ANALYSIS <- "structural"
cosmic <- readRDS(file.path(PROCESSED_DIR, "structural_cosmic_mutations.rds"))
interface <- readRDS(file.path(PROCESSED_DIR, "structural_interface_residues.rds"))
rasa_res <- readRDS(file.path(PROCESSED_DIR, "structural_rasa_conservation_results.rds"))

t1_positions <- T1_DOMAIN_START:T1_DOMAIN_END
mutated_positions <- unique(cosmic$t1_domain$position)
interface_positions <- interface$interface_residues[interface$interface_residues %in% t1_positions]
noninterface_positions <- setdiff(t1_positions, interface_positions)

fisher_2x2 <- function(interface_set, noninterface_set, mutated_set, label) {
  a <- length(intersect(interface_set, mutated_set))       # interface, mutated
  b <- length(interface_set) - a                            # interface, not mutated
  c <- length(intersect(noninterface_set, mutated_set))     # non-interface, mutated
  d <- length(noninterface_set) - c                          # non-interface, not mutated
  tab <- matrix(c(a, c, b, d), nrow = 2,
                dimnames = list(c("Interface", "Non-interface"), c("Mutated", "Not mutated")))
  log_progress(label, " 2x2 table:\n", paste(capture.output(print(tab)), collapse = "\n"))
  ft <- fisher.test(tab)
  log_progress(label, ": OR=", round(unname(ft$estimate), 3),
               " 95%CI=[", round(ft$conf.int[1], 3), ", ", round(ft$conf.int[2], 3), "]",
               " p=", signif(ft$p.value, 3))
  list(table = tab, test = ft)
}

## --- Primary test: all T1-domain positions ---------------------------------
log_progress("--- Primary: interface enrichment, all T1-domain positions ---")
log_progress("T1 domain: ", length(t1_positions), " positions, ", length(interface_positions),
             " interface, ", length(noninterface_positions), " non-interface")
primary <- fisher_2x2(interface_positions, noninterface_positions, mutated_positions, "Primary")

## --- Secondary test: restricted to RASA<0.25 (partially buried) positions --
log_progress("--- Secondary: same test restricted to RASA<0.25 positions ---")
buried <- rasa_res$rasa_t1 %>% filter(RASA < 0.25) %>% pull(position)
interface_buried <- intersect(interface_positions, buried)
noninterface_buried <- intersect(noninterface_positions, buried)
log_progress("Buried (RASA<0.25) T1 positions: ", length(buried), " (", length(interface_buried),
             " interface, ", length(noninterface_buried), " non-interface)")
secondary <- fisher_2x2(interface_buried, noninterface_buried, mutated_positions, "Secondary (RASA<0.25)")

save_table(
  data.frame(
    test = c("Primary (all T1 positions)", "Secondary (RASA<0.25 only)"),
    odds_ratio = c(unname(primary$test$estimate), unname(secondary$test$estimate)),
    ci_low = c(primary$test$conf.int[1], secondary$test$conf.int[1]),
    ci_high = c(primary$test$conf.int[2], secondary$test$conf.int[2]),
    p_value = c(primary$test$p.value, secondary$test$p.value),
    n_interface = c(length(interface_positions), length(interface_buried)),
    n_noninterface = c(length(noninterface_positions), length(noninterface_buried))
  ),
  "interface_enrichment_fisher_tests", ANALYSIS
)

saveRDS(list(primary = primary, secondary = secondary, t1_positions = t1_positions,
             mutated_positions = mutated_positions, interface_positions = interface_positions,
             noninterface_positions = noninterface_positions),
        file.path(PROCESSED_DIR, "structural_interface_fisher_results.rds"))
log_progress("Interface enrichment Fisher's test complete.")
