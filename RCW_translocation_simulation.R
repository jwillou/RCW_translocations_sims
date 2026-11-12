# RED-COCKADED WOODPECKER TRANSLOCATION SIMULATION
#
# MODEL ASSUMPTIONS
# - The ancestral source is large, constant, and initially shared by all donors.
# - Founder gene copies are independent draws from the relevant source pool.
# - All loci are independent, neutral, autosomal, diploid, and biallelic.
# - Generations do not overlap. There is no mutation or selection.
# - TT-to-LP transfers do not remove alleles from TT. Ne is an effective size,
#   not a census count or a limit on the number of birds that can be sampled.
# - No immigration or supplementation occurs after LP establishment.
# - Final samples are independent, unrelated birds under random mating.
#   Use empirical sample sizes AFTER accounting for nestmates/close relatives.
# - This models genetic drift and shared ancestry, not actual donor differences,
#   pedigrees, linkage, individual inbreeding, or sequencing/genotyping errors.
#
# ESTIMATORS
# Hudson: finite-sample correction uses 2 * n_birds (gene copies), NOT n_birds.
# WC84: Weir-Cockerham diploid variance components, including heterozygosity.
# Both estimates are ratios of summed components across loci. Negative estimates
# are retained; neither locus estimates nor final estimates are clipped to zero.
# Confirm which estimator and multilocus aggregation the empirical analysis used.
# References: Bhatia et al. 2013, doi:10.1101/gr.154831.113;
# Weir and Cockerham 1984, doi:10.1111/j.1558-5646.1984.tb05657.x.


#### PARAMETERS ####
random_seed <- 2112
n_replicates <- 1000                
n_loci <- 1000                      # EXAMPLE: loci GENERATED before final filters
n_sample_tt <- 10                   # EXAMPLE: replace with independent TT birds
n_sample_lp <- 10                   # EXAMPLE: replace with independent LP birds
fst_estimator <- "hudson"           # Either "hudson" or "wc84"; both are saved
observed_fst <- NA_real_             # Enter ONLY after estimator/filter matching
# The supplied TT-LP estimate was 0.0547228545; its estimator is not yet verified

# The scenario table uses these values directly; no later block overrides them.
# To expand the grid, replace a value with its commented full-grid option.
founder_models_to_run <- "baseline"  # Full grid: c("baseline", "reduced")
tt_source_fractions <- c(0, 0.5, 1)  # Fraction of LP founders supplied by TT.
Ne_values <- c(10, 25, 50, 100)      # Full grid: c(10, 25, 50, 100)
generations_before_lp <- 2          # Full grid: c(2, 3)
generations_after_lp <- 1           # Full grid: c(1, 2)
frequency_distributions <- "uniform" # Full grid: c("uniform", "beta")

# Founder counts: these are assumed contributors, not estimates of Ne.
# The 12 TT founders include early immigrants, compressed into one episode.
# The 18 LP founders assume all listed 2018-2019 releases contributed equally.
# "reduced" halves both counts. It is a sensitivity assumption, not a finding.
founder_models <- data.frame(
  founder_model = c("baseline", "reduced"),
  founders_tt = c(12, 6),
  founders_lp = c(18, 9),
  stringsAsFactors = FALSE
)

uniform_min <- 0.05
uniform_max <- 0.95
beta_shape1 <- 0.5
beta_shape2 <- 0.5

# Optional filters on the COMBINED final TT + LP sample, not each site separately.
# A locus variable at only one site is kept. A locus fixed for opposite alleles
# in the two sites is also kept. Only loci fixed for the SAME allele are removed
# when the minimum allele count is 1 and the minimum frequency is 0.
# These settings cannot reproduce ascertainment involving other sampled sites.
minimum_pooled_minor_allele_count <- 1
minimum_pooled_minor_allele_frequency <- 0

write_output_files <- TRUE
make_summary_plot <- TRUE
setwd = "/Users/jannawilloughby/Google Drive/My Drive/Willoughby lab/projects - active/red cockaded woodpecker/RCW GBS data analysis/RCW_translocations_sims/"
output_parent <- getwd()            # Outputs go beneath the current R folder.
run_name <- paste0("RCW_simulation_", format(Sys.time(), "%Y%m%d_%H%M%S"))


#### SCENARIO TABLE ####
# Stop early if an input is invalid. seq_len() below also allows zero generations.
stopifnot(
  length(n_replicates) == 1, is.finite(n_replicates),
  n_replicates >= 1, n_replicates == floor(n_replicates),
  length(n_loci) == 1, is.finite(n_loci), n_loci >= 1, n_loci == floor(n_loci),
  length(n_sample_tt) == 1, is.finite(n_sample_tt),
  n_sample_tt >= 2, n_sample_tt == floor(n_sample_tt),
  length(n_sample_lp) == 1, is.finite(n_sample_lp),
  n_sample_lp >= 2, n_sample_lp == floor(n_sample_lp),
  length(fst_estimator) == 1, fst_estimator %in% c("hudson", "wc84"),
  length(observed_fst) == 1,
  is.na(observed_fst) || is.finite(observed_fst),
  length(founder_models_to_run) > 0,
  all(founder_models_to_run %in% founder_models$founder_model),
  all(is.finite(founder_models$founders_tt)),
  all(is.finite(founder_models$founders_lp)),
  all(founder_models$founders_tt >= 1), all(founder_models$founders_lp >= 1),
  all(founder_models$founders_tt == floor(founder_models$founders_tt)),
  all(founder_models$founders_lp == floor(founder_models$founders_lp)),
  length(tt_source_fractions) > 0,
  all(is.finite(tt_source_fractions)),
  all(tt_source_fractions >= 0 & tt_source_fractions <= 1),
  length(Ne_values) > 0, all(is.finite(Ne_values)), all(Ne_values >= 1),
  all(Ne_values == floor(Ne_values)),
  length(generations_before_lp) > 0, all(is.finite(generations_before_lp)),
  all(generations_before_lp >= 0),
  all(generations_before_lp == floor(generations_before_lp)),
  length(generations_after_lp) > 0, all(is.finite(generations_after_lp)),
  all(generations_after_lp >= 0),
  all(generations_after_lp == floor(generations_after_lp)),
  length(frequency_distributions) > 0,
  all(frequency_distributions %in% c("uniform", "beta")),
  uniform_min > 0, uniform_max < 1, uniform_min < uniform_max,
  beta_shape1 > 0, beta_shape2 > 0,
  minimum_pooled_minor_allele_count >= 0,
  minimum_pooled_minor_allele_count == floor(minimum_pooled_minor_allele_count),
  minimum_pooled_minor_allele_count <= n_sample_tt + n_sample_lp,
  minimum_pooled_minor_allele_frequency >= 0,
  minimum_pooled_minor_allele_frequency <= 0.5
)

scenarios <- expand.grid(
  founder_model = unique(founder_models_to_run),
  tt_source_fraction = unique(tt_source_fractions),
  Ne = unique(Ne_values),
  generations_before_lp = unique(generations_before_lp),
  generations_after_lp = unique(generations_after_lp),
  frequency_distribution = unique(frequency_distributions),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

founder_row <- match(scenarios$founder_model, founder_models$founder_model)
scenarios$founders_tt <- founder_models$founders_tt[founder_row]
scenarios$founders_lp <- founder_models$founders_lp[founder_row]
scenarios$scenario_id <- seq_len(nrow(scenarios))

print(scenarios)
cat("\nScenarios:", nrow(scenarios), "\n")
cat("Replicates per scenario:", n_replicates, "\n")
cat("Total replicates:", nrow(scenarios) * n_replicates, "\n")
cat("Example counts must be replaced before an empirical comparison.\n")


#### PREPARE RESULTS ####
set.seed(random_seed)
start_time <- Sys.time()

# One row per replicate; no raw genotype matrices are accumulated in memory.
replicate_results <- data.frame(
  scenario_id = rep(scenarios$scenario_id, each = n_replicates),
  replicate_id = rep(seq_len(n_replicates), times = nrow(scenarios)),
  lp_founders_from_tt = NA_integer_,
  loci_retained = NA_integer_,
  fst_hudson = NA_real_,
  fst_wc84 = NA_real_
)

scenario_summary <- scenarios
scenario_summary$fst_estimator <- fst_estimator
scenario_summary$n_replicates <- n_replicates
scenario_summary$n_valid_fst <- NA_integer_
scenario_summary$median_fst <- NA_real_
scenario_summary$lower_95 <- NA_real_
scenario_summary$upper_95 <- NA_real_
scenario_summary$mean_loci_retained <- NA_real_
scenario_summary$observed_fst <- observed_fst
scenario_summary$observed_within_interval <- NA

# Save chosen settings, including the actual scenario table.
parameters <- list(
  random_seed = random_seed, n_replicates = n_replicates, n_loci = n_loci,
  n_sample_tt = n_sample_tt, n_sample_lp = n_sample_lp,
  fst_estimator = fst_estimator, observed_fst = observed_fst,
  uniform_min = uniform_min, uniform_max = uniform_max,
  beta_shape1 = beta_shape1, beta_shape2 = beta_shape2,
  minimum_pooled_minor_allele_count = minimum_pooled_minor_allele_count,
  minimum_pooled_minor_allele_frequency = minimum_pooled_minor_allele_frequency,
  scenarios = scenarios
)

if (write_output_files) {
  output_directory <- file.path(output_parent, run_name)
  if (dir.exists(output_directory)) {
    stop("That output folder already exists. Change run_name before continuing.")
  }
  dir.create(output_directory, recursive = TRUE)
  write.csv(scenarios, file.path(output_directory, "scenarios.csv"), row.names = FALSE)
  dput(parameters, file = file.path(output_directory, "parameters.txt"))
  writeLines(capture.output(sessionInfo()),
             file.path(output_directory, "R_session_info.txt"))
}

# These values allow the loop body to be inspected manually, line by line.
# Change them to inspect another scenario or another replicate.
scenario_id <- 1
replicate_id <- 1


#### RUN THE SIMULATION ####

for (scenario_id in seq_len(nrow(scenarios))) {
  
  # Read the settings for this one scenario.
  founders_tt <- scenarios$founders_tt[scenario_id]
  founders_lp <- scenarios$founders_lp[scenario_id]
  source_fraction <- scenarios$tt_source_fraction[scenario_id]
  Ne <- scenarios$Ne[scenario_id]
  n_before <- scenarios$generations_before_lp[scenario_id]
  n_after <- scenarios$generations_after_lp[scenario_id]
  source_distribution <- scenarios$frequency_distribution[scenario_id]
  
  for (replicate_id in seq_len(n_replicates)) {
    
    # A. Generate ONE ancestral frequency per locus.
    # The same vector is the source of TT and every non-TT LP founder.
    if (source_distribution == "uniform") {
      p_ancestral <- runif(n_loci, min = uniform_min, max = uniform_max)
    } else {
      p_ancestral <- rbeta(n_loci, shape1 = beta_shape1, shape2 = beta_shape2)
    }
    
    # B. Found TT. Each diploid founder supplies two gene copies at each locus.
    # rbinom() samples the number of alternative-allele copies, not genotypes.
    tt_allele_count <- rbinom(n_loci, size = 2 * founders_tt, prob = p_ancestral)
    p_tt <- tt_allele_count / (2 * founders_tt)
    
    # C. Drift in TT before LP is founded.
    # The next generation has 2*Ne independently sampled gene copies per locus.
    for (generation in seq_len(n_before)) {
      tt_allele_count <- rbinom(n_loci, size = 2 * Ne, prob = p_tt)
      p_tt <- tt_allele_count / (2 * Ne)
    }
    
    # D. Determine how many LP founders come from TT.
    # An integer target is used exactly. A noninteger target is randomly rounded
    # up or down so its average equals the requested fraction.
    # In particular, 50% of 9 means 4 or 5 TT founders, equally likely.
    target_tt_founders <- founders_lp * source_fraction
    lp_founders_from_tt <- floor(target_tt_founders)
    rounding_probability <- target_tt_founders - floor(target_tt_founders)
    if (rounding_probability > 0) {
      lp_founders_from_tt <- lp_founders_from_tt +
        rbinom(1, size = 1, prob = rounding_probability)
    }
    lp_founders_from_pool <- founders_lp - lp_founders_from_tt
    
    # E. Found LP. Draw TT-origin copies at CURRENT TT frequencies and all other
    # copies at the ORIGINAL ancestral frequencies. Sum the counts, then divide
    # by LP's total number of founder gene copies.
    lp_copies_from_tt <- rbinom(n_loci, size = 2 * lp_founders_from_tt, prob = p_tt)
    lp_copies_from_pool <- rbinom(n_loci, size = 2 * lp_founders_from_pool,
                                  prob = p_ancestral)
    lp_allele_count <- lp_copies_from_tt + lp_copies_from_pool
    p_lp <- lp_allele_count / (2 * founders_lp)
    
    # F. Drift independently in TT and LP after LP establishment.
    for (generation in seq_len(n_after)) {
      tt_allele_count <- rbinom(n_loci, size = 2 * Ne, prob = p_tt)
      lp_allele_count <- rbinom(n_loci, size = 2 * Ne, prob = p_lp)
      p_tt <- tt_allele_count / (2 * Ne)
      p_lp <- lp_allele_count / (2 * Ne)
    }
    
    # G. Sample independent birds. Rows are birds; columns are loci.
    # A genotype is 0, 1, or 2 alternative-allele copies. rep(..., each=...)
    # ensures every bird in a column uses that locus's population frequency.
    tt_genotypes <- matrix(
      rbinom(n_sample_tt * n_loci, size = 2,
             prob = rep(p_tt, each = n_sample_tt)),
      nrow = n_sample_tt, ncol = n_loci
    )
    lp_genotypes <- matrix(
      rbinom(n_sample_lp * n_loci, size = 2,
             prob = rep(p_lp, each = n_sample_lp)),
      nrow = n_sample_lp, ncol = n_loci
    )
    
    # H. Obtain sample allele frequencies and observed heterozygosities.
    sample_copies_tt <- colSums(tt_genotypes)
    sample_copies_lp <- colSums(lp_genotypes)
    p_sample_tt <- sample_copies_tt / (2 * n_sample_tt)
    p_sample_lp <- sample_copies_lp / (2 * n_sample_lp)
    h_sample_tt <- colMeans(tt_genotypes == 1)
    h_sample_lp <- colMeans(lp_genotypes == 1)
    
    # I. Optional minor-allele filters, applied to the combined sample.
    total_gene_copies <- 2 * (n_sample_tt + n_sample_lp)
    pooled_alt_count <- sample_copies_tt + sample_copies_lp
    pooled_minor_count <- pmin(pooled_alt_count, total_gene_copies - pooled_alt_count)
    pooled_minor_frequency <- pooled_minor_count / total_gene_copies
    keep_locus <- pooled_minor_count >= minimum_pooled_minor_allele_count &
      pooled_minor_frequency >= minimum_pooled_minor_allele_frequency &
      pooled_minor_count > 0
    loci_retained <- sum(keep_locus)
    
    # Keep the original population frequencies p_tt and p_lp for inspection.
    p1 <- p_sample_tt[keep_locus]
    p2 <- p_sample_lp[keep_locus]
    h1 <- h_sample_tt[keep_locus]
    h2 <- h_sample_lp[keep_locus]
    
    # J. Hudson FST. Correct for sampling GENE COPIES, not number of birds.
    hudson_numerator <- (p1 - p2)^2 -
      p1 * (1 - p1) / (2 * n_sample_tt - 1) -
      p2 * (1 - p2) / (2 * n_sample_lp - 1)
    hudson_denominator <- p1 * (1 - p2) + p2 * (1 - p1)
    fst_hudson <- NA_real_
    if (sum(hudson_denominator) > 0) {
      fst_hudson <- sum(hudson_numerator) / sum(hudson_denominator)
    }
    
    # K. Weir-Cockerham FST for TWO diploid populations (r = 2).
    # n_bar and n_c correct for unequal numbers of sampled birds.
    n_bar <- (n_sample_tt + n_sample_lp) / 2
    n_c <- 2 * n_bar - (n_sample_tt^2 + n_sample_lp^2) / (2 * n_bar)
    p_bar <- (n_sample_tt * p1 + n_sample_lp * p2) / (2 * n_bar)
    h_bar <- (n_sample_tt * h1 + n_sample_lp * h2) / (2 * n_bar)
    s_squared <- (n_sample_tt * (p1 - p_bar)^2 +
                    n_sample_lp * (p2 - p_bar)^2) / n_bar
    
    # Variance among populations (a), among individuals (b), within individuals (c).
    wc_a <- (n_bar / n_c) * (s_squared -
                               (p_bar * (1 - p_bar) - s_squared / 2 - h_bar / 4) / (n_bar - 1))
    wc_b <- (n_bar / (n_bar - 1)) *
      (p_bar * (1 - p_bar) - s_squared / 2 -
         ((2 * n_bar - 1) / (4 * n_bar)) * h_bar)
    wc_c <- h_bar / 2
    fst_wc84 <- NA_real_
    if (sum(wc_a + wc_b + wc_c) > 0) {
      fst_wc84 <- sum(wc_a) / sum(wc_a + wc_b + wc_c)
    }
    
    # L. Store this replicate. Do not replace negative FST estimates with zero.
    result_row <- (scenario_id - 1) * n_replicates + replicate_id
    replicate_results$lp_founders_from_tt[result_row] <- lp_founders_from_tt
    replicate_results$loci_retained[result_row] <- loci_retained
    replicate_results$fst_hudson[result_row] <- fst_hudson
    replicate_results$fst_wc84[result_row] <- fst_wc84
  }
  
  # Summarize this scenario after all its replicates are finished.
  scenario_rows <- ((scenario_id - 1) * n_replicates + 1):(scenario_id * n_replicates)
  if (fst_estimator == "hudson") {
    scenario_fst <- replicate_results$fst_hudson[scenario_rows]
  } else {
    scenario_fst <- replicate_results$fst_wc84[scenario_rows]
  }
  valid_fst <- scenario_fst[is.finite(scenario_fst)]
  scenario_summary$n_valid_fst[scenario_id] <- length(valid_fst)
  scenario_summary$mean_loci_retained[scenario_id] <-
    mean(replicate_results$loci_retained[scenario_rows])
  
  if (length(valid_fst) > 0) {
    scenario_summary$median_fst[scenario_id] <- median(valid_fst)
    scenario_summary$lower_95[scenario_id] <- unname(quantile(valid_fst, 0.025))
    scenario_summary$upper_95[scenario_id] <- unname(quantile(valid_fst, 0.975))
    if (!is.na(observed_fst)) {
      scenario_summary$observed_within_interval[scenario_id] <-
        observed_fst >= scenario_summary$lower_95[scenario_id] &
        observed_fst <= scenario_summary$upper_95[scenario_id]
    }
  }
  
  # Write a checkpoint after each scenario. The raw CSV gets one block appended;
  # summary.csv is refreshed so completed scenarios survive an interrupted run.
  if (write_output_files) {
    write.table(replicate_results[scenario_rows, ],
                file = file.path(output_directory, "replicate_results.csv"),
                sep = ",", row.names = FALSE, col.names = (scenario_id == 1),
                append = (scenario_id != 1), na = "NA")
    write.csv(scenario_summary, file.path(output_directory, "summary.csv"),
              row.names = FALSE)
  }
  
  cat("Completed scenario", scenario_id, "of", nrow(scenarios),
      "| valid estimates:", length(valid_fst),
      "| median FST:", round(scenario_summary$median_fst[scenario_id], 4), "\n")
}


#### SAVE THE RESULTS ####

print(scenario_summary)
elapsed_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat("\nElapsed minutes:", round(elapsed_minutes, 2), "\n")

if (any(scenario_summary$n_valid_fst < n_replicates)) {
  message("Some replicates had no estimable FST. Check n_valid_fst and loci_retained.")
}

if (write_output_files) {
  saveRDS(list(parameters = parameters, scenarios = scenarios,
               summary = scenario_summary, replicates = replicate_results),
          file.path(output_directory, "simulation_results.rds"))
  cat("Output folder:", output_directory, "\n")
}


### PLOT ###
source_fractions <- c(0, 0.5, 1)
source_colors <- c("#59636E", "#0072B2", "#D55E00")
source_symbols <- c(16, 17, 15)
source_lines <- c(3, 1, 1)

# Include all Ne values, holding the other settings constant.
plot_rows <- 
  scenario_summary$founder_model == scenarios$founder_model[1] &
  scenario_summary$generations_before_lp == scenarios$generations_before_lp[1] &
  scenario_summary$generations_after_lp == scenarios$generations_after_lp[1] &
  scenario_summary$frequency_distribution == scenarios$frequency_distribution[1]

plot_data <- scenario_summary[plot_rows, ]
plot_data <- plot_data[is.finite(plot_data$median_fst), ]

plot_ne <- sort(unique(plot_data$Ne))

# Confirm that all three source fractions were simulated.
stopifnot(all(source_fractions %in% plot_data$tt_source_fraction))


### FIGURE 1: GENETIC DIFFERENTIATION ###

if (write_output_files) {
  pdf(file.path(output_directory, "scenario_comparison.pdf"),
      width = 7, height = 5)
}

plot(NA_real_, NA_real_, type = "n",
     xlim = range(plot_ne),
     ylim = range(c(0, plot_data$lower_95, plot_data$upper_95)),
     xaxt = "n",
     xlab = "Assumed effective population size (Ne)",
     ylab = paste0("Multilocus FST (", fst_estimator, ")"),
     main = "Shared ancestry and expected differentiation")

axis(1, at = plot_ne)

for (i in seq_along(source_fractions)) {
  
  source_data <- plot_data[
    plot_data$tt_source_fraction == source_fractions[i], ]
  source_data <- source_data[order(source_data$Ne), ]
  
  # Central 95% simulation intervals.
  segments(source_data$Ne, source_data$lower_95,
           source_data$Ne, source_data$upper_95,
           col = source_colors[i])
  
  # Scenario medians.
  lines(source_data$Ne, source_data$median_fst,
        type = "b", lwd = 2,
        col = source_colors[i],
        pch = source_symbols[i],
        lty = source_lines[i])
}

legend("topright",
       title = "TT source fraction",
       legend = as.character(source_fractions),
       col = source_colors,
       pch = source_symbols,
       lty = source_lines,
       lwd = 2, bty = "n")

if (write_output_files) dev.off()


### CALCULATE REDUCTION RELATIVE TO SOURCE FRACTION 0 ###

effect_data <- plot_data
effect_data$percent_reduction <- NA_real_

for (ne in plot_ne) {
  
  baseline <- effect_data$median_fst[
    effect_data$Ne == ne &
      effect_data$tt_source_fraction == 0
  ]
  
  if (length(baseline) != 1 || baseline <= 0) {
    stop("Each Ne needs one positive median for source fraction 0.")
  }
  
  rows <- effect_data$Ne == ne
  
  effect_data$percent_reduction[rows] <-
    100 * (baseline - effect_data$median_fst[rows]) / baseline
}


#PROPORTIONAL REDUCTION ###

if (write_output_files) {
  pdf(file.path(output_directory, "TT_source_percent_reduction.pdf"),
      width = 7, height = 5)
}

plot(NA_real_, NA_real_, type = "n",
     xlim = range(plot_ne),
     ylim = range(c(0, 100, effect_data$percent_reduction)),
     xaxt = "n",
     xlab = "Assumed effective population size (Ne)",
     ylab = "Reduction in median FST (%)",
     main = "Reduction relative to independent founding")

axis(1, at = plot_ne)

for (i in seq_along(source_fractions)) {
  
  source_data <- effect_data[
    effect_data$tt_source_fraction == source_fractions[i], ]
  source_data <- source_data[order(source_data$Ne), ]
  
  lines(source_data$Ne, source_data$percent_reduction,
        type = "b", lwd = 2,
        col = source_colors[i],
        pch = source_symbols[i],
        lty = source_lines[i])
}

legend("topright",
       title = "TT source fraction",
       legend = as.character(source_fractions),
       col = source_colors,
       pch = source_symbols,
       lty = source_lines,
       lwd = 2, bty = "n")

if (write_output_files) dev.off()

# Plot by Ne
if (make_summary_plot) {
  plot_rows <- scenario_summary$founder_model == scenarios$founder_model[1] &
    scenario_summary$generations_before_lp == scenarios$generations_before_lp[1] &
    scenario_summary$generations_after_lp == scenarios$generations_after_lp[1] &
    scenario_summary$frequency_distribution == scenarios$frequency_distribution[1]
  plot_data <- scenario_summary[plot_rows, ]
  plot_data <- plot_data[is.finite(plot_data$median_fst), ]
  
  if (nrow(plot_data) > 0) {
    plot_ne_values <- sort(unique(plot_data$Ne))
    plot_colors <- rep(c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
                       length.out = length(plot_ne_values))
    plot_symbols <- rep(c(16, 17, 15, 18), length.out = length(plot_ne_values))
    
    # Small horizontal offsets keep overlapping simulation intervals visible.
    # The tick labels still mark the actual TT source fractions: 0, 0.5, and 1.
    plot_offsets <- rep(0, length(plot_ne_values))
    if (length(plot_ne_values) > 1) {
      plot_offsets <- seq(-0.02, 0.02, length.out = length(plot_ne_values))
    }
    
    if (write_output_files) {
      pdf(file.path(output_directory, "scenario_comparison_Ne.pdf"),
          width = 7, height = 5)
    }
    plot_limits <- range(c(0, plot_data$lower_95, plot_data$upper_95,
                           observed_fst), na.rm = TRUE)
    if (diff(plot_limits) == 0) plot_limits <- plot_limits + c(-0.01, 0.01)
    
    # Draw empty axes first, then add one curve per Ne.
    plot(NA_real_, NA_real_, type = "n",
         xlim = c(-0.06, 1.06), ylim = plot_limits, xaxt = "n",
         xlab = "Fraction of LP founders from TT",
         ylab = paste("Multilocus FST (", fst_estimator, ")", sep = ""),
         main = "Shared ancestry and expected differentiation")
    axis(1, at = sort(unique(plot_data$tt_source_fraction)))
    
    for (ne_index in seq_along(plot_ne_values)) {
      ne_plot_data <- plot_data[plot_data$Ne == plot_ne_values[ne_index], ]
      ne_plot_data <- ne_plot_data[order(ne_plot_data$tt_source_fraction), ]
      plot_x <- ne_plot_data$tt_source_fraction + plot_offsets[ne_index]
      
      segments(plot_x, ne_plot_data$lower_95,
               plot_x, ne_plot_data$upper_95, col = plot_colors[ne_index])
      lines(plot_x, ne_plot_data$median_fst, type = "b", lwd = 1.5,
            col = plot_colors[ne_index], pch = plot_symbols[ne_index])
    }
    
    legend_labels <- paste("Ne =", plot_ne_values)
    legend_colors <- plot_colors
    legend_symbols <- plot_symbols
    legend_lines <- rep(1, length(plot_ne_values))
    if (!is.na(observed_fst)) {
      abline(h = observed_fst, lty = 2, col = "black")
      legend_labels <- c(legend_labels, "Observed FST")
      legend_colors <- c(legend_colors, "black")
      legend_symbols <- c(legend_symbols, NA_integer_)
      legend_lines <- c(legend_lines, 2)
    }
    legend("topright", legend = legend_labels, col = legend_colors,
           pch = legend_symbols, lty = legend_lines, bty = "n", cex = 0.8)
    mtext(paste(plot_data$founder_model[1], "|", plot_data$frequency_distribution[1],
                "| generations before/after LP:",
                paste0(plot_data$generations_before_lp[1], "/",
                       plot_data$generations_after_lp[1]),
                "| central 95% simulation intervals"),
          side = 3, line = 0.3, cex = 0.65)
    if (write_output_files) dev.off()
  }
}

# INTERPRETATION
# Compare source fractions WITHIN a demographic/distribution combination.
# Replicates are independent, not paired across source scenarios.
# A value inside a simulation interval is compatible with those assumptions;
# it is not proof of the actual history and is not an estimated migration rate.
# n_loci counts generated loci. Filters and loss of variation can reduce the
# number retained, which is recorded for every replicate. The simulated interval
# assumes independent loci and may be too narrow if empirical loci are linked.
# Before manuscript use, replace demo sample/locus counts, match the empirical
# estimator and ascertainment, and update the Methods to describe what was run.
