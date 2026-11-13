#### PROJECTING CHICK SWAP OUTCOMES ####

#### WD AND DATA ####
setwd("/Users/jannawilloughby/Google Drive/My Drive/Willoughby lab/projects - active/red cockaded woodpecker/RCW GBS data analysis/RCW_translocations_sims/")
panel <- read.csv("RCW_genotype_panel_10000.csv", check.names = FALSE, colClasses = c(sample_id = "character"))
snp_columns <- startsWith(names(panel), "SNP_")
genotypes <- as.matrix(panel[, snp_columns]); rownames(genotypes) <- panel$sample_id

#### PARAMETERS ####
random_seed <- 2112
n_replicates <- 100                  # Replicates for each parameter value set
n_generations <- 20                  # Benefits in offspring can appear from generation 2
years_per_generation <- 5            # A bookkeeping interval; breeding and survival occur ANNUALLY.
n_years <- n_generations * years_per_generation
replicate_seed_stride <- 1000 * (n_generations + 1)
populations <- c("TallTimbers", "LivingstonPlace")
nests_per_year <- c(TallTimbers = 20, LivingstonPlace = 10)  # ENTER ACTUAL COUNT
# The 2025 sheets sampled 6 TT and 4 LP groups
exchange_rates <- c(1, 15, 25, 50)         # TOTAL reciprocal exchanges per five-year window/1 generation, including outside swaps
outside_fractions <- c(0, 0.10, 0.20)  # Sensitivity scenarios
external_populations <- c("Arcadia", "Greenwood")
external_recipient_weights <- c(TallTimbers = 0.2, LivingstonPlace = 0.2)  # Assumed allocation of outside swaps
strategies <- c("no_exchange", "random_exchange", "relatedness_guided")
chicks_per_nest <- 4
adult_survival <- c(M = 0.80, F = 0.75) # Literature-informed scenario values Source/ranges: https://www.fws.gov/species/red-cockaded-woodpecker-dryobates-borealis
eligibility_probability <- 1            # Conditional availability of nestlings for recruitment; upper-bound assumption.
male_probability <- 0.5
avoid_same_rearing_group <- TRUE
pairing_attempts <- 100
n_selection_loci <- 5000
write_output_files <- TRUE
make_plots <- TRUE
output_parent <- getwd()
run_name <- paste0("RCW_annual_projection_", format(Sys.time(), "%Y%m%d_%H%M%S"))
strategy_colors <- c(no_exchange = "#59636E", random_exchange = "#D55E00", relatedness_guided = "#0072B2")
strategy_labels <- c(no_exchange = "No exchange", random_exchange = "Random", relatedness_guided = "Relatedness-guided")
replicate_seed_stride <- 10000 * (n_years + 1)

# define scenarios
scenarios <- data.frame(strategy = "no_exchange", swaps_per_generation = 0, outside_fraction = 0)
for (strategy in strategies[strategies != "no_exchange"]) for (rate in exchange_rates) for (fraction in outside_fractions) {
  scenarios <- rbind(scenarios, data.frame(strategy = strategy, swaps_per_generation = rate, outside_fraction = fraction))
}
print(scenarios)

#### STARTING GENOTYPES AND FIXED EXTERNAL REFERENCE ####
keep <- panel$population %in% c(populations, external_populations) & panel$age_class_recorded == "HY"
initial_G <- genotypes[keep, , drop = FALSE]
initial_meta <- panel[keep, c("sample_id", "population", "sex", "family_group_id")]
names(initial_meta) <- c("id", "population", "sex", "nest")
initial_meta$sex <- ifelse(initial_meta$sex == "Male", "M", ifelse(initial_meta$sex == "Female", "F", NA))
stopifnot(!anyNA(initial_G), all(initial_G %in% 0:2), !anyNA(initial_meta$nest))
stopifnot(n_selection_loci >= 1, n_selection_loci < ncol(initial_G))
focal_rows <- which(initial_meta$population %in% populations)
external_rows <- which(initial_meta$population %in% external_populations)
stopifnot(length(external_rows) > 0, all(populations %in% initial_meta$population))
print(table(initial_meta$population))

# This split is fixed across all alternatives. Marker holdout is not independent-population validation.
set.seed(random_seed)
selection_loci <- sort(sample.int(ncol(initial_G), n_selection_loci))
evaluation_loci <- setdiff(seq_len(ncol(initial_G)), selection_loci)
p_reference <- colMeans(initial_G) / 2
selection_denominator <- 4 * sum(p_reference[selection_loci] * (1 - p_reference[selection_loci]))
evaluation_denominator <- 4 * sum(p_reference[evaluation_loci] * (1 - p_reference[evaluation_loci]))
stopifnot(selection_denominator > 0, evaluation_denominator > 0)
initial_focal_counts <- colSums(initial_G[focal_rows, evaluation_loci, drop = FALSE])
initially_absent_ALT <- initial_focal_counts == 0
initially_absent_REF <- initial_focal_counts == 2 * length(focal_rows)
external_G <- initial_G[external_rows, , drop = FALSE]  # NEVER updated with exports or descendants.
external_meta <- initial_meta[external_rows, , drop = FALSE]
# Each external arrival resamples an empirical chick-genotype proxy from this unchanged pool.
# Proxy genotypes may recur across years; this is a fixed reservoir, not an evolving donor population.
# K below is relative genomic coancestry (half a VanRaden-type relationship), not absolute pedigree F.

result_rows <- list(); swap_rows <- list(); budget_rows <- list(); starting_rows <- list()

#### PROJECTIONS ####
for (replicate_id in seq_len(n_replicates)) {
  seed_base <- random_seed + replicate_id * replicate_seed_stride
  set.seed(seed_base)
  proxy_meta <- initial_meta
  unknown <- which(is.na(proxy_meta$sex))
  proxy_meta$sex[unknown] <- ifelse(runif(length(unknown)) < male_probability, "M", "F")
  donor_meta <- proxy_meta[external_rows, , drop = FALSE]
  
  # Bootstrap observed chick-genotype proxies to the specified census of breeding slots.
  # Copies receive unique IDs but keep the source rearing group; no new initial alleles are invented.
  # This expansion is an explicit assumption, not a reconstruction of unsampled adult genotypes.
  starting_indexes <- integer()
  for (pop_index in seq_along(populations)) {
    site <- populations[pop_index]; set.seed(seed_base + pop_index * 10)
    for (sex in c("M", "F")) {
      candidates <- which(proxy_meta$population == site & proxy_meta$sex == sex)
      if (length(candidates) == 0) stop("No starting ", sex, " genotype proxies for ", site, "; review sexes or pool definition.")
      selected <- candidates[sample.int(length(candidates), nests_per_year[site], replace = TRUE)]
      starting_indexes <- c(starting_indexes, selected)
    }
  }
  starting_G <- initial_G[starting_indexes, , drop = FALSE]
  starting_meta <- proxy_meta[starting_indexes, , drop = FALSE]
  starting_rows[[replicate_id]] <- data.frame(replicate = replicate_id, individual = seq_along(starting_indexes),
                                              source_sample_id = starting_meta$id, population = starting_meta$population, sex = starting_meta$sex,
                                              source_rearing_group = starting_meta$nest, row.names = NULL)
  starting_meta$id <- paste0("r", replicate_id, "_initial_", seq_len(nrow(starting_meta)))
  starting_meta$external_ancestry <- 0; starting_meta$birth_year <- 0
  
  for (scenario_id in seq_len(nrow(scenarios))) {
    strategy <- scenarios$strategy[scenario_id]
    rate <- scenarios$swaps_per_generation[scenario_id]; fraction <- scenarios$outside_fraction[scenario_id]
    parent_G <- starting_G; parent_meta <- starting_meta
    # Combine rows after each scenario, rather than retaining millions of tiny data frames.
    scenario_result_rows <- list(); scenario_budget_rows <- list(); scenario_swap_rows <- list()
    
    # The random and guided policies use exactly the same requested calendar.
    # 15/5 = 3 swaps/year; 25/5 = 5/year. A remainder is placed in randomly selected years.
    # An external swap REPLACES an internal swap in the total budget; it is not an additional release.
    calendar <- data.frame(year = integer(), block = integer(), external = logical(), recipient = character())
    if (rate > 0) {
      rate_index <- match(rate, exchange_rates); fraction_index <- match(fraction, outside_fractions)
      set.seed(seed_base + 100 + rate_index * 10 + fraction_index)
      for (block in seq_len(n_generations)) {
        local_years <- rep(seq_len(years_per_generation), each = floor(rate / years_per_generation))
        remainder <- rate %% years_per_generation
        if (remainder > 0) local_years <- c(local_years, sample.int(years_per_generation, remainder))
        local_years <- sort(local_years)
        is_external <- runif(rate) < fraction
        recipient <- sample(populations, rate, replace = TRUE, prob = external_recipient_weights[populations])
        calendar <- rbind(calendar, data.frame(year = (block - 1) * years_per_generation + local_years,
                                               block = block, external = is_external, recipient = recipient))
      }
    }
    
    for (year in seq_len(n_years)) {
      block <- ceiling(year / years_per_generation)
      pairs <- data.frame(population = character(), father = integer(), mother = integer())
      for (pop_index in seq_along(populations)) {
        site <- populations[pop_index]; set.seed(seed_base + year * 10000 + pop_index * 100 + 1)
        males <- which(parent_meta$population == site & parent_meta$sex == "M")
        females <- which(parent_meta$population == site & parent_meta$sex == "F")
        n_pairs <- min(nests_per_year[site], length(males), length(females)); paired <- FALSE
        while (n_pairs > 0 && !paired) {
          for (attempt in seq_len(pairing_attempts)) {
            fathers <- males[sample.int(length(males), n_pairs)]; mothers <- females[sample.int(length(females), n_pairs)]
            if (!avoid_same_rearing_group || all(parent_meta$nest[fathers] != parent_meta$nest[mothers])) { paired <- TRUE; break }
          }
          if (!paired) n_pairs <- n_pairs - 1
        }
        if (n_pairs > 0) pairs <- rbind(pairs, data.frame(population = site, father = fathers, mother = mothers))
      }
      
      n_chicks <- nrow(pairs) * chicks_per_nest
      child_G <- matrix(NA_integer_, nrow = n_chicks, ncol = ncol(initial_G))
      natal_nests <- character(n_chicks); child_ids <- character(n_chicks)
      if (n_chicks > 0) {
        natal_nests <- rep(paste0("r", replicate_id, "_y", year, "_nest", seq_len(nrow(pairs))), each = chicks_per_nest)
        child_ids <- paste0("r", replicate_id, "_y", year, "_chick", seq_len(n_chicks))
      }
      children <- data.frame(id = child_ids, population = rep(pairs$population, each = chicks_per_nest),
                             sex = rep(NA_character_, n_chicks), nest = natal_nests, natal_nest = natal_nests,
                             moved = rep(FALSE, n_chicks), external = rep(FALSE, n_chicks), external_ancestry = rep(0, n_chicks),
                             birth_year = rep(year, n_chicks))
      birth_population <- children$population
      expected_H <- mate_K <- born_H <- setNames(rep(NA_real_, length(populations)), populations)
      for (pop_index in seq_along(populations)) {
        site <- populations[pop_index]; set.seed(seed_base + year * 10000 + pop_index * 100 + 2)
        pair_rows <- which(pairs$population == site); child_rows <- which(birth_population == site)
        if (length(child_rows) == 0) next
        fathers <- rep(pairs$father[pair_rows], each = chicks_per_nest); mothers <- rep(pairs$mother[pair_rows], each = chicks_per_nest)
        paternal <- parent_G[fathers, , drop = FALSE] / 2; maternal <- parent_G[mothers, , drop = FALSE] / 2
        child_G[child_rows, ] <- matrix(rbinom(length(paternal), 1, paternal), nrow = length(child_rows)) +
          matrix(rbinom(length(maternal), 1, maternal), nrow = length(child_rows))
        children$sex[child_rows] <- ifelse(runif(length(child_rows)) < male_probability, "M", "F")
        children$external_ancestry[child_rows] <- (parent_meta$external_ancestry[fathers] + parent_meta$external_ancestry[mothers]) / 2
        gm <- parent_G[pairs$father[pair_rows], evaluation_loci, drop = FALSE] / 2
        gf <- parent_G[pairs$mother[pair_rows], evaluation_loci, drop = FALSE] / 2
        expected_H[site] <- mean(gm * (1 - gf) + gf * (1 - gm))
        xm <- sweep(2 * gm, 2, 2 * p_reference[evaluation_loci], "-")
        xf <- sweep(2 * gf, 2, 2 * p_reference[evaluation_loci], "-")
        mate_K[site] <- mean(rowSums(xm * xf) / evaluation_denominator)
        born_H[site] <- mean(child_G[child_rows, evaluation_loci, drop = FALSE] == 1)
      }
      
      events <- which(calendar$year == year)
      completed_internal <- completed_external <- 0L
      requested_internal <- sum(!calendar$external[events]); requested_external <- sum(calendar$external[events])
      imported <- external_imported <- setNames(integer(length(populations)), populations)
      used_nests <- character()
      donor_nests <- paste0("external_y", year, "_", donor_meta$nest)
      if (strategy == "relatedness_guided" && n_chicks > 0 && length(events) > 0) {
        X <- sweep(rbind(child_G[, selection_loci, drop = FALSE], external_G[, selection_loci, drop = FALSE]),
                   2, 2 * p_reference[selection_loci], "-")
        selection_K <- tcrossprod(X) / selection_denominator
      }
      set.seed(seed_base + year * 10000 + 9000)
      for (event in events) {
        is_external <- calendar$external[event]
        if (is_external) {
          site_a <- calendar$recipient[event]; site_b <- "External"
          a <- which(children$population == site_a & !children$moved & !children$natal_nest %in% used_nests)
          b <- n_chicks + which(!donor_nests %in% used_nests)
        } else {
          site_a <- populations[1]; site_b <- populations[2]
          a <- which(children$population == site_a & !children$moved & !children$natal_nest %in% used_nests)
          b <- which(children$population == site_b & !children$moved & !children$natal_nest %in% used_nests)
        }
        candidates <- expand.grid(a = a, b = b)
        all_sexes <- c(children$sex, donor_meta$sex)
        candidates <- candidates[all_sexes[candidates$a] == all_sexes[candidates$b], , drop = FALSE]
        if (nrow(candidates) == 0) next
        choice <- sample.int(nrow(candidates), 1)
        
        if (strategy == "relatedness_guided") {
          affected_sites <- if (is_external) site_a else populations
          site_scores <- matrix(NA_real_, nrow(candidates), length(affected_sites))
          # Exact change in the mean potential-mate score. A foster swap replaces a same-sex
          # occupant of a rearing group, so the number of permitted opposite-sex pairings is unchanged.
          # This computes the same score as rebuilding every trial population, with fewer operations.
          for (s in seq_along(affected_sites)) {
            site <- affected_sites[s]
            m <- which(children$population == site & children$sex == "M"); f <- which(children$population == site & children$sex == "F")
            allowed <- matrix(TRUE, length(m), length(f))
            if (avoid_same_rearing_group) allowed <- outer(children$nest[m], children$nest[f], "!=")
            denominator <- sum(allowed)
            if (denominator == 0) next
            base_score <- mean(selection_K[m, f, drop = FALSE][allowed])
            slots <- if (site == site_a) candidates$a else candidates$b
            incoming <- if (site == site_a) candidates$b else candidates$a
            for (sex in c("M", "F")) {
              rows <- which(children$sex[slots] == sex)
              opposite <- if (sex == "M") f else m
              if (length(rows) == 0) next
              allowed_rows <- matrix(TRUE, length(rows), length(opposite))
              if (avoid_same_rearing_group) allowed_rows <- outer(children$nest[slots[rows]], children$nest[opposite], "!=")
              old_sum <- rowSums(selection_K[slots[rows], opposite, drop = FALSE] * allowed_rows)
              new_sum <- rowSums(selection_K[incoming[rows], opposite, drop = FALSE] * allowed_rows)
              site_scores[rows, s] <- base_score + (new_sum - old_sum) / denominator
            }
          }
          scores <- rowMeans(site_scores); scores[!is.finite(scores)] <- Inf
          if (any(is.finite(scores))) { best <- which(scores <= min(scores) + 1e-12); choice <- best[sample.int(length(best), 1)] }
        }
        
        i <- candidates$a[choice]; j <- candidates$b[choice]
        outgoing_id <- children$id[i]; natal_a <- children$natal_nest[i]
        if (is_external) {
          d <- j - n_chicks
          incoming_id <- paste0("r", replicate_id, "_y", year, "_external", event)
          natal_b <- donor_nests[d]; donor_source <- donor_meta$population[d]; donor_proxy <- donor_meta$id[d]
          child_G[i, ] <- external_G[d, ]
          children$id[i] <- incoming_id; children$natal_nest[i] <- natal_b
          children$moved[i] <- TRUE; children$external[i] <- TRUE; children$external_ancestry[i] <- 1
          if (strategy == "relatedness_guided") {
            selection_K[i, ] <- selection_K[j, ]; selection_K[, i] <- selection_K[, j]; selection_K[i, i] <- selection_K[j, j]
          }
          imported[site_a] <- imported[site_a] + 1L; external_imported[site_a] <- external_imported[site_a] + 1L
          completed_external <- completed_external + 1L
        } else {
          incoming_id <- children$id[j]; natal_b <- children$natal_nest[j]; donor_source <- site_b; donor_proxy <- NA_character_
          children$population[c(i, j)] <- c(site_b, site_a); children$nest[c(i, j)] <- rev(children$nest[c(i, j)])
          children$moved[c(i, j)] <- TRUE
          imported[c(site_a, site_b)] <- imported[c(site_a, site_b)] + 1L
          completed_internal <- completed_internal + 1L
        }
        used_nests <- c(used_nests, natal_a, natal_b)
        scenario_swap_rows[[length(scenario_swap_rows) + 1L]] <- data.frame(replicate = replicate_id, strategy = strategy, swaps_per_generation = rate,
                                                                            outside_fraction = fraction, year = year, generation = block, exchange_type = if (is_external) "external" else "TT_LP",
                                                                            site_1 = site_a, site_2 = donor_source, chick_1 = outgoing_id, chick_2 = incoming_id,
                                                                            natal_nest_1 = natal_a, natal_nest_2 = natal_b, sex = children$sex[i], external_proxy_id = donor_proxy, row.names = NULL)
      }
      scenario_budget_rows[[length(scenario_budget_rows) + 1L]] <- data.frame(replicate = replicate_id, strategy = strategy, swaps_per_generation = rate,
                                                                              outside_fraction = fraction, year = year, generation = block, requested_internal = requested_internal,
                                                                              requested_external = requested_external, completed_internal = completed_internal, completed_external = completed_external,
                                                                              requested_total = requested_internal + requested_external, completed_total = completed_internal + completed_external)
      
      # Adults survive between breeding seasons. Yearlings can fill sex-specific vacancies the following year.
      # Unrecruited juveniles are not retained as helpers/floaters; pair bonds are not preserved between years.
      survivor_rows <- recruited_rows <- integer()
      for (pop_index in seq_along(populations)) {
        site <- populations[pop_index]
        set.seed(seed_base + year * 10000 + pop_index * 100 + 3)
        adults <- which(parent_meta$population == site)
        survivors <- adults[runif(length(adults)) < adult_survival[parent_meta$sex[adults]]]
        set.seed(seed_base + year * 10000 + pop_index * 100 + 4)
        resident <- which(children$population == site)
        eligible <- resident[runif(length(resident)) < eligibility_probability]
        set.seed(seed_base + year * 10000 + pop_index * 100 + 5)
        selected <- integer()
        for (sex in c("M", "F")) {
          candidates <- eligible[children$sex[eligible] == sex]
          vacancies <- max(0, nests_per_year[site] - sum(parent_meta$sex[survivors] == sex))
          n_to_recruit <- min(vacancies, length(candidates))
          if (n_to_recruit > 0) selected <- c(selected, candidates[sample.int(length(candidates), n_to_recruit)])
        }
        survivor_rows <- c(survivor_rows, survivors); recruited_rows <- c(recruited_rows, selected)
        next_G <- rbind(parent_G[survivors, evaluation_loci, drop = FALSE], child_G[selected, evaluation_loci, drop = FALSE])
        next_H <- variable_loci <- novel_alleles <- mean_external_ancestry <- NA_real_
        if (nrow(next_G) > 0) {
          next_H <- mean(next_G == 1); counts <- colSums(next_G)
          variable_loci <- sum(counts > 0 & counts < 2 * nrow(next_G))
          novel_alleles <- sum(initially_absent_ALT & counts > 0) + sum(initially_absent_REF & counts < 2 * nrow(next_G))
          mean_external_ancestry <- mean(c(parent_meta$external_ancestry[survivors], children$external_ancestry[selected]))
        }
        scenario_result_rows[[length(scenario_result_rows) + 1L]] <- data.frame(replicate = replicate_id, strategy = strategy,
                                                                                swaps_per_generation = rate, outside_fraction = fraction, year = year, generation = block, population = site,
                                                                                nest_capacity = unname(nests_per_year[site]), breeding_pairs = sum(pairs$population == site),
                                                                                breeding_present = as.integer(any(pairs$population == site)), expected_offspring_H = unname(expected_H[site]),
                                                                                mate_kinship = unname(mate_K[site]), observed_offspring_H = unname(born_H[site]), chicks_born = sum(birth_population == site),
                                                                                imported_chicks = unname(imported[site]), external_imported_chicks = unname(external_imported[site]),
                                                                                surviving_adults = length(survivors), recruited_chicks = length(selected), moved_next_breeders = sum(children$moved[selected]),
                                                                                external_next_breeders = sum(children$external[selected]), next_breeders = nrow(next_G), next_breeder_H = next_H,
                                                                                polymorphic_evaluation_loci = variable_loci, new_alleles_retained = novel_alleles, mean_external_ancestry = mean_external_ancestry)
      }
      parent_G <- rbind(parent_G[survivor_rows, , drop = FALSE], child_G[recruited_rows, , drop = FALSE])
      parent_meta <- rbind(parent_meta[survivor_rows, , drop = FALSE],
                           children[recruited_rows, c("id", "population", "sex", "nest", "external_ancestry", "birth_year"), drop = FALSE])
    }
    result_rows[[length(result_rows) + 1L]] <- do.call(rbind, scenario_result_rows)
    budget_rows[[length(budget_rows) + 1L]] <- do.call(rbind, scenario_budget_rows)
    if (length(scenario_swap_rows) > 0) swap_rows[[length(swap_rows) + 1L]] <- do.call(rbind, scenario_swap_rows)
  }
  cat("Completed replicate", replicate_id, "of", n_replicates, "\n")
}
results <- do.call(rbind, result_rows); annual_budgets <- do.call(rbind, budget_rows); starting_proxies <- do.call(rbind, starting_rows)
swap_log <- if (length(swap_rows) > 0) do.call(rbind, swap_rows) else data.frame()
rm(result_rows, budget_rows, starting_rows, swap_rows)

#### PAIRED COMPARISONS AND SUMMARIES ####
controls <- results[results$strategy == "no_exchange", c("replicate", "year", "population", "expected_offspring_H", "mate_kinship", "next_breeder_H")]
names(controls)[4:6] <- c("control_H", "control_kinship", "control_breeder_H")
comparisons <- merge(results, controls, by = c("replicate", "year", "population"), sort = FALSE)
comparisons$delta_H <- comparisons$expected_offspring_H - comparisons$control_H
comparisons$delta_kinship <- comparisons$mate_kinship - comparisons$control_kinship
comparisons$delta_breeder_H <- comparisons$next_breeder_H - comparisons$control_breeder_H
metrics <- c("expected_offspring_H", "mate_kinship", "delta_H", "delta_kinship", "delta_breeder_H", "next_breeders",
             "breeding_present", "polymorphic_evaluation_loci", "new_alleles_retained", "imported_chicks", "external_imported_chicks",
             "moved_next_breeders", "external_next_breeders", "mean_external_ancestry")
group_columns <- c("strategy", "swaps_per_generation", "outside_fraction", "population", "year", "generation")
group_keys <- do.call(paste, c(comparisons[group_columns], sep = "|"))
groups <- split(seq_len(nrow(comparisons)), group_keys)
summary_rows <- vector("list", length(groups))
for (i in seq_along(groups)) {
  rows <- groups[[i]]; a <- comparisons[rows[1], group_columns, drop = FALSE]
  for (metric in metrics) {
    values <- comparisons[[metric]][rows]; values <- values[is.finite(values)]
    a[[paste0(metric, "_n_valid")]] <- length(values)
    a[[paste0(metric, "_mean")]] <- if (length(values) > 0) mean(values) else NA_real_
    a[[paste0(metric, "_median")]] <- if (length(values) > 0) median(values) else NA_real_
    a[[paste0(metric, "_lower_95")]] <- if (length(values) > 0) unname(quantile(values, 0.025)) else NA_real_
    a[[paste0(metric, "_upper_95")]] <- if (length(values) > 0) unname(quantile(values, 0.975)) else NA_real_
  }
  summary_rows[[i]] <- a
}
summary_results <- do.call(rbind, summary_rows); rownames(summary_results) <- NULL
endpoint_summary <- summary_results[summary_results$year == n_years, ]
endpoint_columns <- c("strategy", "swaps_per_generation", "outside_fraction", "population", "year", "delta_H_mean", "delta_H_lower_95",
                      "delta_H_upper_95", "delta_H_n_valid", "delta_kinship_mean", "delta_breeder_H_mean", "new_alleles_retained_mean",
                      "mean_external_ancestry_mean", "breeding_present_mean", "next_breeders_mean")
print(endpoint_summary[, endpoint_columns], row.names = FALSE)

# Compare realized exchanges with their requested budgets. Unfilled requests are not silently carried forward.
block_budgets <- aggregate(annual_budgets[, c("requested_total", "completed_total", "requested_external", "completed_external")],
                           annual_budgets[, c("replicate", "strategy", "swaps_per_generation", "outside_fraction", "generation")], sum)
block_budgets$completion_fraction <- ifelse(block_budgets$requested_total > 0,
                                            block_budgets$completed_total / block_budgets$requested_total, NA_real_)
block_budgets$realized_external_fraction <- ifelse(block_budgets$completed_total > 0,
                                                   block_budgets$completed_external / block_budgets$completed_total, NA_real_)
budget_keys <- do.call(paste, c(block_budgets[c("strategy", "swaps_per_generation", "outside_fraction")], sep = "|"))
budget_groups <- split(seq_len(nrow(block_budgets)), budget_keys); budget_summary_rows <- vector("list", length(budget_groups))
for (i in seq_along(budget_groups)) {
  rows <- budget_groups[[i]]; a <- block_budgets[rows[1], c("strategy", "swaps_per_generation", "outside_fraction"), drop = FALSE]
  a$requested_per_window <- mean(block_budgets$requested_total[rows]); a$completed_per_window <- mean(block_budgets$completed_total[rows])
  a$completed_external_per_window <- mean(block_budgets$completed_external[rows])
  a$completion_fraction <- if (sum(block_budgets$requested_total[rows]) > 0) sum(block_budgets$completed_total[rows]) / sum(block_budgets$requested_total[rows]) else NA_real_
  a$realized_external_fraction <- if (sum(block_budgets$completed_total[rows]) > 0) sum(block_budgets$completed_external[rows]) / sum(block_budgets$completed_total[rows]) else NA_real_
  budget_summary_rows[[i]] <- a
}
budget_summary <- do.call(rbind, budget_summary_rows); rownames(budget_summary) <- NULL
print(budget_summary, row.names = FALSE)

#### SAVE RESULTS ####
parameters <- list(random_seed = random_seed, n_replicates = n_replicates, n_generations = n_generations, years_per_generation = years_per_generation,
                   n_years = n_years, nests_per_year = nests_per_year, exchange_rates = exchange_rates, outside_fractions = outside_fractions, scenarios = scenarios,
                   external_populations = external_populations, external_recipient_weights = external_recipient_weights, chicks_per_nest = chicks_per_nest,
                   adult_survival = adult_survival, eligibility_probability = eligibility_probability, male_probability = male_probability,
                   avoid_same_rearing_group = avoid_same_rearing_group, pairing_attempts = pairing_attempts,
                   selection_locus_names = colnames(initial_G)[selection_loci], evaluation_locus_names = colnames(initial_G)[evaluation_loci],
                   p_reference = p_reference, initial_metadata = initial_meta, external_proxy_ids = external_meta$id)
if (write_output_files) {
  output_directory <- file.path(output_parent, run_name)
  if (dir.exists(output_directory)) stop("Output folder exists; change run_name before saving.")
  dir.create(output_directory, recursive = TRUE)
  write.csv(comparisons, file.path(output_directory, "projection_replicates.csv"), row.names = FALSE)
  write.csv(summary_results, file.path(output_directory, "projection_summary.csv"), row.names = FALSE)
  write.csv(endpoint_summary[, endpoint_columns], file.path(output_directory, "endpoint_summary.csv"), row.names = FALSE)
  write.csv(annual_budgets, file.path(output_directory, "annual_exchange_counts.csv"), row.names = FALSE)
  write.csv(block_budgets, file.path(output_directory, "five_year_exchange_counts.csv"), row.names = FALSE)
  write.csv(budget_summary, file.path(output_directory, "exchange_budget_summary.csv"), row.names = FALSE)
  write.csv(swap_log, file.path(output_directory, "simulated_swaps.csv"), row.names = FALSE)
  write.csv(starting_proxies, file.path(output_directory, "starting_genotype_proxies.csv"), row.names = FALSE)
  write.csv(scenarios, file.path(output_directory, "exchange_scenarios.csv"), row.names = FALSE)
  saveRDS(list(parameters = parameters, results = comparisons, summary = summary_results, endpoints = endpoint_summary,
               swaps = swap_log, annual_budgets = annual_budgets, block_budgets = block_budgets, budget_summary = budget_summary,
               starting_proxies = starting_proxies), file.path(output_directory, "projection_results.rds"))
  writeLines(capture.output(sessionInfo()), file.path(output_directory, "R_session_info.txt"))
}

#### TRAJECTORIES WITH SHADED 95% SIMULATION BANDS ####
# These are replicate simulation intervals, not confidence intervals for the mean.
# Gaps in estimable values break the bands; extinct/nonbreeding replicates are not assigned H = 0.
plot_metrics <- c("expected_offspring_H", "delta_H", "mate_kinship", "new_alleles_retained", "breeding_present")
plot_labels <- c("Expected offspring heterozygosity", "Change in offspring heterozygosity",
                 "Mean mate coancestry", "Retained alleles new to focal samples", "Replicates with breeding pairs")
if (make_plots) {
  if (write_output_files) pdf(file.path(output_directory, "projection_trajectories.pdf"), width = 12, height = 8.5)
  old_par <- par(no.readonly = TRUE)
  for (metric_index in seq_along(plot_metrics)) for (fraction in outside_fractions) {
    metric <- plot_metrics[metric_index]
    par(mfrow = c(length(populations), length(exchange_rates)), mar = c(4.6, 4.6, 3.2, 1), oma = c(1, 1, 3, 0))
    plot_values <- c(summary_results[[paste0(metric, "_lower_95")]], summary_results[[paste0(metric, "_upper_95")]])
    plot_values <- plot_values[is.finite(plot_values)]
    limits <- if (length(plot_values) > 0) range(plot_values) else c(0, 1)
    if (diff(limits) == 0) limits <- limits + c(-0.001, 0.001)
    for (site in populations) for (rate in exchange_rates) {
      d <- summary_results[summary_results$population == site &
                             ((summary_results$swaps_per_generation == rate & summary_results$outside_fraction == fraction) | summary_results$strategy == "no_exchange"), ]
      plot(NA_real_, NA_real_, xlim = c(1, max(2, n_years)), ylim = limits, xlab = "Year", ylab = plot_labels[metric_index],
           main = paste0(site, ": ", rate, " swaps / ", years_per_generation, " years"), cex.main = 0.9)
      if (metric == "delta_H") abline(h = 0, col = "grey70", lty = 3)
      for (strategy in strategies) {
        a <- d[d$strategy == strategy, ]; a <- a[order(a$year), ]
        center <- a[[paste0(metric, "_median")]]
        if (metric == "breeding_present") center <- a[[paste0(metric, "_mean")]]
        lo <- a[[paste0(metric, "_lower_95")]]; hi <- a[[paste0(metric, "_upper_95")]]
        ok <- is.finite(center) & is.finite(lo) & is.finite(hi)
        runs <- if (any(ok)) split(which(ok), cumsum(c(TRUE, diff(which(ok)) != 1))) else list()
        for (indices in runs) if (length(indices) >= 2) polygon(c(a$year[indices], rev(a$year[indices])),
                                                                c(lo[indices], rev(hi[indices])), col = adjustcolor(strategy_colors[strategy], alpha.f = 0.13), border = NA)
      }
      for (strategy in strategies) {
        a <- d[d$strategy == strategy, ]; a <- a[order(a$year), ]
        center <- a[[paste0(metric, "_median")]]
        if (metric == "breeding_present") center <- a[[paste0(metric, "_mean")]]
        lines(a$year, center, col = strategy_colors[strategy], lwd = 2)
      }
      if (site == populations[1] && rate == exchange_rates[1]) legend("topright", legend = strategy_labels[strategies],
                                                                      col = strategy_colors[strategies], lty = 1, lwd = 2, bty = "n", cex = 0.65)
    }
    mtext(paste0(round(100 * fraction), "% of requested exchanges from a fixed external pool; shaded bands = central 95% simulation intervals"),
          outer = TRUE, side = 3, line = 1, cex = 0.95)
  }
  par(old_par)
  if (write_output_files) dev.off()
  
  #### ENDPOINT RESPONSE TO EXCHANGE INTENSITY ####
  if (write_output_files) pdf(file.path(output_directory, "exchange_intensity_comparison.pdf"), width = 12, height = 8.5)
  old_par <- par(no.readonly = TRUE)
  endpoint_labels <- c(delta_H = "Change in offspring heterozygosity", delta_breeder_H = "Change in breeder heterozygosity",
                       new_alleles_retained = "Retained alleles new to focal samples")
  displayed_rates <- c(0, exchange_rates)  # Equally spaced treatment categories keep 0 and 1 readable.
  for (metric in c("delta_H", "delta_breeder_H", "new_alleles_retained")) {
    par(mfrow = c(length(populations), length(outside_fractions)), mar = c(4.6, 4.6, 3, 1), oma = c(1, 0, 3, 0))
    plot_values <- c(endpoint_summary[[paste0(metric, "_lower_95")]], endpoint_summary[[paste0(metric, "_upper_95")]])
    plot_values <- plot_values[is.finite(plot_values)]
    limits <- if (length(plot_values) > 0) range(plot_values) else c(0, 1)
    if (diff(limits) == 0) limits <- limits + c(-0.001, 0.001)
    for (site in populations) for (fraction in outside_fractions) {
      plot(NA_real_, NA_real_, xlim = c(0.8, length(displayed_rates) + 0.2), ylim = limits, xaxt = "n", xlab = paste0("Exchange budget (swaps / ", years_per_generation, " years)"),
           ylab = endpoint_labels[metric], main = paste0(site, "; ", round(100 * fraction), "% outside"))
      axis(1, at = seq_along(displayed_rates), labels = displayed_rates); abline(h = 0, lty = 3, col = "grey70")
      for (strategy in strategies[strategies != "no_exchange"]) {
        a <- endpoint_summary[endpoint_summary$population == site & endpoint_summary$outside_fraction == fraction & endpoint_summary$strategy == strategy, ]
        a <- rbind(endpoint_summary[endpoint_summary$population == site & endpoint_summary$strategy == "no_exchange", ], a)
        a <- a[order(a$swaps_per_generation), ]
        x <- match(a$swaps_per_generation, displayed_rates)
        center <- a[[paste0(metric, "_median")]]; lo <- a[[paste0(metric, "_lower_95")]]; hi <- a[[paste0(metric, "_upper_95")]]
        ok <- is.finite(center) & is.finite(lo) & is.finite(hi)
        runs <- if (any(ok)) split(which(ok), cumsum(c(TRUE, diff(which(ok)) != 1))) else list()
        for (indices in runs) if (length(indices) >= 2) polygon(c(x[indices], rev(x[indices])),
                                                                c(lo[indices], rev(hi[indices])), col = adjustcolor(strategy_colors[strategy], alpha.f = 0.13), border = NA)
        lines(x, center, type = "b", pch = 16, lwd = 2, col = strategy_colors[strategy])
      }
      if (site == populations[1] && fraction == outside_fractions[1]) legend("topright", legend = strategy_labels[strategies[-1]],
                                                                             col = strategy_colors[strategies[-1]], lty = 1, pch = 16, bty = "n", cex = 0.7)
    }
    mtext(paste0("Outcomes at year ", n_years, "; central 95% simulation intervals"), outer = TRUE, side = 3, line = 1)
  }
  par(old_par)
  if (write_output_files) dev.off()
}

#### CLARIFIED ALLELE-RETENTION FIGURE ####
d <- summary_results[summary_results$year == max(summary_results$year), ]
rates <- sort(unique(d$swaps_per_generation))
fractions <- sort(unique(d$outside_fraction))
sites <- c("TallTimbers", "LivingstonPlace")
policies <- c("random_exchange", "relatedness_guided")
colors <- c(random_exchange = "#D55E00", relatedness_guided = "#0072B2")
labels <- c(random_exchange = "Random", relatedness_guided = "Relatedness-guided")
upper <- max(1, d$new_alleles_retained_upper_95, na.rm = TRUE)

pdf(file.path(output_directory, "allele_retention_clarified.pdf"), width = 12, height = 8.5)
old_par <- par(no.readonly = TRUE)
par(mfrow = c(length(sites), length(fractions)),
    mar = c(4.5, 5, 3, 1), oma = c(3, 0, 3, 0))

for (site in sites) for (fraction in fractions) {
  plot(NA_real_, NA_real_, xlim = c(0.8, length(rates) + 0.2),
       ylim = c(0, upper), xaxt = "n",
       xlab = "Exchange budget (swaps / 5 years)",
       ylab = "Retained alleles absent from\ninitial pooled TT + LP samples",
       main = paste0(site, "; ", 100 * fraction, "% outside"))
  axis(1, at = seq_along(rates), labels = rates)
  
  # Draw shaded simulation intervals first.
  for (policy in policies) {
    a <- d[d$population == site &
             ((d$strategy == policy & d$outside_fraction == fraction) |
                d$strategy == "no_exchange"), ]
    a <- a[order(a$swaps_per_generation), ]
    x <- match(a$swaps_per_generation, rates)
    polygon(c(x, rev(x)),
            c(a$new_alleles_retained_lower_95,
              rev(a$new_alleles_retained_upper_95)),
            col = adjustcolor(colors[policy], alpha.f = 0.13), border = NA)
  }
  
  # Add median lines over the shading.
  for (policy in policies) {
    a <- d[d$population == site &
             ((d$strategy == policy & d$outside_fraction == fraction) |
                d$strategy == "no_exchange"), ]
    a <- a[order(a$swaps_per_generation), ]
    lines(match(a$swaps_per_generation, rates),
          a$new_alleles_retained_median,
          type = "b", pch = 16, lwd = 2, col = colors[policy])
  }
  
  if (site == sites[1] && fraction == fractions[1])
    legend("topleft", legend = labels[policies], col = colors[policies],
           lty = 1, pch = 16, bty = "n", cex = 0.8)
}

mtext(paste("Year", max(d$year), "— central 95% simulation intervals"),
      side = 3, outer = TRUE, line = 1)
mtext("Zero does not imply that no alleles were introduced to an individual recipient population.",
      side = 1, outer = TRUE, line = 1, cex = 0.85)
par(old_par)
dev.off()

#### MAIN TEXT FIGURES ####
trajectory_budget <- 25
outside_to_plot <- c(0, 0.20)
endpoint_year <- 100
window_years <- 5
rates_to_plot <- c(0, 1, 15, 25, 50)
sites_to_plot <- c("TallTimbers", "LivingstonPlace")
site_labels <- c(TallTimbers = "Tall Timbers", LivingstonPlace = "Livingston Place")

# Population colors; darker shades distinguish the relatedness-guided benchmark.
population_colors <- c(TallTimbers = "dodgerblue3", LivingstonPlace = "firebrick3",
                       Arcadia = "goldenrod2", Greenwood = "seagreen4")
benchmark_colors <- c(TallTimbers = "dodgerblue4", LivingstonPlace = "firebrick4")

figure_directory <- file.path(if (exists("output_directory")) output_directory else getwd(), "main_text_panels")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

#EXPECTED OFFSPRING INBREEDING THROUGH TIME ####
curves <- data.frame(strategy = c("no_exchange", "random_exchange", "random_exchange"),
                     rate = c(0, trajectory_budget, trajectory_budget),
                     fraction = c(0, outside_to_plot), lty = c(1, 1, 2),
                     label = c("No exchange", "Random: TT-LP only",
                               paste0("Random: ", 100 * outside_to_plot[2], "% outside")))

trajectory_data <- summary_results[summary_results$population %in% sites_to_plot &
                                     (summary_results$strategy == "no_exchange" |
                                        (summary_results$strategy == "random_exchange" & summary_results$swaps_per_generation == trajectory_budget &
                                           summary_results$outside_fraction %in% outside_to_plot)), ]

# Common y-axis limits for the two population panels.
trajectory_limits <- range(c(trajectory_data$mate_kinship_lower_95,
                             trajectory_data$mate_kinship_upper_95), na.rm = TRUE)
if (diff(trajectory_limits) == 0) trajectory_limits <- trajectory_limits + c(-0.01, 0.01)

for (site_index in seq_along(sites_to_plot)) {
  site <- sites_to_plot[site_index]
  panel_colors <- c(no_exchange = "#59636E",
                    random_exchange = unname(population_colors[site]),
                    relatedness_guided = unname(benchmark_colors[site]))
  pdf(file.path(figure_directory, paste0("Projection_Fig1", LETTERS[site_index], "_", site, ".pdf")),
      width = 5, height = 5, useDingbats = FALSE)
  par(mar = c(4.7, 4.8, 3.7, 1), mgp = c(2.9, 0.8, 0), las = 1, cex = 0.9)
  plot(NA_real_, NA_real_, xlim = c(0, max(trajectory_data$year)), ylim = trajectory_limits,
       xlab = "Year", ylab = "Expected offspring inbreeding",
       main = paste0(LETTERS[site_index], "  ", site_labels[site]))
  mtext(paste0(trajectory_budget, " reciprocal swaps per ", window_years, "-year window"),
        side = 3, line = 0.4, cex = 0.8)
  abline(v = c(20, 50), col = "grey88", lty = 3)
  
  # Draw shaded bands first, then median lines.
  for (layer in c("band", "line")) for (i in seq_len(nrow(curves))) {
    a <- trajectory_data[trajectory_data$population == site &
                           trajectory_data$strategy == curves$strategy[i] &
                           trajectory_data$swaps_per_generation == curves$rate[i] &
                           trajectory_data$outside_fraction == curves$fraction[i], ]
    a <- a[order(a$year), ]
    stopifnot(nrow(a) > 0)
    center <- a$mate_kinship_median
    lo <- a$mate_kinship_lower_95; hi <- a$mate_kinship_upper_95
    color <- panel_colors[curves$strategy[i]]
    
    if (layer == "band") {
      ok <- is.finite(center) & is.finite(lo) & is.finite(hi)
      runs <- if (any(ok)) split(which(ok), cumsum(c(TRUE, diff(which(ok)) != 1))) else list()
      for (j in runs) if (length(j) >= 2) polygon(c(a$year[j], rev(a$year[j])),
                                                  c(lo[j], rev(hi[j])), col = adjustcolor(color, alpha.f = 0.12), border = NA)
    } else lines(a$year, center, col = color, lty = curves$lty[i], lwd = 2)
  }
  
  legend("topleft", legend = curves$label, col = panel_colors[curves$strategy],
         lty = curves$lty, lwd = 2, bty = "n", cex = 0.78)
  dev.off()
}

# REDUCTION IN EXPECTED INBREEDING ACROSS EXCHANGE BUDGETS ####
endpoint_data <- summary_results[summary_results$year == endpoint_year &
                                   summary_results$population %in% sites_to_plot &
                                   (summary_results$strategy == "no_exchange" |
                                      (summary_results$strategy %in% c("random_exchange", "relatedness_guided") &
                                         summary_results$swaps_per_generation %in% rates_to_plot &
                                         summary_results$outside_fraction %in% outside_to_plot)), ]

# Benefit = no exchange minus managed.
# Reverse signs AND interval endpoints when transforming the paired differences.
endpoint_data$reduction_median <- -endpoint_data$delta_kinship_median
endpoint_data$reduction_lower <- -endpoint_data$delta_kinship_upper_95
endpoint_data$reduction_upper <- -endpoint_data$delta_kinship_lower_95

# Common y-axis limits for all four response panels.
response_limits <- range(c(0, endpoint_data$reduction_lower,
                           endpoint_data$reduction_upper), na.rm = TRUE)
if (diff(response_limits) == 0) response_limits <- response_limits + c(-0.01, 0.01)
policies <- c("random_exchange", "relatedness_guided")
policy_labels <- c("Random", "Relatedness-guided (benchmark)")
panel_number <- 0

for (site in sites_to_plot) for (fraction in outside_to_plot) {
  panel_colors <- c(no_exchange = "#59636E",
                    random_exchange = unname(population_colors[site]),
                    relatedness_guided = unname(benchmark_colors[site]))
  panel_number <- panel_number + 1
  line_type <- if (fraction == 0) 1 else 2
  pdf(file.path(figure_directory, paste0("Projection_Fig2", LETTERS[panel_number], "_", site, "_",
                                         round(100 * fraction), "pct_outside.pdf")),
      width = 5, height = 5, useDingbats = FALSE)
  par(mar = c(4.7, 5.8, 3.7, 1), mgp = c(3.2, 0.8, 0), las = 1, cex = 0.9)
  
  # Budgets are displayed as equally spaced treatment categories.
  plot(NA_real_, NA_real_, xlim = c(0.8, length(rates_to_plot) + 0.2),
       ylim = response_limits, xaxt = "n",
       xlab = paste0("Reciprocal swaps per ", window_years, "-year window"),
       ylab = "Reduction in expected offspring inbreeding\n(relative to no exchange)",
       main = paste0(LETTERS[panel_number], "  ", site_labels[site]))
  axis(1, at = seq_along(rates_to_plot), labels = rates_to_plot)
  mtext(paste0("Year ", endpoint_year, "; ", round(100 * fraction), "% outside exchanges"),
        side = 3, line = 0.4, cex = 0.8)
  abline(h = 0, col = "grey65", lty = 3)
  
  for (layer in c("band", "line")) for (policy in policies) {
    a <- endpoint_data[endpoint_data$population == site &
                         ((endpoint_data$strategy == policy & endpoint_data$outside_fraction == fraction) |
                            endpoint_data$strategy == "no_exchange"), ]
    a <- a[order(a$swaps_per_generation), ]
    stopifnot(all(rates_to_plot %in% a$swaps_per_generation))
    x <- match(a$swaps_per_generation, rates_to_plot)
    center <- a$reduction_median
    lo <- a$reduction_lower; hi <- a$reduction_upper
    
    if (layer == "band") {
      ok <- is.finite(center) & is.finite(lo) & is.finite(hi)
      runs <- if (any(ok)) split(which(ok), cumsum(c(TRUE, diff(which(ok)) != 1))) else list()
      for (j in runs) if (length(j) >= 2) polygon(c(x[j], rev(x[j])),
                                                  c(lo[j], rev(hi[j])), col = adjustcolor(panel_colors[policy], alpha.f = 0.12), border = NA)
    } else lines(x, center, type = "b", pch = 16, col = panel_colors[policy],
                 lty = line_type, lwd = 2)
  }
  
  points(match(0, rates_to_plot), 0, pch = 16, col = panel_colors["no_exchange"])
  legend("topleft", legend = policy_labels, col = panel_colors[policies],
         lty = line_type, pch = 16, lwd = 2, bty = "n", cex = 0.75)
  dev.off()
}


#### SUPPLEMENTARY PROJECTION FIGURES ####
#if reloading needed run
#saved <- readRDS("~/Library/CloudStorage/GoogleDrive-janna.willoughby@gmail.com/My Drive/Willoughby lab/projects - active/red cockaded woodpecker/RCW GBS data analysis/RCW_translocations_sims/RCW_annual_projection_20261112_101721/projection_results.rds")
#summary_results <- saved$summary
#comparisons <- saved$results

#### SETTINGS AND SAVED RESULTS ####
supp_rates <- c(1, 15, 25, 50)
supp_outside <- c(0, 0.10, 0.20)
supp_years <- c(20, 50, 100) # Retained in the paired-results CSV, not separate plots
supp_endpoint <- 100
supp_trajectory_budget <- 25
supp_ps3_outside <- 0.20
supp_sites <- c("TallTimbers", "LivingstonPlace")
supp_site_labels <- c(TallTimbers = "Tall Timbers", LivingstonPlace = "Livingston Place")
supp_site_colors <- c(TallTimbers = "dodgerblue3", LivingstonPlace = "firebrick3")
supp_guided_colors <- c(TallTimbers = "dodgerblue4", LivingstonPlace = "firebrick4")
supp_policies <- c("random_exchange", "relatedness_guided")
supp_policy_labels <- c("Random", "Relatedness-guided (benchmark)")
supp_band_alpha <- 0.08
supp_directory <- file.path(if (exists("output_directory")) output_directory else getwd(), "supplementary_compact")
dir.create(supp_directory, recursive = TRUE, showWarnings = FALSE)

# In a fresh session, first read projection_results.rds into saved with readRDS().
if (exists("saved") && all(c("results", "summary") %in% names(saved))) {
  supp_summary <- saved$summary
  supp_raw <- saved$results
} else {
  stopifnot(exists("summary_results"), exists("comparisons"))
  supp_summary <- summary_results
  supp_raw <- comparisons
}
supp_raw <- supp_raw[, c("replicate", "year", "population", "strategy", "swaps_per_generation", "outside_fraction", "mate_kinship", "imported_chicks", "moved_next_breeders")]

#### PS1: TRAJECTORIES AT A SINGLE EXCHANGE BUDGET ####
ps1_budget <- supp_trajectory_budget
ps1_outside <- supp_outside
ps1_sites <- c("TallTimbers", "LivingstonPlace")
ps1_site_labels <- c(TallTimbers = "Tall Timbers", LivingstonPlace = "Livingston Place")
ps1_colors <- c(TallTimbers = "dodgerblue3", LivingstonPlace = "firebrick3")
ps1_metrics <- c("mate_kinship", "expected_offspring_H")
ps1_ylabels <- c(mate_kinship = "Expected offspring inbreeding", expected_offspring_H = "Expected offspring heterozygosity")
ps1_summary <- supp_summary

ps1_data <- ps1_summary[ps1_summary$population %in% ps1_sites & (ps1_summary$strategy == "no_exchange" | (ps1_summary$strategy == "random_exchange" & ps1_summary$swaps_per_generation == ps1_budget & ps1_summary$outside_fraction %in% ps1_outside)), ]
stopifnot(all(c("no_exchange", "random_exchange") %in% ps1_data$strategy))
ps1_curves <- data.frame(strategy = c("no_exchange", rep("random_exchange", length(ps1_outside))), outside = c(0, ps1_outside), lty = c(1, 1, 3, 2), label = c("No exchange", "Random: TT-LP only", "Random: 10% outside", "Random: 20% outside"))
# Match the main figures: 0% outside is solid and 20% outside is dashed.
# The additional 10% scenario is dotted. Outside fractions are scenario settings.

pdf(file.path(supp_directory, "PS1_trajectories.pdf"), width = 10, height = 10, useDingbats = FALSE)
par(mfrow = c(2, 2))
ps1_panel <- 0
for (metric in ps1_metrics) {
  # Common y-axis limits across the two populations for each outcome.
  limits <- unlist(ps1_data[, paste0(metric, c("_lower_95", "_upper_95"))], use.names = FALSE)
  limits <- limits[is.finite(limits)]
  stopifnot(length(limits) > 0)
  ylimits <- range(limits); padding <- max(diff(ylimits) * 0.05, 0.005)
  ylimits <- ylimits + c(-padding, padding)
  for (site in ps1_sites) {
    ps1_panel <- ps1_panel + 1
    a <- ps1_data[ps1_data$population == site, ]
    stopifnot(any(a$strategy == "no_exchange"))
    stopifnot(all(ps1_outside %in% a$outside_fraction[a$strategy == "random_exchange"]))
    colors <- c("#59636E", rep(unname(ps1_colors[site]), length(ps1_outside)))
    par(mar = c(4.7, 4.9, 4.0, 1), mgp = c(3.0, 0.8, 0), las = 1, cex = 0.85)
    plot(NA, xlim = range(a$year), ylim = ylimits, xlab = "Year", ylab = ps1_ylabels[metric])
    title(main = paste0(LETTERS[ps1_panel], ". ", ps1_site_labels[site]), line = 2, cex.main = 1.1)
    mtext(paste0(ps1_budget, " reciprocal swaps per 5-year window"), side = 3, line = 0.5, cex = 0.8)
    for (layer in c("band", "line")) for (j in seq_len(nrow(ps1_curves))) {
      b <- a[a$strategy == ps1_curves$strategy[j] & a$outside_fraction == ps1_curves$outside[j], ]
      b <- b[order(b$year), ]
      x <- b$year; middle <- b[[paste0(metric, "_median")]]
      lower <- b[[paste0(metric, "_lower_95")]]; upper <- b[[paste0(metric, "_upper_95")]]
      if (layer == "band") {
        ok <- which(is.finite(lower) & is.finite(upper))
        runs <- if (length(ok)) split(ok, cumsum(c(TRUE, diff(ok) != 1))) else list()
        for (k in runs) if (length(k) >= 2) polygon(c(x[k], rev(x[k])), c(lower[k], rev(upper[k])), col = adjustcolor(colors[j], alpha.f = 0.08), border = NA)
      } else lines(x, middle, col = colors[j], lty = ps1_curves$lty[j], lwd = 2)
    }
    legend("topleft", legend = ps1_curves$label, col = colors, lty = ps1_curves$lty, lwd = 2, bty = "n", cex = 0.73)
  }
}
dev.off()
# A/B: expected inbreeding at TT/LP. C/D: expected heterozygosity at TT/LP.

#### PS2: PAIRED BENEFITS OF GUIDANCE ####
supp_keys <- c("replicate", "year", "population", "swaps_per_generation", "outside_fraction")
a <- supp_raw[supp_raw$strategy == "random_exchange" & supp_raw$year %in% supp_years & supp_raw$population %in% supp_sites & supp_raw$swaps_per_generation %in% supp_rates & supp_raw$outside_fraction %in% supp_outside, c(supp_keys, "mate_kinship")]
b <- supp_raw[supp_raw$strategy == "relatedness_guided" & supp_raw$year %in% supp_years & supp_raw$population %in% supp_sites & supp_raw$swaps_per_generation %in% supp_rates & supp_raw$outside_fraction %in% supp_outside, c(supp_keys, "mate_kinship")]
stopifnot(!anyDuplicated(a[, supp_keys]), !anyDuplicated(b[, supp_keys]))
supp_guidance_pairs <- merge(a, b, by = supp_keys, suffixes = c("_random", "_guided"))
supp_guidance_pairs$guidance_benefit <- supp_guidance_pairs$mate_kinship_random - supp_guidance_pairs$mate_kinship_guided

supp_guidance_summary <- expand.grid(population = supp_sites, outside_fraction = supp_outside, year = supp_years, swaps_per_generation = supp_rates, stringsAsFactors = FALSE)
for (column in c("n_valid", "median", "lower_95", "upper_95", "fraction_guidance_better")) supp_guidance_summary[[column]] <- NA_real_
for (i in seq_len(nrow(supp_guidance_summary))) {
  a <- supp_guidance_summary[i, ]
  values <- supp_guidance_pairs$guidance_benefit[supp_guidance_pairs$population == a$population & supp_guidance_pairs$outside_fraction == a$outside_fraction & supp_guidance_pairs$year == a$year & supp_guidance_pairs$swaps_per_generation == a$swaps_per_generation]
  values <- values[is.finite(values)]
  supp_guidance_summary$n_valid[i] <- length(values)
  if (length(values)) {
    supp_guidance_summary[i, c("median", "lower_95", "upper_95")] <- as.list(quantile(values, c(0.5, 0.025, 0.975), names = FALSE))
    supp_guidance_summary$fraction_guidance_better[i] <- mean(values > 0)
  }
}
write.csv(supp_guidance_pairs, file.path(supp_directory, "PS2_paired_guidance_results.csv"), row.names = FALSE)
write.csv(supp_guidance_summary, file.path(supp_directory, "PS2_guidance_summary.csv"), row.names = FALSE)

supp_guidance_plot <- supp_guidance_summary[supp_guidance_summary$year == supp_endpoint, ]
all_limits <- unlist(supp_guidance_plot[, c("lower_95", "upper_95")], use.names = FALSE)
ylimits <- range(c(0, all_limits[is.finite(all_limits)]))
padding <- max(diff(ylimits) * 0.08, 0.005); ylimits <- ylimits + c(-padding, padding)
supp_source_labels <- c("TT-LP only", "10% outside", "20% outside")
supp_source_lty <- c(1, 3, 2)
pdf(file.path(supp_directory, "PS2_guidance.pdf"), width = 10, height = 5, useDingbats = FALSE)
par(mfrow = c(1, 2))
for (site_index in seq_along(supp_sites)) {
  site <- supp_sites[site_index]
  color <- unname(supp_guided_colors[site]); x <- seq_along(supp_rates)
  par(mar = c(4.7, 5.3, 4.1, 1), mgp = c(3.4, 0.8, 0), las = 1, cex = 0.85)
  plot(NA, xlim = c(0.8, length(x) + 0.2), ylim = ylimits, xaxt = "n", xlab = "Reciprocal swaps per 5-year window", ylab = "Additional reduction in expected inbreeding")
  axis(1, at = x, labels = supp_rates)
  title(main = paste0(LETTERS[site_index], ". ", supp_site_labels[site]), line = 2, cex.main = 1.1)
  mtext(paste0("Year ", supp_endpoint, "; random minus relatedness-guided"), side = 3, line = 0.5, cex = 0.8)
  abline(h = 0, col = "grey55", lty = 3)
  for (layer in c("band", "line")) for (j in seq_along(supp_outside)) {
    a <- supp_guidance_plot[supp_guidance_plot$population == site & supp_guidance_plot$outside_fraction == supp_outside[j], ]
    a <- a[match(supp_rates, a$swaps_per_generation), ]
    if (layer == "band") {
      ok <- which(is.finite(a$lower_95) & is.finite(a$upper_95))
      runs <- if (length(ok)) split(ok, cumsum(c(TRUE, diff(ok) != 1))) else list()
      for (k in runs) if (length(k) >= 2) polygon(c(x[k], rev(x[k])), c(a$lower_95[k], rev(a$upper_95[k])), col = adjustcolor(color, alpha.f = supp_band_alpha), border = NA)
    } else {
      lines(x, a$median, col = color, lwd = 2, lty = supp_source_lty[j]); points(x, a$median, col = color, pch = c(16, 17, 15)[j])
    }
  }
  legend("topleft", legend = supp_source_labels, col = color, lty = supp_source_lty, lwd = 2, pch = c(16, 17, 15), bty = "n", cex = 0.73)
}
dev.off()

#### PS3: RECRUITMENT AND EXTERNAL GENETIC CONTRIBUTIONS ####
# Recruitment is pooled over years WITHIN each replicate before calculating intervals.
# moved_next_breeders counts exchanged chicks entering the next breeding pool.
# Do not average annual recruitment fractions or count repeated adult survival as recruitment.
a <- supp_raw[supp_raw$strategy %in% supp_policies & supp_raw$year > 0 & supp_raw$year <= supp_endpoint & supp_raw$population %in% supp_sites & supp_raw$swaps_per_generation %in% supp_rates & supp_raw$outside_fraction %in% supp_outside, ]
supp_recruitment <- aggregate(a[, c("imported_chicks", "moved_next_breeders")], by = a[, c("replicate", "population", "strategy", "swaps_per_generation", "outside_fraction")], FUN = sum)
supp_recruitment$recruitment_fraction <- with(supp_recruitment, ifelse(imported_chicks > 0, moved_next_breeders / imported_chicks, NA_real_))
write.csv(supp_recruitment, file.path(supp_directory, "PS3_recruitment_by_replicate.csv"), row.names = FALSE)

supp_recruitment_summary <- expand.grid(population = supp_sites, strategy = supp_policies, outside_fraction = supp_outside, swaps_per_generation = supp_rates, stringsAsFactors = FALSE)
for (column in c("n_valid", "median", "lower_95", "upper_95", "pooled_fraction")) supp_recruitment_summary[[column]] <- NA_real_
for (i in seq_len(nrow(supp_recruitment_summary))) {
  a <- supp_recruitment_summary[i, ]
  b <- supp_recruitment[supp_recruitment$population == a$population & supp_recruitment$strategy == a$strategy & supp_recruitment$outside_fraction == a$outside_fraction & supp_recruitment$swaps_per_generation == a$swaps_per_generation, ]
  values <- b$recruitment_fraction[is.finite(b$recruitment_fraction)]
  supp_recruitment_summary$n_valid[i] <- length(values)
  if (length(values)) supp_recruitment_summary[i, c("median", "lower_95", "upper_95")] <- as.list(quantile(values, c(0.5, 0.025, 0.975), names = FALSE))
  if (sum(b$imported_chicks) > 0) supp_recruitment_summary$pooled_fraction[i] <- sum(b$moved_next_breeders) / sum(b$imported_chicks)
}
write.csv(supp_recruitment_summary, file.path(supp_directory, "PS3_recruitment_summary.csv"), row.names = FALSE)

# Combine the recruitment summaries with endpoint allele-retention and ancestry summaries.
# "New alleles" means absent from the COMBINED initial TT + LP sample, not newly mutated.
# External ancestry means ancestry from arrivals during the simulation, not historical ancestry.
supp_ps3 <- supp_recruitment_summary[, c("population", "strategy", "outside_fraction", "swaps_per_generation", "n_valid", "median", "lower_95", "upper_95")]
supp_ps3$metric <- "recruitment_fraction"
for (metric in c("new_alleles_retained", "mean_external_ancestry")) {
  a <- supp_summary[supp_summary$year == supp_endpoint & supp_summary$population %in% supp_sites & supp_summary$strategy %in% supp_policies & supp_summary$swaps_per_generation %in% supp_rates & supp_summary$outside_fraction %in% supp_outside, c("population", "strategy", "outside_fraction", "swaps_per_generation", paste0(metric, c("_n_valid", "_median", "_lower_95", "_upper_95")))]
  names(a)[5:8] <- c("n_valid", "median", "lower_95", "upper_95")
  a$metric <- metric
  supp_ps3 <- rbind(supp_ps3, a)
}
write.csv(supp_ps3, file.path(supp_directory, "PS3_panel_summaries.csv"), row.names = FALSE)

supp_ps3_labels <- c(recruitment_fraction = "Exchanged chicks recruiting (%)", new_alleles_retained = "Retained alleles absent from initial TT + LP", mean_external_ancestry = "Ancestry from simulated external arrivals (%)")
supp_ps3_plot <- supp_ps3[supp_ps3$outside_fraction == supp_ps3_outside, ]
pdf(file.path(supp_directory, "PS3_recruitment_and_external_contributions.pdf"), width = 10, height = 15, useDingbats = FALSE)
par(mfrow = c(3, 2))
supp_panel <- 0
for (metric in unique(supp_ps3_plot$metric)) {
  scale_factor <- if (metric == "new_alleles_retained") 1 else 100
  all_limits <- supp_ps3_plot$upper_95[supp_ps3_plot$metric == metric] * scale_factor
  all_limits <- all_limits[is.finite(all_limits)]
  ymax <- if (length(all_limits)) max(all_limits) else 0
  ylimits <- c(0, max(ymax * 1.12, if (scale_factor == 100) 1 else 5))
  for (site in supp_sites) {
    supp_panel <- supp_panel + 1
    outside <- supp_ps3_outside
    a <- supp_ps3_plot[supp_ps3_plot$population == site & supp_ps3_plot$metric == metric, ]
    colors <- c(unname(supp_site_colors[site]), unname(supp_guided_colors[site]))
    x <- seq_along(supp_rates)
    par(mar = c(4.7, 5.3, 4.1, 1), mgp = c(3.4, 0.8, 0), las = 1, cex = 0.85)
    plot(NA, xlim = c(0.8, length(x) + 0.2), ylim = ylimits, xaxt = "n", xlab = "Reciprocal swaps per 5-year window", ylab = supp_ps3_labels[metric])
    axis(1, at = x, labels = supp_rates)
    title(main = paste0(LETTERS[supp_panel], ". ", supp_site_labels[site]), line = 2, cex.main = 1.05)
    timing <- if (metric == "recruitment_fraction") sprintf("Years 1-%d", supp_endpoint) else sprintf("Year %d", supp_endpoint)
    mtext(sprintf("%s; %g%% outside exchanges", timing, 100 * outside), side = 3, line = 0.5, cex = 0.85)
    for (layer in c("band", "line")) for (j in seq_along(supp_policies)) {
      b <- a[a$strategy == supp_policies[j], ]; b <- b[match(supp_rates, b$swaps_per_generation), ]
      middle <- b$median * scale_factor; lower <- b$lower_95 * scale_factor; upper <- b$upper_95 * scale_factor
      if (layer == "band") {
        ok <- which(is.finite(lower) & is.finite(upper))
        runs <- if (length(ok)) split(ok, cumsum(c(TRUE, diff(ok) != 1))) else list()
        for (k in runs) if (length(k) >= 2) polygon(c(x[k], rev(x[k])), c(lower[k], rev(upper[k])), col = adjustcolor(colors[j], alpha.f = supp_band_alpha), border = NA)
      } else {
        lines(x, middle, col = colors[j], lwd = 2, lty = if (outside == 0) 1 else 2); points(x, middle, col = colors[j], pch = c(16, 17)[j])
      }
    }
    legend("topleft", legend = supp_policy_labels, col = colors, lty = if (outside == 0) 1 else 2, lwd = 2, pch = c(16, 17), bty = "n", cex = 0.73)
  }
}

dev.off()


