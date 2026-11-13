# Red-cockaded woodpecker reintroduction and nestling-exchange simulations

R code for two complementary analyses of red-cockaded woodpecker (*Dryobates borealis*) conservation in the Red Hills region of northern Florida and southwestern Georgia:

1. **Founding history and expected differentiation:** how founder sampling, shared ancestry, and genetic drift affect expected differentiation between Tall Timbers and Livingston Place.
2. **Nestling-exchange projections:** how exchange frequency, donor access, and chick-selection strategy affect relatedness among future breeding pairs and expected inbreeding in their offspring.

## Scripts and requirements

| File | Purpose | Required input |
| --- | --- | --- |
| `RCW_translocation_simulation.R` | Simulate founding, drift, and multilocus FST | Parameter values; no empirical genotype file |
| `RCW_projections.R` | Project annual reproduction, survival, recruitment, and nestling exchanges; generate main and supplementary figures | `RCW_genotype_panel_10000.csv` |

## 1. Founding history and expected differentiation

### Model

Each replicate generates a common ancestral allele-frequency vector at independent, neutral, diploid biallelic loci. Tall Timbers is founded by sampling alleles from this pool and experiences drift before Livingston Place is founded. A specified fraction of LP founders is sampled from contemporary simulated TT frequencies; the remainder is sampled from the original ancestral pool. Both populations then experience additional drift.

Final samples are drawn from the simulated populations, and multilocus FST is calculated using both Hudson and Weir–Cockerham estimators. Both are ratios of summed components across loci. Negative estimates are retained. The selected `fst_estimator` determines the scenario summaries and plots.

### Supplied settings

| Parameter | Supplied value |
| --- | --- |
| Random seed | 2112 |
| Replicates per scenario | 1,000 |
| Generated loci per replicate, before filtering | 1,000 |
| TT / LP sampled individuals | 10 / 10; example counts to replace for an empirical comparison |
| Founder model | Baseline: 12 TT and 18 LP founders |
| Fraction of LP founders originating from TT | 0, 0.5, 1 |
| Effective population size, assigned equally to TT and LP | 10, 25, 50, 100 |
| Generations before / after LP establishment | 2 / 1 |
| Ancestral frequency distribution | Uniform between 0.05 and 0.95 |
| Summary estimator | Hudson |
| Observed FST comparison | Disabled (`observed_fst <- NA_real_`) |

### Outputs

Results are saved beneath `RCW_simulation_<timestamp>/`:

| Output | Contents |
| --- | --- |
| `scenarios.csv` | Parameter combinations and scenario identifiers |
| `replicate_results.csv` | Founder allocations, retained locus counts, and both FST estimates |
| `summary.csv` | Scenario medians and central 95% simulation intervals |
| `simulation_results.rds` | Parameters, scenarios, summaries, and replicate results |
| `parameters.txt`, `R_session_info.txt` | Settings and software environment |
| `scenario_comparison.pdf` | FST versus effective population size, by TT source fraction |
| `TT_source_percent_reduction.pdf` | Percentage reduction in median FST relative to source fraction zero |
| `scenario_comparison_Ne.pdf` | FST versus TT source fraction, by effective population size |

## 2. Nestling-exchange projections

### Input format

`RCW_genotype_panel_10000.csv` must contain one row per sampled individual and the following columns:

| Column | Meaning |
| --- | --- |
| `sample_id` | Sample identifier, read as character |
| `population` | `TallTimbers`, `LivingstonPlace`, `Arcadia`, or `Greenwood` for study individuals |
| `age_class_recorded` | The script selects hatch-year birds coded `HY` |
| `sex` | `Male` or `Female`; other or missing values are assigned a sex within each replicate |
| `family_group_id` | Recorded population–year–cluster/rearing-group identifier |
| `SNP_...` | Numeric genotype dosage: 0, 1, or 2 alternative-allele copies |

The 10,000 SNPs are divided once into 5,000 selection loci and 5,000 outcome-evaluation loci. Genomic coancestry uses allele frequencies pooled across the starting empirical individuals and held fixed throughout the simulation. The marker split is held constant across scenarios and replicates; it is not independent-population validation.

### Model and supplied settings

| Parameter | Supplied value |
| --- | --- |
| Random seed | 2112 |
| Replicates per scenario | 100 |
| Projection length | 100 years: 20 five-year budget windows |
| Target breeding pairs per year | TT: 20; LP: 10 |
| Chicks per breeding pair | 4 |
| Annual adult survival | Male: 0.80; female: 0.75 |
| Offspring recruitment eligibility probability | 1, conditional on available vacancies |
| Requested reciprocal swaps per five-year window | 1, 15, 25, 50 |
| Outside-source probabilities | 0, 0.10, 0.20 |
| Strategies | No exchange, random exchange, relatedness-guided exchange |

The design contains **25 scenarios**: one no-exchange control plus two exchange strategies × four budgets × three outside-source probabilities. Pair capacities are scenario inputs and are not inferred from the sampled nest counts.

The model runs **annually with adult survival and overlapping age classes**. Despite their names, `n_generations` and `years_per_generation` organize five-year exchange-budget windows; they do not impose complete breeder replacement every five years. The output column `generation` identifies the budget window.

Initial breeding slots are filled by sex-specific resampling of observed chick genotypes. These are breeding-population proxies, not reconstructed parental genotypes. Birds pair randomly each year, excluding shared rearing groups. Offspring inherit independently sampled parental alleles at each locus. Surviving adults remain, and eligible offspring fill sex-specific vacancies. Unrecruited juveniles are not retained as helpers or floaters, and pair bonds are not preserved.

Internal exchanges move one same-sex chick in each direction between TT and LP. Each natal nest can contribute at most one moved chick per year. External exchanges replace an internal event within the requested budget; they do not increase that budget. Incoming external genotypes are sampled from the unchanged Arcadia–Greenwood reference pool. Exports and descendants do not alter that pool. Equal external-recipient weights imply equal probabilities of selecting TT or LP, not exactly equal realized counts.

Random selection samples feasible chick pairs. Relatedness-guided selection chooses the feasible exchange minimizing mean genomic coancestry among permitted opposite-sex pairings in the current chick cohort, averaging the affected population scores for internal exchanges. It does not optimize against every surviving adult or known future mate. Exchanges are attempted even when no candidate is predicted to improve on no exchange. Unfilled requests are recorded and are not carried forward.

### Outputs and interpretation

Results are saved beneath `RCW_annual_projection_<timestamp>/`:

| Output | Contents |
| --- | --- |
| `projection_replicates.csv` | Annual population outcomes and paired differences from no exchange |
| `projection_summary.csv` | Means, medians, valid counts, and central 95% simulation intervals |
| `endpoint_summary.csv` | Selected final-year summary columns |
| `projection_results.rds` | Complete saved result objects, parameters, and outcome summaries |
| `annual_exchange_counts.csv`, `five_year_exchange_counts.csv` | Requested and completed exchange counts |
| `exchange_budget_summary.csv`, `exchange_scenarios.csv` | Budget summaries and scenario definitions |
| `simulated_swaps.csv` | Completed exchange events |
| `starting_genotype_proxies.csv` | Empirical individuals used to fill initial breeding slots |
| `R_session_info.txt` | Software environment |

Key metric definitions:

- **`mate_kinship`:** genomic coancestry between actual breeding partners, interpreted as estimated expected offspring inbreeding relative to the fixed reference.
- **`expected_offspring_H`:** expected offspring heterozygosity from parental allele-transmission probabilities, evaluated at the reserved SNPs.
- **`new_alleles_retained`:** alleles present in the subsequent breeding pool but absent from the **combined initial TT–LP empirical sample**. Zero does not imply that no allele was introduced to an individual recipient population.
- **`mean_external_ancestry`:** mean ancestry among subsequent breeders attributable to external arrivals during the projection, not historical ancestry inferred from genotypes.
- **`moved_next_breeders`:** exchanged chicks recruited into the subsequent breeding pool; this does not repeatedly count their later adult survival.

### Figures

The initial plotting blocks generate `projection_trajectories.pdf`, `exchange_intensity_comparison.pdf`, and `allele_retention_clarified.pdf`. Appended publication-figure blocks produce:

| Location | Figure |
| --- | --- |
| `main_text_panels/Projection_Fig1*.pdf` | Two separate 5 × 5-inch population panels: inbreeding trajectories at 25 swaps per five years |
| `main_text_panels/Projection_Fig2*.pdf` | Four separate 5 × 5-inch panels: final-year reductions across exchange budgets |
| `supplementary_compact/PS1_trajectories.pdf` | One page, four panels: inbreeding and heterozygosity trajectories under random exchange at 25 swaps, with 0%, 10%, or 20% outside sources |
| `supplementary_compact/PS2_guidance.pdf` | One page, two panels: paired guidance benefits at year 100 across budgets and source scenarios |
| `supplementary_compact/PS3_recruitment_and_external_contributions.pdf` | One page, six panels: recruitment, retained alleles, and external ancestry across budgets at 20% outside sources |

## Estimator references

- Bhatia, G., Patterson, N., Sankararaman, S., and Price, A. L. 2013. Estimating and interpreting FST: The impact of rare variants. *Genome Research* 23:1514–1521. https://doi.org/10.1101/gr.154831.113
- Weir, B. S., and Cockerham, C. C. 1984. Estimating F-statistics for the analysis of population structure. *Evolution* 38:1358–1370. https://doi.org/10.1111/j.1558-5646.1984.tb05657.x
- VanRaden, P. M. 2008. Efficient methods to compute genomic predictions. *Journal of Dairy Science* 91:4414–4423. https://doi.org/10.3168/jds.2007-0980
