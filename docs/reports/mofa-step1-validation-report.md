# MOFA Step 1 Validation & Factor Audit Report

**Date:** November 19, 2026 (PDBP rerun with library-size normalization; supersedes the October 7, October 6, October 5 and September 25, 2026 versions)  
**Discovery cohort:** PPMI baseline ($N = 1,329$ participants with factor scores)  
**Validation cohort:** PDBP baseline ($N = 1,248$ projected participants; 107 with proteomics)  
**Main model (Oct 7):** Multi-group MOFA (`mofapy2` 0.7.5, `groups = sex`, 10 factors, `scale_views=True`, `center_groups=True` by default) on **RNA-seq + CSF proteomics**. Trained until converged, with proteomics zeros treated as missing. File: `ppmi_mofa_outputs_noplasma/model_seed_0.hdf5`.  
**Views:** whole blood RNA-seq (5,000 genes) and CSF proteomics (249 proteins). Plasma proteomics was dropped (Section 5).  
**Sensitivity model:** the same model with plasma as a third view (159 proteins after zero QC), `ppmi_mofa_outputs/model_seed_0.hdf5`. Its factors, variance explained, associations and stability match the main model to within 0.1 (Section 2).

---

## 0. Changes from earlier versions

- **Batch correction was rerun** (commit `409a2c9`, Oct 5). `plate` has no association with any factor in either cohort, so the Sep 25 finding that Factors 4 and 6 were dominated by plate no longer applies.
- **The "58% of total model variance" figure for Factor 1 (Sep 25) was wrong.** That figure is the variance Factor 1 explains *within the RNA-seq view*. Per-view R² values can't be added up to a model total.
- **The variance-explained replication code was fixed** (Section 4). The Oct 5 morning figure (Factor 1 ≈ 120%) is superseded.
- **The model was retrained on Oct 6** with three changes:
  - proteomics zeros treated as missing, plus zero-based QC (Section 5);
  - training run until converged (the earlier model stopped at 20 iterations; all 5 seeds now converge at 220);
  - Tau and AlphaZ saved, so the PDBP projection uses them.

  **Factor numbers changed.** The aging axis was Factor 4 and is now **Factor 10**. The PD axis was Factor 8 and is now **Factor 7**. The CSF factor is still **Factor 2**. Factor 10's negative R² in CSF, a sign of incomplete training, is gone.
- **Seeds aren't independent.** mofapy2 starts every run from a PCA of the data, so the 5 seeds give the same solution: final ELBOs within 0.07 of −5,017,117.7, and every factor matched with |r| = 1.0. Stability is now assessed by subsampling (Section 2).
- **PDBP RNA-seq is now normalized for library size (Nov 19).** Before this, the PDBP notebook cell applied log1p to raw counts, while PPMI used DESeq2 size factors. `preprocess_rnaseq.py` now takes `n_top_genes=None`, which applies size factors estimated within PDBP to all genes, and PDBP was reprocessed, batch-corrected, projected and re-evaluated for the main and single-group models. PPMI values are unchanged. PDBP changes:
  - **Factor 4's** PDBP RNA-seq R² fell from 2.9–3.0% to 0.8–0.9%, and its PD association now meets the replication rule. It is probably technical, though (Section 3, technical check).
  - **Factor 1's** PDBP RNA-seq R² rose from 41–44% to 44–49%.
  - The RNA-seq per-factor R² rank agreement moved from 0.95 / 0.70 (women / men) to 0.85 / 0.78.
  - **No conclusion on the sex question changed.**
  
  Pre-fix PDBP results are kept in `results/mofa/pre_normalization_fix/`. The three-view sensitivity model's PDBP results were not rerun and still predate the fix.

## 1. Variance explained per view and sex (main model: RNA-seq + CSF)

PPMI values come from the model. PDBP values come from the projection, using the same formula. The model file contains the saved Tau and AlphaZ, which `project_pdbp.py` uses when present; the no-plasma projection log hasn't been checked to confirm this. PDBP values are from the Nov 19 rerun with library-size normalization. Figure: `pdbp_mofa_outputs_noplasma/cohort_variance_replication_comparison.png` on Terra. Tables: `results/mofa/main_noplasma/cohort_variance_replication_comparison.csv`, `results/mofa/main_noplasma/variance_explained_total.csv`, `results/mofa/main_noplasma/pdbp_variance_explained_total.csv`.

**Total R² per view (%)**

| View | PPMI Female | PPMI Male | PDBP Female | PDBP Male | PDBP n (F / M) |
|---|---|---|---|---|---|
| RNA-seq | 66.4 | 69.9 | 52.4 | 57.7 | 554 / 694 |
| CSF proteomics | 24.3 | 27.5 | 27.1 | 26.3 | 41 / 66 |

**Per-factor R² for factors above 0.5% in some view (%, Female / Male)**

| Factor | View | PPMI | PDBP | PDBP before the Nov 19 fix |
|---|---|---|---|---|
| Factor 1 | RNA-seq | 56.9 / 61.7 | 43.5 / 48.8 | 40.7 / 44.3 |
| Factor 2 | CSF | 23.5 / 26.5 | 26.9 / 26.3 | 26.9 / 26.2 |
| Factor 3 | RNA-seq | 3.7 / 2.8 | 3.1 / 2.5 | 3.1 / 2.6 |
| Factor 4 | RNA-seq | 1.4 / 1.5 | 0.8 / 0.9 | 2.9 / 3.0 |
| Factor 5 | RNA-seq | 1.0 / 0.9 | 1.0 / 1.0 | 0.9 / 1.0 |
| Factor 6 | CSF / RNA-seq | 0.4 / 0.6 (CSF); 0.5 / 0.5 (RNA) | 0.2 / 0.0 (CSF); 0.6 / 0.5 (RNA) | 0.3 / −0.1 (CSF); 0.6 / 0.4 (RNA) |
| Factor 7 | RNA-seq | 1.0 / 0.7 | 1.0 / 1.1 | 1.0 / 1.1 |
| Factor 8 | RNA-seq | 0.9 / 0.9 | 0.9 / 0.9 | 1.1 / 0.6 |
| Factor 9 | RNA-seq | 0.6 / 0.5 | 0.3 / 0.3 | 0.7 / 0.6 |
| Factor 10 | RNA-seq | 0.5 / 0.4 | 0.7 / 0.7 | 0.7 / 0.7 |

Spearman agreement of per-factor R² (women / men): RNA-seq 0.85 / 0.78; CSF 0.30 / 0.52. The CSF values are low because only Factor 2 explains more than 1% of CSF, so the other factors' near-zero values rank at random.

**What this shows:**
- **CSF replicates closely.** Factor 2 explains about 24–27% of CSF in both cohorts. Treating CSF's few zeros as missing roughly doubled it from the earlier model's 12%.
- **RNA-seq replicates.** Factor 1 explains less in PDBP (44–49% vs 57–62%), as expected for data the model wasn't trained on. Factors 3, 5–8 and 10 replicate at similar levels; Factors 4 and 9 explain about half as much in PDBP. Before library-size normalization, Factor 4 explained twice as much in PDBP (2.9–3.0%), because it had partly absorbed sequencing depth.
- **No factor is specific to one sex.** Every factor explains similar shares in women and men, in both views and both cohorts. The sex question rests on the sex × diagnosis interaction test (Section 3).
- **Removing plasma changed nothing else.** Every per-factor R² in RNA-seq and CSF is within 0.01 percentage points of the three-view model, and the factor numbering is the same.

## 2. Factor–metadata associations and stability (main model)

PPMI values use deduplicated metadata (one row per participant) and the off-medication baseline MDS-UPDRS III exam (rerun Oct 7). Compared with the earlier duplicated merge, they moved by at most 1.7 in −log10 p, and no conclusion changed.

Associations are −log10 p (Pearson for numeric covariates, ANOVA for categorical), with both sexes pooled. Tables: `results/mofa/main_noplasma/factor_metadata_associations_ppmi.csv` and `results/mofa/main_noplasma/factor_metadata_associations_pdbp.csv`. PDBP values are from the Nov 19 rerun.

Stability comes from `subsample_stability.py`: 10 subsets of 80% of PPMI participants, drawn within each sex and each retrained. It reports the |r| of factor scores between the full model and each subset. A factor is stable if it's reproduced with |r| ≥ 0.8 in at least 90% of subsets. Variance explained isn't part of the rule: factors explaining 0.5–1% of 5,000 genes carry the strongest replicated clinical signals.

| Factor | Main signal (PPMI / PDBP) | Max R² in PPMI (%) | Subset \|r\| median / min | Reproduced in | Weights \|r\| / top-5% overlap (median) | Reading |
|---|---|---|---|---|---|---|
| **Factor 7** | Diagnosis 4.5 / 12.5; MDS-UPDRS III 9.1 / 9.6; age 9.5 / 13.1; MoCA 5.4 / 2.2 | 1.0 | 0.97 / 0.88 | 10/10 | 0.88 / 0.37 | **PD axis.** Replicates and is stable. Also tracks age |
| **Factor 10** | Age 25.4 / 28.4; MoCA 3.5 / 10.0; MDS-UPDRS III 6.6 / 1.1 | 0.5 | 0.99 / 0.72 | 9/10 | 0.78 / 0.50 | **Aging axis.** Replicates and is stable |
| Factor 2 | Age 7.1 / 9.2; MDS-UPDRS III 3.5 / 2.7; MoCA 2.1 / 2.6 | 26.5 (CSF) | 1.00 / 0.21 | 9/10 | 0.99 / 0.85 | The CSF axis, mildly age-related |
| Factor 3 | Age 4.9 / 0.1; diagnosis 2.3 / 0.8 | 3.7 | 1.00 / 1.00 | 10/10 | 0.92 / 0.80 | Stable, but its clinical signal doesn't replicate |
| **Factor 4** | Diagnosis ≤ 2.6 / 2.9 (was 0.2 before the Nov 19 fix) | 1.5 | ≥ 0.99 / ≥ 0.99 | 10/10 | ≥ 0.94 / ≥ 0.76 | **Probably technical.** Its PD association meets the replication rule but doesn't survive adjustment for genes detected in PPMI (Section 3) |
| Factor 5 | Age 2.6 / 3.1; MDS-UPDRS III 0.6 / 4.7 | 1.0 | 0.97 / 0.88 | 10/10 | 0.91 / 0.48 | Weak |
| Factor 6 | Diagnosis 0.1 / 7.4 | 0.6 | 0.95 / 0.75 | 8/10 | 0.92 / 0.45 | PDBP only, and **not stable** |
| Factor 9 | MoCA 0.0 / 3.5; diagnosis 1.9 / 0.2 | 0.6 | 1.00 / 0.98 | 10/10 | 0.98 / 0.71 | Doesn't replicate |
| Factors 1, 8 | No clinical association in PPMI (all ≤ 2.6); PDBP age 0.0 / 3.5 | 61.7 / 0.9 | ≥ 0.99 / ≥ 0.99 | 10/10 | ≥ 0.94 / ≥ 0.76 | Stable. Factor 1 is likely blood cell composition. Factor 8's adjusted PD effect replicates in men only (Section 3) |
| any | plate 0.0 / 0.0; sex ≤ 0.1 / ≤ 0.2 | | | | | No batch signal. Sex ≈ 0 by design after centering within each sex (Section 3) |

**Weakest subsets:**
- **Factor 2 (CSF)** reached only |r| = 0.21 in subset 8, the same subset in both the main and three-view runs, with |r| ≥ 0.99 in all other subsets. Checked in the three-view run (Oct 7): subset 8 still has a CSF factor, but it explains half as much CSF variance (10.3 / 11.4% vs 21.5 / 26.8% in subset 0), and its scores don't match the main CSF axis. That fit settled on a secondary CSF pattern. The same participants give the same result without plasma, so this reflects that particular subset of about 283 CSF participants rather than chance in training. One such subset in 10 still meets the stability rule. A likely cause: the starting factors come from a PCA dominated by the 5,000 RNA-seq genes.
- **Factor 6** blends with Factor 5 (subset 6: 0.75; subset 3: 0.78). Both explain about 0.6–1%, which is why Factor 6 isn't stable.
- **Factor 10 (aging)** blended with Factor 9 in subset 6 only (0.72), and reached |r| ≥ 0.97 elsewhere.
- **Factor 7 (PD)** reached |r| ≥ 0.88 in every subset.
- **If changes are wanted later, test in this order:** (1) 15 factors instead of 10, as `config.yaml` already specifies (`mofa_num_factors`); (2) `weight_views=True` in `set_train_options`, which weights the objective by the number of features per view so RNA-seq doesn't dominate CSF. Option 2 changes what the model optimizes and would need to be documented as a design choice.

**Weight agreement across subsets:** the median |r| of the weights is 0.78–0.99, and the median top-5% overlap is 0.37–0.87. Both are higher than in the three-view run, where plasma's noisy weights pulled the medians down. Factor 7 (PD) has the lowest top-5% overlap (0.37), so pathway enrichment should use rank-based methods (for example GSEA on the full weight vector) rather than a hard top-5% cut-off.

**Candidates for the sex × diagnosis interaction test: Factor 7 (PD) and Factor 10 (aging).** Factor 7 also tracks age (14.0 in PDBP), so the interaction model must adjust for age.

**Three-view sensitivity model:** per-factor R² in RNA-seq and CSF, associations (largest difference 0.1 in −log10 p) and subsampling verdicts all match the main model. The plasma view itself explained 1.2–1.5% in PPMI and −0.5 / −0.2% in PDBP. Removing plasma doesn't change any conclusion.

Caveats:
- **PPMI rows are duplicated in the merge.** 1,329 participants produced 1,432 merged rows, which inflates PPMI significance somewhat. Deduplicate before reporting final p-values.
- **Clinical associations are not adjusted for age.** Factor 7's and Factor 10's associations with MDS-UPDRS and MoCA should be refit with age as a covariate.
- **The plate ANOVA tests only differences in means.** It reads close to zero after batch correction by design.
- **The diagnosis variable may have more than two levels:** it's tested by ANOVA. Confirm its categories (for example PD / control / prodromal) before setting up the interaction test.

## 3. Interpreting sex in a multi-group MOFA model

**What the model does.** With `groups = sex`, mofapy2 centers every feature within each sex before training (`center_groups=True` by default). The weights W are shared by both sexes. Each sex gets its own factor scores Z and its own factor-activity prior, so a factor can be active (explain variance) in one sex and not the other.

**Consequences:**
1. **PPMI factor means can't differ by sex.** Within-group centering removes mean differences between men and women, so each factor has a mean of about 0 in both groups. The sex ANOVA in PPMI (−log10 p ≤ 0.3) is therefore uninformative by construction. It does not show that the sexes are similar.
2. **The PDBP sex associations come from preprocessing, not biology.** The projection (`project_pdbp.py`) pools both sexes and does no centering within each sex. Any mean differences between men and women in the PDBP data pass straight through the shared weights. This is especially true of sex-chromosome genes (XIST, RPS4Y1 and others), which are not filtered out of the RNA-seq view. That's why sex shows up strongly in PDBP (Factor 4: 27.9) but not in PPMI. Comparing sex associations across the two cohorts is therefore not valid. *Resolved Oct 5:* after re-projecting with centering within each sex, the PDBP sex associations drop to ≤ 0.1, matching PPMI.
3. **"Sex-specific loading patterns" can't be tested.** The weights are the same for both sexes, so the spec's item "(3) conservation of sex-specific factor loading patterns" can't be measured as written.

**Does this invalidate the premise?** No, but it narrows the question MOFA can answer. Multi-group MOFA is designed to find which *axes of variation* are shared by both sexes and which are active in only one. It is not designed to find *differences in mean level* between sexes (the MOFA2 documentation says this explicitly for multi-group models). The ADR decision ("separate shared factors from sex-specific factors") still holds if "sex-specific" means a factor whose variance is present in one sex and absent in the other. It does not hold if it means "men and women score differently on this factor."

**Why the grouping stays as sex, not diagnosis.** With `groups = diagnosis`, every feature would be centered separately in PD and in controls. That would remove the average PD-vs-control difference, which is the main disease signal, before training. The model would then describe variation within PD versus within controls, which isn't a sex question. A four-group model (sex × diagnosis) would remove both the sex and the disease differences in mean level, and the proteomics groups would be small. Centering within each sex keeps the PD-vs-control differences inside each sex, so `groups = sex` is the right primary design.

**The right readouts for the sex question:**
- **Factor activity per sex.** Compare R² per factor between the female and male groups (already in `variance_explained/r2_per_factor`). Factors with substantial R² in one group and near-zero in the other are the sex-specific ones.
- **The sex × diagnosis interaction test (primary result).** For each factor, fit
  `factor ~ diagnosis × sex + age` (plus MDS-UPDRS III × sex and MoCA × sex for clinical severity), using one row per participant. The interaction term tests whether the factor tracks PD *differently* in men and women. Report the diagnosis effect within each sex next to the interaction p-value, and correct for the number of factors (Benjamini–Hochberg). Because centering within each sex keeps the PD-vs-control differences inside each sex, this test is valid on the current model's factors. Run it in PPMI, then in the re-projected PDBP, and count a factor as replicated when the interaction goes in the same direction in both.
- **Single-group sensitivity check.** Train a second MOFA on PPMI with a single group (no sex grouping), the same views and the same preprocessing, then run the same `diagnosis × sex + age` interaction test on its factors. Average sex differences are not removed in this model, so a factor whose PD effect differs by sex can show up directly in the interaction test rather than only through per-sex factor activity. Match its factors to the sex-grouped model's factors by correlating their weights. Sex-differential PD axes found in both models are the strongest result. Here, decide explicitly whether to drop sex-chromosome features (XIST, RPS4Y1 and others), because in a single-group model they will form a pure sex factor.
- **Projection that matches training.** Before projecting PDBP, center and scale each view within each sex the same way PPMI was. This is now implemented in `project_pdbp.py` (Section 4). Then compare per-sex R² and per-sex clinical associations across cohorts.
- **Imbalance between PD and sex.** PD is more common in men, and the PPMI PD group is mostly male. Centering within each sex therefore subtracts averages taken over different PD/control mixes. Including diagnosis and age as covariates in the interaction model handles this. Raw factor-score comparisons don't.

### Sex × diagnosis interaction results (main model; PPMI Oct 7, PDBP rerun Nov 19)

From `mofa/scripts/clinical_models.py`, with diagnosis PD (cases) vs Ctrl (controls). PPMI's 21 "Other" participants are excluded; PDBP has none. That leaves PPMI with 753 PD / 555 controls (F 316 / 290, M 437 / 265) and PDBP with 767 PD / 481 controls (F 281 / 273, M 486 / 208). Factors are standardized within each cohort and models are adjusted for age. A result counts as replicated when PPMI q < 0.05 (BH across factors), PDBP p < 0.05, and the signs agree.

- **No factor's PD effect differs by sex.** None of the 30 interaction tests (10 factors × diagnosis, MDS-UPDRS III, MoCA) replicates, and none reaches PPMI q < 0.05 (smallest q = 0.21). Two PDBP-only nominal results don't hold in PPMI: Factor 3 × UPDRS III (p = 0.039) and Factor 8 × MoCA (p = 0.045). Before the Nov 19 fix, Factor 1 × diagnosis was also nominal in PDBP (p = 0.037; now 0.074).
- **How small an effect the test can detect:** the standard error of the diagnosis × sex term is about 0.11–0.12 SD in each cohort. The test therefore had about 80% power to detect sex differences in the PD effect of about 0.3 SD or more (before multiple-testing correction). Smaller sex differences can't be ruled out.
- **Factor 7 is the main replicated PD axis.** PD is associated with lower Factor 7 scores, adjusted for age and sex: −0.19 log-odds per SD in PPMI (q = 0.001) and −0.44 in PDBP (p < 0.001; −0.46 before the Nov 19 fix). The PD effect is similar in each sex, age-adjusted:

  | | PPMI | PDBP |
  |---|---|---|
  | Women | −0.16 SD (p = 0.07) | −0.44 SD (p < 0.001) |
  | Men | −0.21 SD (p = 0.004) | −0.36 SD (p < 0.001) |

  Only the men's estimate meets the replication rule, because PPMI has more men. The interaction term is null (PPMI p = 0.62, PDBP p = 0.66). Factor 7 is therefore a replicated, sex-invariant multi-omic PD signature.
  - **Caveat:** the PD effect is about twice as large in PDBP as in PPMI. PPMI enrolls recently diagnosed PD; PDBP includes longer disease duration and people on dopaminergic treatment, either of which could change the blood transcriptome. Factor 7 may partly track disease duration or medication. Check before interpreting it biologically, for example by adding disease duration and medication status as covariates in PDBP.
- **Factor 4 meets the replication rule but is probably technical (new since the Nov 19 fix).** PD is associated with higher Factor 4 scores: +0.12 log-odds per SD in PPMI (q = 0.03) and +0.19 in PDBP (p = 0.002). Within each sex: women +0.10 (PPMI p = 0.22) / +0.18 (PDBP p = 0.03); men +0.14 (p = 0.09) / +0.20 (p = 0.02). The interaction is null (PPMI p = 0.82, PDBP p = 0.64). Before the fix this association didn't replicate (PDBP unadjusted diagnosis −log10 p 0.2, now 2.9).
- **Technical check (Nov 19, `mofa/scripts/factor_technical_check.py`; output in `ppmi_mofa_outputs_noplasma/technical_check/` on Terra).** Library size (total raw counts) and genes detected (count > 0) per participant, correlated with every factor, and logit(PD) ~ factor + age + male refitted with both added:

  | | ρ genes detected (PPMI / PDBP) | ρ log library size | PD log-odds per SD, unadjusted → adjusted (PPMI; PDBP) |
  |---|---|---|---|
  | Factor 1 | 0.82 / 0.59 | −0.24 / −0.16 | +0.13 → −0.07 (p = 0.76); −0.02 → +0.18 (p = 0.50) |
  | Factor 4 | 0.48 / 0.30 | 0.06 / −0.01 | +0.12 → +0.05 (p = 0.43); +0.19 → +0.21 (p = 0.001) |
  | Factor 7 | 0.06 / −0.03 | 0.02 / −0.03 | −0.19 → −0.20 (p < 0.001); −0.44 → −0.44 (p < 0.001) |
  | Factor 8 | −0.14 / −0.25 | −0.10 / −0.08 | −0.18 → −0.16 (p = 0.006); −0.11 → −0.10 (p = 0.09) |

  - **Factor 7 is robust:** uncorrelated with either covariate, and its PD association is unchanged.
  - **Factor 4 looks technical.** Its top positive RNA-seq genes are long, mostly neuronal genes with little expression in blood (CSMD2, HMCN2, PAX6, DNAH9, KIF1A, GRIN2B, GRIN2A, LAMA1, SLC6A1, LRP2, DSCAM, NOS1, KCNC1); the negative side is mostly lncRNAs and low-expressed genes (MIR210HG, LINC02156, TOMM22-DT, IER5L-AS1). It tracks genes detected, which points to low-level background reads (pre-mRNA or genomic DNA). Its PPMI PD association disappears after adjustment; in PDBP it doesn't. Its top CSF proteins (IGFBP2, APOD, SPP1, CLU vs PLXNB2, CTSD, HPX, TF) carry little weight, since Factor 4 explains 0.2% of CSF. Don't interpret it as PD biology; RIN or intronic read fraction, not in the current metadata, would settle whether it's genomic DNA.
  - **Caveat:** genes detected also tracks Factor 1 (ρ = 0.82), probably blood cell composition, which PD alters (Craig et al., 2021). Adjusting for it is therefore conservative.
- **Factor 8 replicates in men only, without a sex difference.** PPMI −0.18 log-odds per SD (q = 0.002), PDBP −0.11 (p = 0.06) across both sexes. Within men: −0.23 (q = 0.02) / −0.19 (p = 0.03); within women: −0.15 (p = 0.09) / −0.03 (p = 0.74). The interaction is null (PPMI p = 0.23, PDBP p = 0.19), so this is not evidence of a sex-specific effect.
- **Factor 10 (aging)** has no diagnosis association once age is adjusted for. It's an aging axis, not a PD axis.
- **Clinical severity within PD:** no factor's association with MDS-UPDRS III or MoCA within PD cases replicates, and no severity × sex interaction does either (smallest PPMI q = 0.21). Factor 7 tracks severity within PD **in PPMI only**: higher scores go with milder motor signs (−0.12 SD of MDS-UPDRS III per SD, q = 0.001) and better cognition (+0.25 SD of MoCA, q = 0.005). In PDBP these are null (p = 0.82 and 0.87). PPMI's early, untreated PD may let a severity signal show that PDBP's mixed disease duration and treatment obscure, which matches the Factor 7 caveat above. Factor 7's pooled association with MDS-UPDRS III (−log10 p 9.1 / 9.6) mostly reflects the PD-vs-control difference.
- **MDS-UPDRS III baseline exams:** PPMI had repeated baseline exams for 162 participants, all on PD medication. `query_comprehensive_metadata.py` now keeps one exam per participant (rerun Oct 7): off medication for 137, state not recorded for 19, on medication for 6 (no off exam available). The diagnosis and interaction results were unchanged by this fix, since diagnosis, sex and age didn't differ between duplicate rows.

**Implication for the study question:** in this model, the multi-omic PD signal (Factor 7) looks the same in women and men within the precision available. The single-group sensitivity model below confirms this in a design that keeps average sex differences.

### Single-group sensitivity model (Nov 19, issue 06)

Same inputs as the main model (RNA-seq + CSF, PPMI), trained with one group instead of grouping by sex, so average sex differences are kept. Sex-chromosome genes were kept as well. It converged after 415 iterations, and PDBP was projected from the same inputs as the main model. Factors were matched to the main model with `compare_models.py` (factor-score |r| on the shared participants).

- **A sex factor forms, as expected.** Single-group Factor 3 (4.5% of RNA-seq, 1.0% of CSF) separates the sexes almost completely: mean score +2.2 in women vs −1.9 in men, sex −log10 p > 300 in both cohorts. It has no counterpart in the main model (|r| = 0.02). In PDBP, its unadjusted associations with diagnosis (11.3) and MDS-UPDRS III (13.2) reflect PDBP's male-heavy PD group (486 men / 281 women, vs 208 / 273 among controls); they don't hold once sex is adjusted for.
- **The PD and aging axes carry over.**
  - Main-model Factor 7 ↔ single-group Factor 8 (|r| = 0.94). PD association replicates (log-odds per SD −0.20 PPMI, q < 0.001; −0.50 PDBP, p < 0.001), with similar effects in each sex and a null interaction (PPMI q = 0.77, PDBP p = 0.64).
  - Main-model Factor 4 ↔ single-group Factor 5 (|r| = 0.99). PD association meets the replication rule (+0.12 PPMI, q = 0.045; +0.19 PDBP, p = 0.002), as in the main model, and is probably the same technical axis.
  - Main-model Factor 8 ↔ single-group Factor 7 (|r| = 0.99). As in the main model, it replicates within men only (−0.22, q = 0.02 / −0.18, p = 0.03) with a null interaction.
  - Main-model Factor 10 ↔ single-group Factor 10 (|r| = 0.93): age 30.7 / 36.9. It's now also sex-associated (24.9 / 19.5), because average sex differences are kept.
  - Main-model Factors 1, 3, 5 and 9 match single-group Factors 1, 4, 6 and 9 with |r| = 0.99–1.00; main Factor 6 matches single-group Factor 2 only partly (|r| = 0.80). Matching: `results/mofa/single_group/factor_matching_*.csv`.
- **No replicated sex interaction.** None of the 30 interaction tests replicates. The only PPMI q < 0.05 is the sex factor × MoCA (q = 0.02), which is null in PDBP (p = 0.82) and mostly encodes sex itself. The sex factor's within-men PD effect also meets the replication rule (+0.015 / +0.012 SD), but that factor's SD is set by the gap between the sexes, so within-sex effects on it are negligible.
- **The CSF axis isn't stable across designs.** The main model's CSF factor (Factor 2, 26.5% of CSF) has no match in the single-group model (best |r| = 0.27). The single-group CSF factor (its Factor 2) explains 13.5% of CSF, and CSF total R² falls from 24–27% to 14.4%. It is also sex-associated (15.5 / 12.6), which suggests that proteins with average sex differences pulled the fit toward a different, weaker CSF pattern. This matches the subset-8 drop-out (Section 2). The CSF factor is reproducible within one design but depends on the design. If CSF becomes a focus, test 15 factors or `weight_views=True`.
- **Exploratory only:** single-group Factor 2 (CSF) goes with better MoCA within PD in both cohorts (PPMI +0.22 SD, q = 0.012; PDBP +0.09 SD, p = 0.012). PDBP has CSF for only 107 participants, so for most of PDBP this score rests on very little data.
- **PDBP variance explained (Nov 19 rerun):** 57.5% of RNA-seq and 11.7% of CSF (107 participants); per-factor Spearman agreement with PPMI 0.84 (RNA-seq) and 0.45 (CSF).
- **Subsampling stability:** all 10 factors are stable, each reproduced in 10 of 10 subsets (minimum |r| 0.88–1.00), including the single-group CSF factor (minimum |r| = 0.97). Subsets were drawn from the whole cohort, since there is one group. `results/mofa/single_group/subset_stability_summary.csv`.

**Conclusion on the sex question:** neither MOFA design (factors centered within each sex, or a single group keeping average sex differences) finds a multi-omic PD axis whose effect differs by sex. The replicated PD signatures are the same in women and men in both designs.

## 4. Known issues with the variance-explained replication

The "PPMI vs PDBP R² per factor" comparison (`cohort_variance_replication_comparison.png`) is not interpretable yet. PPMI Factor 1 shows ≈ 120%, and PPMI Factor 10 is negative.

1. **The view names are lost** (`extract_factors.py`). mofapy2 stores view names in a dataset at `views/views`. Because `f['views']` is a group, the script reads its key list (`['views']`), finds one entry where three are expected, and falls back to `View_0`, `View_1`, `View_2`. PPMI view names therefore never match PDBP's `rna_seq`, `csf_proteomics` and `plasma_proteomics`. As a result, `compare_variance_replication` silently switches to the per-factor fallback.
2. **R² values are summed where they shouldn't be** (`evaluate_factors.py`). The fallback sums R² over every row: 3 views × 2 sex groups in PPMI, and 3 views × 1 pooled group in PDBP. Factor 1's View_2 value (≈ 59% per group) counted twice gives ≈ 119%. R² from different views or groups uses different denominators and can't be added. The two cohorts also don't have the same number of rows being summed.
3. **The two cohorts use different R² formulas.** PPMI uses MOFA's $1 - SS_{res}/SS_{tot}$, calculated on data centered within each sex. PDBP (`project_pdbp.py`) uses $\sum \hat{X}_k^2 / SS_{tot}$ on data that is neither centered within each sex nor scaled by view. These are different quantities.
4. **The projection doesn't match how training treated the data:**
   - **No view scaling.** Training used `scale_views=True`, but PDBP isn't scaled the same way, so its views reach the shared weights on a different scale.
   - **Missing values are filled with 0.** For a participant without proteomics, a block of zeros is treated as real data, which pulls their factor scores toward 0. MOFA would instead leave those features out for that participant.
   - **Views aren't weighted by their noise.** All features are stacked into one unweighted least-squares fit, so the 5,000 RNA-seq features dominate the ~450 protein features. MOFA weights each view by its estimated noise precision.

**Fix (implemented and rerun on Terra Oct 5; results in Section 1):**
- **`extract_factors.py`:** reads view names from `views/views` and stops with an error if the counts don't match. It also exports the model's total R² per view and group (`variance_explained_total.csv`).
- **`project_pdbp.py`:**
  - Takes `--metadata` (with participant ID and sex).
  - Applies the same step mofapy2 applies in training: center each feature within each sex, then scale each view to unit variance.
  - Projects each participant using only the features they actually have, weighted by each feature's noise precision. The precision is estimated from the PPMI training residuals stored in the HDF5, because mofapy2 doesn't save Tau by default.
  - Calculates R² per view and group with mofapy2's own formula.
  - Writes the `group` with the factor scores and reports how many model features are actually observed in PDBP.
- **`evaluate_factors.py`:** compares R² per factor side by side for each view and group, never summing. It writes `cohort_variance_replication_comparison.csv` and a grid plot with one panel per view and group, and prints the Spearman agreement of each per-factor R² profile. If the views or groups don't match, it stops with an error instead of falling back silently.
- **Verification:**
  - A small mofapy2 0.7.5 model trained with the project's settings (two sexes with built-in mean differences, missing proteomics views, missing features).
  - Re-projecting the training cohort reproduces MOFA's own factor scores (r = 0.97–1.00) and per-factor R² to within 4 percentage points.
  - A held-out cohort gives matching R² for each view and sex, recovers a simulated female-only factor (female R² 7–12%, male < 1%), and shows no remaining sex association (−log10 p < 0.05).
  - Unit tests covering these cases have been added to `mofa/tests/test_pdbp_projection.py` and `mofa/tests/test_evaluate_factors.py`.
- **Saved Tau and AlphaZ (added Oct 5):** when the model file has them (`train_mofa.py` now saves both), `project_pdbp.py` uses the saved noise precision and each sex's factor-activity precision as the prior. A factor that is inactive in one sex is then shrunk toward zero for that sex, as in training. Older models fall back to residual-estimated noise and an N(0, I) prior, and the log says which path was used. On the toy model, the saved values brought re-projected factor scores slightly closer to MOFA's own (r ≥ 0.97 for all factors, ≥ 0.996 for active ones; R² within 3.7 points vs 3.9 with the fallback).
- **Remaining approximation:** MOFA's own factor update also accounts for uncertainty in the weights (their second moment under the spike-and-slab prior), which mofapy2 doesn't save. This accounts for the remaining small gap.

## 5. Proteomics non-detects (zeros) and how they are now handled

**Diagnostics** (`mofa/scripts/plasma_diagnostics.py`, run on Terra Oct 5; output in `plasma_diagnostics/`). The raw PDIA intensity matrices encode non-detects as exact zeros. The pipeline transformed them with `log2(x + 1)` into 0, about 17 log2 units below a typical protein. After centering they became extreme outliers. Values are for PPMI / PDBP.

| | Plasma | CSF |
|---|---|---|
| Samples with data | 167 / 107 | 354 / 107 |
| Entries that are exact zeros | 14.4% / 7.9% | 0.02% / 0.05% |
| Proteins with > 30% zeros (> 50%) | 45 (13) / 12 (3) | 0 / 0 |
| Worst sample, % zeros | 65% / 21% | 0.8% / 0.8% |
| Share of the processed view's variance coming from zero entries | **63% / 71%** | 16% / 26% |
| Per-sample zero fraction vs plate (−log10 p, PPMI) | 0.2 | 0.1 |
| PCA, zeros kept: PPMI PC1 variance; PDBP replication ratio (top 1 / top 5 PCs) | 13.8%; 0.49 / 0.46 | 22.9%; 0.88 / 0.74 |
| PCA, zeros as missing: PPMI PC1 variance; PDBP replication ratio (top 1 / top 5 PCs) | 16.3%; 0.40 / 0.69 | 27.1%; 0.95 / 0.83 |

The replication ratio is the PDBP variance explained by the top-k PPMI PCs, divided by the variance explained by PDBP's own top-k PCs.

**What this shows:**
- **The plasma view mostly encoded which proteins were detected, not how abundant they were.** Zeros made up 63–71% of its variance, and the main proteins driving PPMI PC1 were the zero-heavy ones (34–39% zeros). The zeros aren't a plate effect. They differ in rate between cohorts (14% vs 8%), so a zero pattern can't transfer from PPMI to PDBP. This explains the plasma result in Section 1.
- **With zeros treated as missing, plasma has real abundance structure** (PC1 16%, driven by low-zero proteins), and replication across the top 3–5 PCs improves to a ratio of 0.69. Under either treatment, no plasma PC associates with age, diagnosis, sex or any current MOFA factor (all −log10 p ≤ 2.4). Plasma may therefore end up as a plasma-only signal unrelated to the clinical question.
- **CSF has the same artifact on a smaller scale.** Zeros are rare, but each one is an extreme outlier, together making up 16–26% of CSF variance. Treating them as missing strengthens CSF PC1 (23% → 27–30%) and its replication (0.88 → 0.95). CSF PC1 corresponds to MOFA Factor 2 (−log10 p = 66), which also links to Factor 4 (7.5), sex (5.6) and diagnosis (3.0).

**What changed in preprocessing (implemented Oct 5; PPMI preprocessing and training need rerunning):**
- **`proteomics_qc.normalize_proteomics(..., zeros_as_missing=True)`:** marks exact zeros as missing before `log2(x + 1)`. MOFA treats missing values natively, and the projection uses only observed features (Section 4).
- **`proteomics_qc.filter_proteomics(...)`:** a new QC step on the raw matrix. It drops samples with more than 50% zeros, then proteins with more than 30% zeros among the remaining samples. Both thresholds are adjustable. Notebook cell 25 now applies it to PPMI CSF and plasma before normalizing.
- **`preprocess_pdbp.py`:** now uses the same `proteomics_qc` functions, so the two cohorts can no longer drift apart.
  - It applies the same sample QC (`--max_sample_zero_frac`, default 0.50) and treats zeros as missing (`--keep_zeros` turns this off).
  - It does **not** filter proteins. The protein set comes from the PPMI-trained model, so PDBP plays no part in feature selection.
- **Tests:** added to `mofa/tests/test_pdbp_preprocessing.py`.

**Expected effect:**
- **Plasma:** about 45 of 202 proteins and the worst PPMI samples will be dropped.
- **CSF:** no proteins or samples are dropped. Its few zeros become missing.

**Caveat:** non-detects are missing because the protein was low, not at random, which MOFA's model of missingness doesn't account for. Filtering out zero-heavy proteins keeps the remaining missingness small. As a sensitivity check, impute zeros at half each protein's minimum detected value and confirm the main factors don't change.

**Decision (Oct 6): drop plasma from the main model.** The rule set before retraining was to keep plasma only if it gained a factor that replicates in PDBP. In the retrained model, no factor explains more than 0.7% of plasma variance in PPMI, the plasma total is 1.2–1.5% (lower than before the zero fix), and the PDBP total is below zero (−0.5 / −0.2%). Plasma does have structure of its own (PC1 16% once zeros are masked), but it isn't shared with RNA-seq or CSF and doesn't track any clinical variable.
- **Main model:** RNA-seq + CSF, built with `validate_and_intersect_views(..., plasma_df=None)` into `ppmi_mofa_inputs_noplasma/` and trained into `ppmi_mofa_outputs_noplasma/` (new notebook cells after the three-view model). PDBP goes into `pdbp_mofa_inputs_noplasma/` and `pdbp_mofa_outputs_noplasma/`.
- **Sensitivity analysis:** the three-view model (`ppmi_mofa_outputs/`).
- **Limitation to report:** plasma proteomics (159 proteins; 165 PPMI and 107 PDBP participants) showed no shared or clinically associated structure.

## 6. Next steps

- [x] Fix the R² extraction and comparison, and the projection preprocessing (Section 4).
- [x] Rerun `extract_factors.py`, `project_pdbp.py` (with `--metadata data/metadata_m0_PDBP.csv`) and `evaluate_factors.py` on Terra, regenerate the replication figure, and confirm the view order (Oct 5).
- [x] Report per-sex factor activity (R² for each sex group) for PPMI and the re-projected PDBP (Section 1). No factor is clearly specific to one sex in the current model; recheck after retraining.
- [x] Diagnose plasma (`plasma_diagnostics.py`, Oct 5): non-detect zeros dominated the view (Section 5).
- [x] Treat proteomics zeros as missing and add zero-based QC (samples > 50%, proteins > 30%) in PPMI and PDBP preprocessing (Section 5).
- [x] Rerun PPMI proteomics preprocessing and `preprocess_pdbp.py` (Oct 6). PPMI plasma: 2 of 167 samples and 43 of 202 proteins dropped; CSF: none. PDBP: no samples dropped.
- [x] Retrain on the new proteomics inputs (Oct 6): all 5 seeds converged at 220 iterations, with identical solutions.
- [x] Decide on plasma (Oct 6): drop it from the main model; keep the three-view model as a sensitivity analysis (Section 5).
- [x] Fix the stability rule in `subsample_stability.py` (Oct 6): stability is reproducibility alone, and R² is reported as a separate `active` flag. Also reports each factor's weakest subset and what it matched there.
- [x] Make plasma optional in `prep_utils.validate_and_intersect_views` and add no-plasma cells to the notebook (Oct 6).
- [x] On Terra, rerun `subsample_stability.py` on the three-view model with its existing `--output`. It reuses the 10 saved subset models and only updates the summary; check which factor replaced Factor 2 in its weakest subset.
- [x] (Oct 7) On Terra, build and train the no-plasma main model, then extract, preprocess and project PDBP, evaluate both cohorts, and run subsampling (commands in the notebook's raw cell after the no-plasma training). Factor numbers may change again, so re-identify the PD and aging axes.
- [ ] Optional: half-minimum imputation sensitivity check for CSF zeros.
- [ ] Replace the all-factor Spearman summary with one restricted to active factors (for example PPMI R² ≥ 1% in that view and sex).
- [x] Update `train_mofa.py` (Oct 5): 1,000 iterations by default (`n_iter`), each HDF5 now also saves `Tau` and `AlphaZ`, and the log says whether a seed converged or stopped at the cap.
- [x] Extend `select_best_model.py` (Oct 5): reports convergence per seed (converged vs stopped at the cap, last ELBO change) and seed stability (best-matching factor in every other seed by |r| of factor scores). It writes `seed_convergence_summary.csv` and `seed_factor_matching*.csv` next to the models.
- [x] **Seed restarts are not independent.** mofapy2 0.7.5 initializes factors from a PCA of the data and weights at zero, so the seed only matters through scikit-learn's randomized PCA. The multi-seed ELBO check in Task 1.4 is therefore weak evidence of stability. Add a subsampling check: retrain on several random 80% subsets of PPMI participants, then match factors back to the full model.
- [x] Run `mofa/scripts/subsample_stability.py` on the three-view model (Oct 6): 10 subsets of 80% of participants, drawn within each sex. Results in Section 2; Factors 1–5 and 7–10 are stable, Factor 6 isn't. Repeat for the no-plasma main model, and carry forward only stable factors.
- [x] Deduplicate the PPMI metadata merge (1,329 participants become 1,432 rows): `evaluate_factors.dedupe_metadata` (Oct 7), to be rerun on Terra.
- [x] **Sex × diagnosis interaction test** (`clinical_models.py`, run Oct 7, PD vs Ctrl, rerun with the MDS-UPDRS III exam fix and deduplicated metadata): no replicated interaction; Factor 7 is the replicated, sex-invariant PD axis (Section 3). Ticket `pdbp-projection` 04 closed.
- [x] **Single-group sensitivity check** (Nov 19; Section 3): sex-chromosome genes kept and a sex factor formed; the PD and aging axes carried over; no replicated interaction. Subsampling: all 10 factors stable.
- [x] **PDBP RNA-seq library-size normalization** (Nov 19; Section 0): PDBP reprocessed, projected and re-evaluated for the main and single-group models. Factor 4 now replicates as a PD axis; no change to the sex conclusion.
- [ ] Optional: rerun the three-view model's PDBP projection and evaluation on the normalized PDBP inputs, for consistency.
- [x] Factor 4 technical check (Nov 19; Section 3): probably technical (long neuronal genes, tracks genes detected, PPMI PD association lost after adjustment). Factor 7 is robust to the same adjustment.
- [ ] Optional: get RNA quality metrics (RIN, intronic or intergenic read fraction) from AMP-PD to confirm what Factor 4 captures.
- [ ] Revise the spec's replication metric (3) to cover per-sex factor activity and the interaction test rather than sex-specific loadings.
- [x] **Pathway enrichment, GSEA part (Nov 19; issue 05):** Factor 7 contrasts neutrophil granule genes (PD side) with B-cell and immunoglobulin genes, so it largely tracks blood cell composition; Factor 10 has little pathway structure (3 gene sets), with cytotoxic T/NK genes vs naive-lymphocyte genes among its top weights. Details in `docs/mofa_writeup.md` and `results/mofa/main_noplasma/enrichment/`. Remaining: optional STRING network views.
- [ ] (original item) Then move on to pathway enrichment (issue 05) for stable, replicated factors: the aging and PD axes (Factors 10 and 7 in the main model; not Factor 4, which is probably technical) and any sex-differential PD factors that agree across both models. Use rank-based enrichment on the full weight vector (for example GSEA), since top-5% feature overlap across subsets is only moderate (Section 2).
