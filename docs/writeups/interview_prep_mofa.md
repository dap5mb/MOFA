# Interview prep: sex differences in PD with multi-group MOFA (PPMI → PDBP)

Talking points, not a script. Numbers are from `mofa_writeup.md` and the validation report (as of 2026-11-19). Sections marked **[you]** need your own words about what you owned.

## 2-minute version

Women and men differ in how often PD occurs, when it starts and how it presents. I asked whether that has a molecular counterpart in blood and CSF of living patients. Most molecular sex-difference work in PD comes from post-mortem brain. I used multi-group MOFA on baseline whole-blood RNA-seq and CSF proteomics in PPMI (1,329 people), grouped by sex, and projected the factors into PDBP (1,248 people) for replication.

The interesting part was that the method could not answer the question as first posed. Grouping by sex centers every feature within each sex, so it removes average sex differences by design. I reframed the question as "does any factor relate to PD differently in women and men?" and tested it with sex × diagnosis interaction models in both cohorts, plus a second MOFA design without sex grouping.

Along the way I found and fixed several problems that would have produced false positives: a variance-explained calculation that added up incompatible quantities, a projection that didn't match training, proteomics non-detects that dominated the plasma signal, and a missing library-size normalization in the validation cohort.

The result was a well-bounded negative. One PD factor replicated, with the same effect in women and men. Enrichment showed it is mostly the neutrophil vs B-cell balance in blood. No sex interaction replicated in either design, and the tests had about 80% power for differences of 0.3 SD or more.

## 1. The question and why it mattered

- Sex differences in PD incidence (about 1.5× in men), age at onset and presentation are well documented clinically (Wooten 2004; Haaxma 2007).
- Molecular sex differences come mostly from post-mortem brain transcriptomics (López-Cerdán 2022; Tranchevent 2023). Whether they show up in accessible tissues is open, and that matters for biomarkers and stratified trials.
- Multi-omics, because a shared axis across blood and CSF would be more robust than any single assay, and MOFA can find axes active in one group only.

## 2. My role and the decisions I owned **[you]**

Fill in: what you designed vs inherited, compute environment (Terra, AMP-PD data use agreement), who you reported to. Decisions to claim if they were yours:
- Keeping sex as the MOFA grouping, and adding the interaction test and single-group model once the limitation was clear.
- The replication rule, fixed before looking at PDBP: PPMI q < 0.05 (BH across factors), PDBP p < 0.05, same sign.
- Dropping plasma by a rule set before retraining (keep only if it gains a factor that replicates).
- Replacing seed restarts with subsampling stability, and defining "stable" as reproducibility alone.
- Not interpreting Factor 4 after the technical check, even though it met the replication rule.

## 3. Analytical strategy and why

- **Unsupervised factor model, not a supervised one.** Supervised multi-omics (e.g. DIABLO / sPLS-DA on PD vs control) would find a PD axis by construction and overfit, especially with a few hundred CSF samples. I wanted the main axes of variation first, then to ask which relate to PD and whether that differs by sex.
- **Multi-group MOFA with groups = sex.** Weights are shared across groups and factor activity is per group, so a factor can be present in one sex only. That is the "sex-specific axis" question.
- **Discovery → projection → replication.** Train in PPMI only, hold the weights fixed and project PDBP. No PDBP data are used to select features or fit anything.
- **Interaction models as the primary readout:** factor ~ PD × sex + age; MDS-UPDRS III and MoCA × sex within PD; plus the PD effect within each sex.

## 4. The hard or ambiguous part (pick one or two to tell in depth)

**a. The method couldn't test the question as stated.** mofapy2 centers features within each group by default, so factor means cannot differ by sex in PPMI. The sex ANOVA there is uninformative by construction, and the PDBP "sex associations" were artifacts of a projection that didn't center within sex. Grouping by diagnosis would be worse: it removes the PD vs control difference, which is the disease signal. Resolution: keep groups = sex for activity, add interaction tests (valid because centering within sex keeps PD-vs-control differences inside each sex), and add a single-group model that keeps average sex differences.

**b. Variance explained looked like 120%.** View names were lost when reading the HDF5, so the comparison fell back to summing R² across 3 views × 2 groups, which have different denominators. PPMI and PDBP also used different R² formulas. The projection didn't scale views, filled missing proteomics with zeros (pulling scores to 0) and didn't weight views by noise. I rewrote the projection as the posterior mean z = (WᵀTW + A)⁻¹WᵀTy over each person's observed features, using the model's saved noise precision (Tau) and factor prior (AlphaZ), with mofapy2's own R² formula in both cohorts.

**c. Non-detects in proteomics.** Exact zeros mark undetected proteins. After log transformation each is an outlier about 17 log2 units low. They made up 63–71% of plasma variance, and plasma PC1 was just "which proteins were detected". Setting zeros to missing fixed CSF (Factor 2 rose to 24–27% of CSF variance and replicated) but left plasma with no replicating structure, so plasma was dropped.

**d. The late surprise: library size.** Near the end I found PDBP RNA-seq had been log-transformed without size-factor normalization. After fixing it, a new factor (Factor 4) met the replication rule. Before accepting it, I checked it against technical covariates. Its top genes are long neuronal genes that blood barely expresses (CSMD2, GRIN2B, DSCAM). It tracked genes detected per sample (ρ = 0.48). Its PPMI PD association vanished after adjustment (+0.12 → +0.05, p = 0.43). I concluded it is probably technical background (pre-mRNA or genomic DNA) and did not interpret it. Lesson: passing a replication rule is not enough for a 1%-variance factor, because a technical difference between cases and controls replicates if both cohorts are processed the same way.

## 5. Alternatives and tradeoffs

| Choice | Alternatives | Why I chose this |
|---|---|---|
| groups = sex | groups = diagnosis; 4 groups (sex × dx); single group | Diagnosis grouping removes the disease signal; 4 groups removes both and leaves tiny CSF groups; single group kept as a sensitivity model |
| Interaction tests on factors | Gene-level sex × dx models (limma / DESeq2) | Factors give a few well-powered tests across omics; gene-level is better for small effects, which I note as the next step |
| Zeros as missing | Keep zeros; half-minimum imputation | Zeros dominated variance; missing is closer to MOFA's model, but non-detects are missing-not-at-random (low abundance). Half-min imputation listed as a sensitivity check |
| Drop plasma | Keep it as a third view | No replicating structure after QC (max 0.7% of plasma variance per factor; PDBP R² < 0). The three-view model is kept as a sensitivity analysis and changes nothing |
| 10 factors, default view weighting | 15 factors; `weight_views=True` | CSF factor depends on design and dropped out in one subsample, likely because 5,000 genes dominate the PCA initialization over 249 proteins. Listed as the first thing to try if CSF becomes the focus |
| Subsampling stability | Seed restarts | mofapy2 initializes from PCA, so 5 seeds gave identical solutions (|r| = 1.0); seeds say nothing about stability |
| GSEA on full weights | Over-representation on top 5% | Top-5% overlap across subsamples was only 0.37 for the PD factor; rank-based uses the whole vector. STRING on top genes as a secondary view |

## 6. How I validated it

- **Replication in an independent cohort** with the weights fixed and a rule set in advance.
- **Stability:** retrained on 10 subsets of 80% within sex; stable = |r| ≥ 0.8 in ≥ 9/10. Nine of ten factors stable; the PD factor reproduced in 10/10 (|r| ≥ 0.88). Single-group model: 10/10 factors stable.
- **Two model designs** (sex-grouped and single-group) matched by factor-score correlation; the PD factor carried over (|r| = 0.94).
- **Technical checks:** plate not associated with any factor after batch correction; library size and genes detected for every factor (the PD factor was unaffected, |ρ| ≤ 0.06).
- **Power:** SE of the PD × sex term about 0.11–0.12 SD in each cohort, so about 80% power for 0.3 SD differences.
- **Code checks:** unit tests for each script, including simulated data with a planted female-only factor that the projection recovered (female R² 7–12%, male < 1%).
- **References** retrieved from PubMed and checked automatically (record, DOI, abstract quote).

## 7. Conclusion and impact

- **No replicated sex-specific multi-omic PD axis** in either design. The replicated PD factor (Factor 7: log-odds per SD −0.19 in PPMI, q = 0.001; −0.44 in PDBP, p < 0.001) has the same effect in women and men (interaction p = 0.62 / 0.66).
- **What Factor 7 is:** neutrophil granule genes on the PD side vs B-cell and immunoglobulin genes on the other (156 gene sets at FDR < 0.05; an 18-gene neutrophil cluster in STRING). That matches the higher neutrophil and lower lymphocyte signal reported in PD whole blood in these cohorts (Craig 2021), and suggests it is mostly cell composition.
- **Scope of the negative:** excludes large sex-specific axes that dominate blood/CSF variation, not sex differences in PD biology. Effects smaller than 0.3 SD, or confined to a few genes, would not form a factor.
- **Decision it supported:** close MOFA for the sex question; next steps are gene-level interaction tests, adjusting for cell composition and for disease duration / medication (PDBP's effect is twice PPMI's, which enrolls recently diagnosed PD).

## Likely follow-ups and short answers

- **Why not just run sex-stratified differential expression?** It would be the right tool for small, gene-level effects. MOFA asks a different question, about coordinated axes across omics; I'd frame them as complementary and recommend the gene-level test next.
- **Isn't a negative result just low power?** Partly bounded: I can say what size I could detect (≥ 0.3 SD). PPMI also has more men, so "significant in men only" (Factors 7 and 8) is a power difference; the interaction term is the correct test and is null.
- **Why is the PDBP effect twice as large?** PPMI is recently diagnosed; PDBP has longer duration and treatment. Untested; adjusting for duration and medication would check it.
- **How do you know Factor 7 isn't technical too?** Uncorrelated with library size and genes detected; association unchanged after adjustment; stable in 10/10 subsets; coherent biology (granule-protein cluster). It could still be composition rather than cell-intrinsic change.
- **Missing data in MOFA?** Participants with any view are included; missing views are left out of the likelihood, not imputed. The projection does the same per participant.
- **Why is CSF so different?** Only 354 PPMI and 107 PDBP participants. Factor 2 (24–27% of CSF variance) replicates but depends on design: it dropped out in one subsample and was not recovered by the single-group model.
- **What would you do differently?** Check normalization parity between cohorts and non-detect handling before modeling; define technical-covariate checks as part of the replication rule; plan cell-composition estimates from the start.
- **Sex-chromosome genes?** Kept, consistent with not regressing out biology before MOFA. In the single-group model they formed a pure sex factor (Factor 3) with no counterpart in the grouped model, which acted as a positive control.
