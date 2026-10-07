# Erlotinib response prediction: frozen analysis plan

Prepared 2026-10-06 for Hamza Ahmed with AI assistance, before model fitting.

## Question and scope
Can pretreatment gene expression predict erlotinib drug-response AUC in held-out cancer cell lines, and improve on mean-response and tissue-only baselines? This is an internal, pan-cancer prediction benchmark, not patient treatment guidance or external clinical validation. Erlotinib was selected for its EGFR target and relevance to a lung-cancer research portfolio, before inspecting model performance. No other drugs will be screened for a more favorable result.

## Inputs and eligibility
- GDSC2 release 8.5, 27 October 2023 fitted response workbook: erlotinib, DRUG_ID 1168.
- Cell Model Passports merged Sanger/Broad RNA-seq FPKM archive, 29 March 2021.
- Cell Model Passports model identities/tissues/related-model metadata, 20 November 2024.
- Use one observation per Sanger model ID. Duplicate responses or ambiguous metadata joins must stop preparation for review.
- AUC is the fitted dose-response area, not a classifier ROC AUC. Lower AUC indicates greater inhibition over the screened concentration range. AUC depends on that range.
- Initial source audit found 955 responses: 705 at MIN_CONC 0.002001 and MAX_CONC 2 micromolar. Restrict the primary cohort to this largest exact concentration window; exclude other windows without considering performance. Require finite AUC in [0,1], finite fit RMSE <=0.3, matching expression and metadata, a cell-line model, and less than 10% missing expression.
- Keep unique source gene IDs, with symbols for readability. Do not merge distinct gene IDs merely because their symbols agree. Reject negative expression or infinite values. Use log2(FPKM+1), not a second log of already transformed data. The archive member explicitly contains FPKM.
- Expression processing predates modeling and may include cohort-wide source normalization. We do not claim raw-read independence. The mixed-source archive and old annotations limit transferability.

## Split and preprocessing
- Connect models sharing a public donor/patient ID, sample ID, or parent-model link; keep all connected models in the same group. Missing donor IDs use available relationships, otherwise model ID. Grouping reduces known related-model leakage but cannot exclude undocumented relationships.
- Deterministic 80/20 GroupShuffleSplit, random_state=42, by connected group. Freeze the manifest before model fitting. Sort model IDs first. Do not choose the split on its performance.
- Five-fold shuffled GroupKFold on development data only, random_state=42. Save each model's fold assignment and verify disjoint groups.
- All missing-value imputation (training medians), variance selection and standardization are fitted inside each training fold. Remove training-all-missing and constant genes; select the 1,000 largest-variance log-expression genes using only the fold's training data.
- Test labels are not used in feature selection, fitting, tuning, model-family selection, or reporting choices. Cohort QC inspects validity and availability, not predictive associations.

## Models and prespecified tuning
1. Training-mean prediction.
2. Tissue-only ridge regression: one-hot tissue, unknown tissues mapped to all-zero; alpha in [0.1,1,10,100].
3. Expression Elastic Net: standardized selected genes, alpha in [0.003,0.01,0.03,0.1], l1_ratio in [0.1,0.5,0.9], max_iter=20000, tol=0.0001.
4. Expression Random Forest: 300 trees, max_features in [0.1,0.5], min_samples_leaf in [2,5,10], random_state=42.

Clip every prediction to [0,1], including in cross-validation scoring. Choose tuning settings by mean five-fold validation RMSE. Designate the expression model with lower development CV RMSE as primary before accessing test outcomes. Keep results for every prespecified model. If convergence fails, resolve numerical issues without using test results; document any change.

## Evaluation and uncertainty
- Primary metric: held-out RMSE in drug-response AUC units. Also MAE, R-squared, and Spearman correlation. Undefined constant-prediction correlations are missing, not zero.
- Report development CV means and fold variability as model-selection diagnostics, not independent estimates of the chosen model's generalization performance.
- Bootstrap held-out donor/model groups 2,000 times (seed 2026), retaining all models in a selected group. Report percentile 95% intervals for RMSE and paired RMSE differences (expression minus each baseline). Intervals are conditional on one split and fitted model and do not include training or selection uncertainty.
- Display observed versus predicted AUC, model-error comparison with intervals, residuals, response distribution, and tissue composition. Describe regression toward the mean and skewed response distributions when present.
- Report exploratory lung-only held-out metrics only if at least 10 test models, with sample size and no clinical claim. No test-informed retraining.
- Show Elastic Net coefficients as exploratory model weights, not validated biomarkers or causal effects; correlated genes can substitute for one another.
- Preserve unsuccessful or negative results. No accuracy percentage or clinical response prediction claim.

## Reproducibility and deliverables
Record source URLs, dates, archive CRC and local SHA-256 checksums; save cohort exclusions, frozen split, CV scores, chosen parameters, test predictions, metrics, model coefficients, figures, and environment versions. Deliver executable Python code, an executed notebook, an HTML notebook preview, plain-language summary and learning guide. A second independent dataset is a future extension; this project does not claim external validation.
