# Explain this project in your own words

## The research question
We asked whether gene expression measured before treatment helps predict how cancer cell lines respond to erlotinib. A cell line is a laboratory model. It is not a patient, and performance here does not establish clinical usefulness.

## What the measurements mean
- **FPKM:** a normalized RNA expression measure. We apply log2(FPKM+1) to compress its large range.
- **Drug-response AUC:** area under a fitted viability-versus-dose curve. Lower values mean more inhibition across the tested doses. It is not classification ROC AUC, and higher is not automatically better.
- **RMSE:** the typical size of prediction errors, emphasizing large errors. Lower is better, and it has the same units as the response AUC.
- **R-squared:** comparison with the test-set-mean reference, which is different from our training-mean baseline. It can be negative.
- **Spearman correlation:** whether models rank sensitive and resistant cell lines similarly to the observations. A constant baseline has no meaningful rank correlation.

## Why four models?
The training-mean baseline predicts the same response for everyone. The tissue-only baseline asks how far we can get just by knowing where the cancer originated. Elastic Net fits a regularized linear combination of gene expression. Random Forest averages many decision trees and can represent nonlinear patterns.

## How we avoided common mistakes
We restricted all responses to the same tested concentration window. We froze a test set before training. Related model groups stay together. Gene selection, imputation and scaling learn only from development/training-fold data. We chose the expression model using development cross-validation, then evaluated it on the held-out data. We kept the results even if they were disappointing.

## What the uncertainty means
The bootstrap repeatedly resamples held-out research-model groups. Its interval describes uncertainty conditional on these trained models and this one split. It does not include all uncertainty from training, changing datasets, or selecting a model. A negative expression-minus-baseline RMSE difference favors expression; an interval crossing zero means the direction remains uncertain by this analysis.

## Questions to practice
1. Why would mixing assay concentration ranges make AUC comparisons misleading?
2. Why is selecting genes before splitting the data a leakage risk?
3. What does the tissue-only baseline protect us from overlooking?
4. Why might a model have low overall RMSE but miss unusually sensitive cell lines?
5. What do your actual test metrics and paired intervals show?
6. What would independent validation require?

## Honest attribution
This research portfolio project was prepared with AI assistance for code, analysis and documentation. Before describing it in an interview, inspect the source, rerun it, and explain the methodological decisions and results yourself. Describe your own contribution accurately. No publication, wet-lab experimentation, clinical validation or independent biomarker discovery is implied.
