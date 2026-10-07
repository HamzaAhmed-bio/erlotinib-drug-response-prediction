# Erlotinib response prediction — Project 2

Completed October 6, 2026. Prepared for Hamza Ahmed with AI assistance.

Random Forest was selected using development cross-validation. Its held-out RMSE was 0.0419, compared with 0.0543 for the training-mean baseline and 0.0499 for tissue-only ridge. The expression-minus-tissue RMSE difference was -0.0080 (95% conditional bootstrap interval -0.0119 to -0.0039); this interval excludes zero.

## Study at a glance
- 701 eligible cell lines from 26 tissues; 560 development and 141 test.
- Uniform 0.002001–2 micromolar screening range; GDSC2 release 8.5 and frozen 2021 FPKM expression.
- 37,263 input gene IDs; 1,000 high-variance genes selected inside each training fold.
- Four models: training mean, tissue-only ridge, Elastic Net, Random Forest.
- Primary expression family chosen using five-fold development CV, not test performance.
- Held-out R-squared for Random Forest: 0.404; Spearman correlation: 0.545.

## Interpretation
The model still compresses the response range: the test cell line with the lowest observed AUC had AUC 0.703, but Random Forest predicted 0.889. Improved aggregate error does not mean accurate prediction for every sensitive cell line.

This is internal validation on laboratory cell models. It is not external validation, a patient prediction tool, or a validated biomarker. Response AUC is not ROC AUC. Bootstrap intervals are conditional on this one split and fitted models. Inspect errors in unusually sensitive models as well as aggregate RMSE.

Open `Research_Report.html` for the readable notebook export, or `Erlotinib_Response_Study.ipynb` for code and saved outputs. `docs/learning_guide.md` explains the methods in plain language.
