"""Create a reviewable notebook companion, execute its checks/plots, export HTML."""
import json
from pathlib import Path
import sys
import os
import nbformat as nbf
from nbclient import NotebookClient
from nbconvert import HTMLExporter
from jupyter_client import KernelManager
from bs4 import BeautifulSoup

ROOT=Path(__file__).resolve().parents[1]
R=ROOT/'results'
RUNTIME=Path(os.environ.get('DRUG_DATA_DIR',ROOT/'data'))/'notebook_runtime'
RUNTIME.mkdir(parents=True,exist_ok=True)
for variable,subdir in [('IPYTHONDIR','ipython'),('JUPYTER_RUNTIME_DIR','runtime'),('JUPYTER_CONFIG_DIR','config')]:
    os.environ[variable]=str(RUNTIME/subdir)

def main():
    import pandas as pd
    c=json.loads((R/'cohort_summary.json').read_text())
    sel=json.loads((R/'model_selection.json').read_text())
    winner=sel['selected_expression_model']
    scores=pd.read_csv(R/'test_metrics.csv').set_index('model')
    diffs=pd.read_csv(R/'paired_rmse_differences.csv')
    row=scores.loc[winner]
    improvement=(1-row.RMSE/scores.loc['Mean baseline','RMSE'])*100
    tissue_diff=diffs[(diffs.model==winner)&(diffs.baseline=='Tissue-only ridge')].iloc[0]
    test_predictions=pd.read_csv(R/'test_predictions.csv')
    sensitive=test_predictions.loc[test_predictions.AUC.idxmin()]
    limitation=f"The model still compresses the response range: the test cell line with the lowest observed AUC had AUC {sensitive.AUC:.3f}, but {winner} predicted {sensitive[winner]:.3f}. Improved aggregate error does not mean accurate prediction for every sensitive cell line."
    direction=('lower' if tissue_diff.RMSE_difference<0 else 'higher')
    uncertainty=('excludes zero' if tissue_diff.CI_high<0 or tissue_diff.CI_low>0 else 'includes zero')
    conclusion=f"{winner} was selected using development cross-validation. Its held-out RMSE was {row.RMSE:.4f}, compared with {scores.loc['Mean baseline','RMSE']:.4f} for the training-mean baseline and {scores.loc['Tissue-only ridge','RMSE']:.4f} for tissue-only ridge. The expression-minus-tissue RMSE difference was {tissue_diff.RMSE_difference:+.4f} (95% conditional bootstrap interval {tissue_diff.CI_low:+.4f} to {tissue_diff.CI_high:+.4f}); this interval {uncertainty}."
    nb=nbf.v4.new_notebook()
    nb.metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'language_info':{'name':'python','version':sys.version.split()[0]}}
    cells=[]
    def md(s): cells.append(nbf.v4.new_markdown_cell(s))
    def code(s): cells.append(nbf.v4.new_code_cell(s))
    md(f'''# Predicting erlotinib response from gene expression
**A reproducible cancer-cell-line benchmark | Project 2 | October 6, 2026**  
Prepared for **Hamza Ahmed**, with AI assistance.

## tl;dr
{conclusion}

The study includes **{c['eligible_models']} cell lines**, with **{c['development_models']} development** and **{c['test_models']} held-out test** models. These are internal research-model results. They do not establish prediction of patient benefit or a validated biomarker.

{limitation}

**How to read this notebook:** the default run independently checks and visualizes saved pipeline outputs. Set `RUN_FULL_ANALYSIS=True` below to regenerate preparation and model fitting from the downloaded source files. The separate `reproduce.ps1` runs the full workflow, including source downloads.''')
    md('''## Context & Methods
**Question:** does pretreatment gene expression improve erlotinib-response prediction beyond a training-mean and a tissue-only baseline?

**Outcome:** fitted GDSC2 drug-response AUC, from 0 to 1. Lower values indicate more inhibition across the tested doses. **This is not ROC AUC.** The primary error metric is predictive RMSE; lower is better.

### Key Assumptions
- The primary cohort uses one assay window: **0.002001–2 micromolar**. Of 955 erlotinib records, 250 used other windows and four of the remaining 705 lacked expression.
- The source FPKM matrix is transformed as log2(FPKM+1). Training-only variance ranking selects 1,000 genes from 37,263 available gene IDs. Missing-value handling and standardization stay inside training folds.
- An 80/20 group split and five development folds are frozen at seed 42. Public donor, sample and parent links define groups. All 701 eligible models happened to belong to separate known groups; undocumented relationships remain possible.
- Elastic Net and Random Forest are compared using development CV only. The chosen family and tuning parameters are saved before evaluating test outcomes. All predictions are clipped to [0,1].
- Source measurements are already normalized and come from a merged Sanger/Broad archive. We did not reprocess raw reads or independently remove all batch effects.

Full prespecified choices: [analysis plan](docs/analysis_plan.md). Source definitions and download URLs: [sources](docs/sources.md).''')
    code('''from pathlib import Path
import json, sys, subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display

ROOT = Path.cwd()
assert (ROOT / 'results' / 'split_manifest.csv').exists(), 'Run from the project directory'
RUN_FULL_ANALYSIS = False
if RUN_FULL_ANALYSIS:
    for script in ['02_prepare.py', '03_model.py', '04_validate.py']:
        subprocess.run([sys.executable, str(ROOT / 'scripts' / script)], check=True)
R = ROOT / 'results'
FIG = R / 'figures'
FIG.mkdir(exist_ok=True)
manifest = pd.read_csv(R / 'split_manifest.csv')
predictions = pd.read_csv(R / 'test_predictions.csv')
metrics = pd.read_csv(R / 'test_metrics.csv')
cv = pd.read_csv(R / 'cv_summary.csv')
differences = pd.read_csv(R / 'paired_rmse_differences.csv')
primary = json.loads((R / 'model_selection.json').read_text())['selected_expression_model']
receipt = json.loads((R / 'validation_receipt.json').read_text())
assert receipt['status'] == 'passed'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
    'axes.spines.right':False,'axes.titleweight':'bold','figure.dpi':120,'savefig.dpi':180,
    'axes.labelcolor':'#25324A','text.color':'#25324A','axes.edgecolor':'#A3AAB7'})
colors = {'Mean baseline':'#8A94A6','Tissue-only ridge':'#D5A33A','Elastic Net':'#2A6496','Random Forest':'#9A5678'}
def finish(fig, filename):
    fig.tight_layout()
    fig.savefig(FIG / (filename + '.png'), bbox_inches='tight')
    fig.savefig(FIG / (filename + '.svg'), bbox_inches='tight')
    plt.show()
print('Saved results loaded; independent validation status:', receipt['status'])''')
    md('''## Data
The available data cover several cancer tissues. A good result on a mixed-tissue split does not demonstrate generalization to a new tissue or another laboratory. Inspect the response distribution: a narrow, resistant majority can make overall error look small while sensitive outliers remain difficult.''')
    code('''counts = pd.crosstab(manifest.tissue, manifest.split)
ordered = counts.sum(axis=1).sort_values(ascending=False).index
top = counts.loc[ordered[:12]].copy()
top.loc['Other tissues'] = counts.loc[ordered[12:]].sum()
fig, ax = plt.subplots(figsize=(10,6))
ypos = np.arange(len(top))
ax.barh(ypos, top.development, label='Development', color='#2A6496')
ax.barh(ypos, top.test, left=top.development, label='Held-out test', color='#D5A33A')
ax.set(yticks=ypos, yticklabels=top.index, xlabel='Number of cell lines', title='Cohort composition by tissue and split')
ax.invert_yaxis(); ax.legend(frameon=False)
finish(fig, '01_tissue_composition')

fig, ax = plt.subplots(figsize=(10,4))
bins = np.linspace(0,1,31)
for split,color,style in [('development','#2A6496','-'),('test','#B27022','--')]:
    values=manifest.loc[manifest.split.eq(split),'AUC']
    ax.hist(values,bins=bins,weights=np.ones(len(values))/len(values),histtype='step',
            linewidth=2,linestyle=style,label=f'{split.title()} (n={len(values)})',color=color)
ax.set(xlabel='Observed drug-response AUC (lower = more inhibition)',ylabel='Fraction of cell lines',
       title='Response distribution within the same assay window',xlim=(0,1))
ax.legend(frameon=False)
finish(fig, '02_response_distribution')
display(manifest.groupby('split').AUC.agg(['count','min','median','mean','max']).round(4))''')
    md('''## Results
### Development selection
These CV scores were used to choose settings and the primary expression model. They are selection diagnostics, not independent test estimates. The fold standard deviation describes variability among these five folds; it is not a confidence interval.''')
    code('''display(cv.round({'CV_RMSE':4,'CV_fold_SD':4}))
print('Primary expression model, selected before test evaluation:', primary)''')
    md('''### Held-out prediction error
The bars below show **test RMSE**, with 95% percentile intervals from 2,000 bootstrap resamples of held-out groups. These intervals condition on the fitted models and split. They do not measure uncertainty from repeating model development. Use the paired differences below, rather than overlap of separate intervals, to assess the comparison.''')
    code('''for row in metrics.itertuples():
    error = predictions.AUC.to_numpy() - predictions[row.model].to_numpy()
    assert np.isclose(np.sqrt(np.mean(error**2)),row.RMSE)
    assert np.isclose(np.mean(np.abs(error)),row.MAE)
fig, ax = plt.subplots(figsize=(10,4.5))
ypos = np.arange(len(metrics))
errs=np.vstack([metrics.RMSE-metrics.RMSE_CI_low,metrics.RMSE_CI_high-metrics.RMSE])
ax.barh(ypos,metrics.RMSE,xerr=errs,color=[colors[x] for x in metrics.model],capsize=4)
ax.set(yticks=ypos,yticklabels=metrics.model,xlabel='Test RMSE in drug-response AUC units (lower is better)',
       title=f'Held-out model comparison | {len(predictions)} cell lines',xlim=(0,None))
ax.invert_yaxis()
finish(fig,'03_test_error')
display(metrics.round(4))
display(differences.round(4))''')
    md(f'''The primary model's error was **{direction} than the tissue-only baseline** in this split. The paired interval **{uncertainty}**. Negative differences favor expression. The final interpretation must preserve this uncertainty; a favorable point estimate alone is not proof of reproducible improvement.''')
    md('''### What the predictions miss
The diagonal is perfect agreement. Points far from it are large errors. The residual is **predicted minus observed** AUC, so a positive residual means the model predicts less inhibition than was observed. Predictions clustered near the resistant majority indicate regression toward the mean.''')
    code('''fig, axes = plt.subplots(1,2,figsize=(12,5))
y = predictions.AUC.to_numpy(); p = predictions[primary].to_numpy()
axes[0].scatter(y,p,s=28,alpha=.7,color=colors[primary],edgecolors='white',linewidths=.3)
low=max(0,min(y.min(),p.min())-.04); high=min(1,max(y.max(),p.max())+.02)
axes[0].plot([low,high],[low,high],'--',color='#636B76',label='Perfect agreement')
axes[0].set(xlim=(low,high),ylim=(low,high),xlabel='Observed drug-response AUC',ylabel='Predicted drug-response AUC',title=primary+' | test predictions')
axes[0].legend(frameon=False)
axes[1].scatter(y,p-y,s=28,alpha=.7,color=colors[primary],edgecolors='white',linewidths=.3)
axes[1].axhline(0,color='#636B76',linestyle='--')
axes[1].set(xlabel='Observed drug-response AUC',ylabel='Prediction error (predicted - observed)',title='Errors across observed sensitivity')
finish(fig,'04_predictions_and_errors')
worst=predictions[['model_id','tissue','AUC',primary]].copy()
worst['absolute_error']=(worst[primary]-worst.AUC).abs()
display(worst.nlargest(5,'absolute_error').round(4))''')
    md('''### Exploratory model weights
These are Elastic Net coefficients per one development-set standard deviation of expression. Positive coefficients predict higher response AUC (less inhibition); negative coefficients predict lower AUC. These are correlated predictive features, not causal effects or independently validated biomarkers. This view is exploratory, even if Random Forest is the primary model.''')
    code('''coef=pd.read_csv(R/'elastic_net_coefficients.csv')
top=coef[coef.absolute_coefficient.gt(0)].head(15).sort_values('coefficient_per_training_SD')
if len(top):
    labels=top.symbol.fillna(top.gene_id).astype(str)+' | '+top.gene_id.astype(str)
    fig,ax=plt.subplots(figsize=(11,6))
    ax.barh(labels,top.coefficient_per_training_SD,color=['#2A6496' if v<0 else '#B27022' for v in top.coefficient_per_training_SD])
    ax.axvline(0,color='#636B76',linewidth=1)
    ax.set(xlabel='AUC coefficient per training-set expression SD',title='Largest Elastic Net weights | exploratory')
    finish(fig,'05_elastic_net_weights')
else:
    print('Regularization retained no nonzero gene coefficients.')
lung=pd.read_csv(R/'lung_exploratory_metrics.csv')
if len(lung):
    print('Exploratory lung subset: small sample; no separate model fitting or clinical inference.')
    display(lung.round(4))
else:
    print('Fewer than 10 lung test models; no prespecified lung-subset metrics reported.')''')
    md(f'''## Takeaways
{conclusion}

- This experiment measures internal generalization to held-out cell lines from the same assembled resource.
- {limitation}
- Gene coefficients are exploratory. Tissue composition, source processing, related models not captured by metadata, one random split and limited sensitive examples remain limitations.
- The next scientific extension is a **preplanned independent cell-line screening dataset** with compatible response endpoints and model identities. It should be kept separate from this internal test set.
- No clinical treatment decisions, patient-response claims, or admission guarantees follow from these results.

### Reproducibility and attribution
The full download/preparation/modeling/validation workflow is in `scripts/` and `reproduce.ps1`. Frozen source hashes, split identities, every tuning-fold score, predictions and bootstrap draws are retained. Raw source files and fitted binary models are excluded from the shareable ZIP. The notebook is a readable companion to those artifacts.

Read the [learning guide](docs/learning_guide.md) before presenting the work. Code, analysis and documentation were prepared with AI assistance; describe your contribution accurately.''')
    code('''print('Validation checks:')
for check in receipt['checks']:
    print(' -',check)
print()
print('Verification limits:')
for limit in receipt['limits']:
    print(' -',limit)''')
    nb.cells=cells
    path=ROOT/'Erlotinib_Response_Study.ipynb'
    nbf.validate(nb); nbf.write(nb,path)
    km=KernelManager(kernel_name='python3')
    km.kernel_spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}']
    client=NotebookClient(nb,km=km,timeout=600,resources={'metadata':{'path':str(ROOT)}})
    try:
        client.execute()
    finally:
        if km.has_kernel: km.shutdown_kernel(now=True)
    nbf.validate(nb); nbf.write(nb,path)
    exporter=HTMLExporter(template_name='lab',exclude_input=True)
    body,_=exporter.from_notebook_node(nb,resources={'metadata':{'name':'Erlotinib Response Study'}})
    html=BeautifulSoup(body,'html.parser')
    html.title.string='Erlotinib Response Study - Hamza Ahmed'
    # All evidence is static; remove scripts so the export is fully offline.
    for script in html.find_all('script'): script.decompose()
    (ROOT/'Research_Report.html').write_text(str(html),encoding='utf-8')
    summary=f'''# Erlotinib response prediction — Project 2

Completed October 6, 2026. Prepared for Hamza Ahmed with AI assistance.

{conclusion}

## Study at a glance
- {c['eligible_models']} eligible cell lines from {c['tissues']} tissues; {c['development_models']} development and {c['test_models']} test.
- Uniform 0.002001–2 micromolar screening range; GDSC2 release 8.5 and frozen 2021 FPKM expression.
- {c['genes']:,} input gene IDs; 1,000 high-variance genes selected inside each training fold.
- Four models: training mean, tissue-only ridge, Elastic Net, Random Forest.
- Primary expression family chosen using five-fold development CV, not test performance.
- Held-out R-squared for {winner}: {row.R2:.3f}; Spearman correlation: {row.Spearman:.3f}.

## Interpretation
{limitation}

This is internal validation on laboratory cell models. It is not external validation, a patient prediction tool, or a validated biomarker. Response AUC is not ROC AUC. Bootstrap intervals are conditional on this one split and fitted models. Inspect errors in unusually sensitive models as well as aggregate RMSE.

Open `Research_Report.html` for the readable notebook export, or `Erlotinib_Response_Study.ipynb` for code and saved outputs. `docs/learning_guide.md` explains the methods in plain language.
'''
    (ROOT/'Research_Summary.md').write_text(summary,encoding='utf-8')
    print('Executed notebook and HTML export saved.',flush=True)

if __name__=='__main__': main()
