"""Independent numerical, provenance, identity, and train/test checks."""
import hashlib
import json
import os
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from modeling import predict

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get('DRUG_DATA_DIR', ROOT/'data'))
OUT = ROOT/'results'

def main():
    checks=[]
    for record in json.loads((ROOT/'metadata/source_receipt.json').read_text())['files']:
        with (DATA/record['file']).open('rb') as f: digest=hashlib.file_digest(f,'sha256').hexdigest()
        assert digest==record['sha256'], record['file']
    checks.append('All four downloaded source SHA-256 hashes match the receipt')
    manifest=pd.read_csv(OUT/'split_manifest.csv')
    train=manifest.split.eq('development').to_numpy(); test=~train
    assert manifest.model_id.is_unique
    assert not set(manifest.loc[train,'group_id']) & set(manifest.loc[test,'group_id'])
    assert (manifest.loc[test,'cv_fold']==-1).all()
    for f in range(5):
        a=manifest[train & manifest.cv_fold.ne(f)]; b=manifest[train & manifest.cv_fold.eq(f)]
        assert not set(a.group_id)&set(b.group_id)
    assert manifest.MIN_CONC.eq(.002001).all() and manifest.MAX_CONC.eq(2).all()
    checks.append('Unique model identities; one assay window; test and CV group disjointness')
    loaded=np.load(DATA/'prepared_expression.npz',allow_pickle=False)
    X=loaded['X']
    assert np.array_equal(loaded['model_ids'],manifest.model_id.to_numpy())
    # Independently read 10 source model columns and 10 genes, then compare transformation.
    rng=np.random.default_rng(101)
    cols=rng.choice(len(manifest),10,replace=False)
    rows=rng.choice(X.shape[1],10,replace=False)
    raw=pd.read_csv(DATA/'rnaseq_fpkm_20210329.csv',skiprows=[1,2,3],usecols=['model_id']+manifest.iloc[cols].model_id.tolist()).set_index('model_id')
    values=raw.loc[loaded['gene_ids'][rows],manifest.iloc[cols].model_id].to_numpy(dtype=np.float32).T
    np.testing.assert_allclose(X[np.ix_(cols,rows)],np.log2(values+1),rtol=1e-6)
    checks.append('100 sampled expression values independently match log2(source FPKM+1)')
    pred=pd.read_csv(OUT/'test_predictions.csv')
    assert pred.model_id.tolist()==manifest.loc[test,'model_id'].tolist()
    np.testing.assert_allclose(pred.AUC,manifest.loc[test,'AUC'])
    scores=pd.read_csv(OUT/'test_metrics.csv').set_index('model')
    y=pred.AUC.to_numpy()
    for name in scores.index:
        p=pred[name].to_numpy()
        assert np.isfinite(p).all() and ((p>=0)&(p<=1)).all()
        rmse=np.sqrt(np.mean((y-p)**2)); mae=np.mean(np.abs(y-p)); r2=1-np.sum((y-p)**2)/np.sum((y-y.mean())**2)
        np.testing.assert_allclose([rmse,mae,r2],scores.loc[name,['RMSE','MAE','R2']].to_numpy(dtype=float),rtol=1e-10,atol=1e-12)
        model=joblib.load(ROOT/'models'/f'{name.replace(" ","_")}.joblib')
        inp=manifest.loc[test,['tissue']] if name=='Tissue-only ridge' else np.zeros((test.sum(),1)) if name=='Mean baseline' else X[test]
        np.testing.assert_allclose(p,predict(model,inp),rtol=1e-9,atol=1e-12)
        if name in ['Elastic Net','Random Forest']:
            variance=np.var(X[train].astype(np.float64),axis=0)
            expected=np.argsort(-variance,kind='stable')[:1000]
            np.testing.assert_array_equal(model.named_steps['features'].indices_,expected)
            np.testing.assert_allclose(model.named_steps['features'].medians_,np.median(X[train][:,expected],axis=0))
            if name=='Elastic Net':
                np.testing.assert_allclose(model.named_steps['scale'].mean_,X[train][:,expected].astype(np.float64).mean(axis=0))
    np.testing.assert_allclose(pred['Mean baseline'],manifest.loc[train,'AUC'].mean())
    checks.append('Every reported RMSE, MAE, R2 independently recomputed; all saved model predictions reproduced')
    checks.append('Feature selection, imputation, standardization and mean baseline match development-only statistics')
    cv=pd.read_csv(OUT/'cv_all_scores.csv')
    assert cv.RMSE.notna().all()
    summary=pd.read_csv(OUT/'cv_summary.csv').set_index('model')
    for name,frame in cv.groupby('model'):
        average=frame.groupby('parameters').RMSE.mean()
        assert np.isclose(average.min(),summary.loc[name,'CV_RMSE'])
    selection=json.loads((OUT/'model_selection.json').read_text())
    expected=summary.loc[['Elastic Net','Random Forest'],'CV_RMSE'].idxmin()
    assert expected==selection['selected_expression_model']
    digest=hashlib.sha256((OUT/'split_manifest.csv').read_bytes()).hexdigest()
    assert selection['split_sha256']==digest==json.loads((OUT/'cohort_summary.json').read_text())['split_sha256']
    checks.append('Primary-model choice matches development CV; frozen split hashes reconcile')
    boot=pd.read_csv(OUT/'bootstrap_rmse.csv')
    assert len(boot)==2000
    for name in scores.index:
        np.testing.assert_allclose(np.quantile(boot[name],[.025,.975]),scores.loc[name,['RMSE_CI_low','RMSE_CI_high']].to_numpy(dtype=float))
    for row in pd.read_csv(OUT/'paired_rmse_differences.csv').itertuples():
        np.testing.assert_allclose(np.quantile(boot[row.model]-boot[row.baseline],[.025,.975]),[row.CI_low,row.CI_high])
    checks.append('All bootstrap intervals and paired comparisons match saved replicates')
    audit=pd.read_csv(OUT/'cohort_audit.csv')
    assert audit.eligible.sum()==len(manifest)
    checks.append('Cohort counts reconcile to the source eligibility audit')
    receipt={'status':'passed','checks':checks,'limits':['Expression source spot-check covers 100 values, not every cell','No external or clinical validation','Conditional bootstrap excludes training and model-selection uncertainty'],
             'code_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'scripts').glob('*.py')}}
    (OUT/'validation_receipt.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt,indent=2))

if __name__=='__main__': main()
