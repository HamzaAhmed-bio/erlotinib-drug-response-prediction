"""Tune on development groups, freeze a winner, then evaluate the held-out set."""
import hashlib
import json
import os
from pathlib import Path
import warnings
from datetime import datetime, timezone
import joblib
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import ParameterGrid
from threadpoolctl import threadpool_limits
from modeling import ExpressionFeatures, make_model, predict, metrics
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get('DRUG_DATA_DIR', ROOT/'data'))
OUT = ROOT/'results'
GRID = {'Mean baseline':[{}], 'Tissue-only ridge':{'alpha':[.1,1,10,100]},
        'Elastic Net':{'alpha':[.003,.01,.03,.1],'l1_ratio':[.1,.5,.9]},
        'Random Forest':{'max_features':[.1,.5],'min_samples_leaf':[2,5,10]}}

def main():
    warnings.filterwarnings('error', category=ConvergenceWarning)
    manifest = pd.read_csv(OUT/'split_manifest.csv')
    loaded = np.load(DATA/'prepared_expression.npz',allow_pickle=False)
    X, ids, gene_ids = loaded['X'], loaded['model_ids'],loaded['gene_ids']
    assert np.array_equal(ids,manifest.model_id.to_numpy())
    dev = manifest.split.eq('development').to_numpy()
    test = ~dev
    groups = manifest.group_id.to_numpy()
    assert not set(groups[dev]) & set(groups[test])
    # Hold-out labels are only accessed after the model selection record is written.
    ydev = manifest.loc[dev,'AUC'].to_numpy()
    Xdev = X[dev]
    tissue = manifest[['tissue']].fillna('Unknown')
    folds = manifest.loc[dev,'cv_fold'].to_numpy()
    fold_data = []
    for fold in range(5):
        tr,va = np.flatnonzero(folds!=fold), np.flatnonzero(folds==fold)
        assert not set(groups[dev][tr]) & set(groups[dev][va])
        selector = ExpressionFeatures().fit(Xdev[tr])
        xtr,xva = selector.transform(Xdev[tr]),selector.transform(Xdev[va])
        scaler = StandardScaler().fit(xtr)
        fold_data.append((tr,va,xtr,xva,scaler.transform(xtr),scaler.transform(xva)))
    rows, chosen, cv_summary = [],{},[]
    for name, grid in GRID.items():
        print('Tuning:',name,flush=True)
        candidates = grid if isinstance(grid,list) else list(ParameterGrid(grid))
        candidate_scores = []
        for params in candidates:
            scores=[]
            for fold,(tr,va,xtr,xva,ztr,zva) in enumerate(fold_data):
                pipeline = make_model(name,params)
                if name in ['Elastic Net','Random Forest']:
                    model = clone(pipeline.named_steps['model'])
                    a,b = (ztr,zva) if name=='Elastic Net' else (xtr,xva)
                elif name=='Tissue-only ridge':
                    model=pipeline
                    a,b=tissue.loc[dev].iloc[tr],tissue.loc[dev].iloc[va]
                else:
                    model=pipeline
                    a,b=np.zeros((len(tr),1)),np.zeros((len(va),1))
                model.fit(a,ydev[tr])
                score=metrics(ydev[va],predict(model,b))['RMSE']
                scores.append(score)
                rows.append({'model':name,'parameters':json.dumps(params,sort_keys=True),'fold':fold,'n_validation':len(va),'RMSE':score})
            candidate_scores.append((float(np.mean(scores)),float(np.std(scores,ddof=1)),params))
        best=min(candidate_scores,key=lambda a:a[0])
        chosen[name]=best[2]
        cv_summary.append({'model':name,'CV_RMSE':best[0],'CV_fold_SD':best[1],'parameters':json.dumps(best[2],sort_keys=True)})
        print(name,'best CV RMSE',round(best[0],5),best[2],flush=True)
    cv=pd.DataFrame(cv_summary)
    pd.DataFrame(rows).to_csv(OUT/'cv_all_scores.csv',index=False)
    cv.to_csv(OUT/'cv_summary.csv',index=False)
    winner=cv[cv.model.isin(['Elastic Net','Random Forest'])].sort_values('CV_RMSE').iloc[0].model
    selection={'selected_expression_model':winner,'parameters':chosen,'criterion':'lowest mean development 5-fold grouped CV RMSE',
               'selected_utc':datetime.now(timezone.utc).isoformat(),
               'split_sha256':hashlib.sha256((OUT/'split_manifest.csv').read_bytes()).hexdigest(),
               'test_results_used_for_selection':False}
    (OUT/'model_selection.json').write_text(json.dumps(selection,indent=2))
    print('Frozen primary model:',winner,flush=True)
    predictions=manifest.loc[test,['model_id','group_id','tissue','AUC']].copy()
    predictions['development_mean']=float(ydev.mean())
    models={}
    (ROOT/'models').mkdir(exist_ok=True)
    for name,params in chosen.items():
        model=make_model(name,params)
        if name=='Tissue-only ridge': a,b=tissue.loc[dev],tissue.loc[test]
        elif name=='Mean baseline': a,b=np.zeros((dev.sum(),1)),np.zeros((test.sum(),1))
        else: a,b=X[dev],X[test]
        model.fit(a,ydev)
        predictions[name]=predict(model,b)
        models[name]=model
        joblib.dump(model,ROOT/'models'/f'{name.replace(" ","_")}.joblib',compress=3)
    predictions.to_csv(OUT/'test_predictions.csv',index=False)
    y=predictions.AUC.to_numpy()
    eval_rows=[]
    for name in chosen:
        eval_rows.append({'model':name,'n_test':len(y),**metrics(y,predictions[name].to_numpy())})
    evaluation=pd.DataFrame(eval_rows)
    rng=np.random.default_rng(2026)
    unique_groups=predictions.group_id.unique()
    group_rows={g:np.flatnonzero(predictions.group_id.to_numpy()==g) for g in unique_groups}
    boot=[]
    for i in range(2000):
        ix=np.concatenate([group_rows[g] for g in rng.choice(unique_groups,len(unique_groups),replace=True)])
        row={'replicate':i}
        for name in chosen: row[name]=float(np.sqrt(np.mean((y[ix]-predictions[name].to_numpy()[ix])**2)))
        boot.append(row)
    boot=pd.DataFrame(boot)
    boot.to_csv(OUT/'bootstrap_rmse.csv',index=False)
    for i,row in evaluation.iterrows():
        evaluation.loc[i,'RMSE_CI_low'],evaluation.loc[i,'RMSE_CI_high']=np.quantile(boot[row.model],[.025,.975])
    evaluation.to_csv(OUT/'test_metrics.csv',index=False)
    diffs=[]
    for name in ['Elastic Net','Random Forest']:
        for baseline in ['Mean baseline','Tissue-only ridge']:
            delta=boot[name]-boot[baseline]
            point=float(evaluation.set_index('model').loc[name,'RMSE']-evaluation.set_index('model').loc[baseline,'RMSE'])
            lo,hi=np.quantile(delta,[.025,.975])
            diffs.append({'model':name,'baseline':baseline,'RMSE_difference':point,'CI_low':lo,'CI_high':hi,'interpretation':'negative favors expression model'})
    pd.DataFrame(diffs).to_csv(OUT/'paired_rmse_differences.csv',index=False)
    subgroup=[]
    for label,frame in predictions.groupby('tissue'):
        if label=='Lung' and len(frame)>=10:
            for name in chosen: subgroup.append({'tissue':label,'model':name,'n':len(frame),**metrics(frame.AUC,frame[name])})
    pd.DataFrame(subgroup,columns=['tissue','model','n','RMSE','MAE','R2','Spearman']).to_csv(OUT/'lung_exploratory_metrics.csv',index=False)
    annotation=pd.read_csv(OUT/'gene_annotation.csv').set_index('gene_id')
    en=models['Elastic Net']
    genes=gene_ids[en.named_steps['features'].indices_]
    coef=pd.DataFrame({'gene_id':genes,'coefficient_per_training_SD':en.named_steps['model'].coef_})
    coef['symbol']=coef.gene_id.map(annotation.symbol)
    coef['absolute_coefficient']=coef.coefficient_per_training_SD.abs()
    coef.sort_values('absolute_coefficient',ascending=False).to_csv(OUT/'elastic_net_coefficients.csv',index=False)
    summary={'selected_expression_model':winner,'test_models':len(predictions),'test_groups':len(unique_groups),
        'bootstrap_replicates':2000,'convergence_warnings':0,'all_model_predictions_finite':bool(np.isfinite(predictions[list(chosen)].to_numpy()).all()),
        'elastic_net_nonzero_coefficients':int((coef.absolute_coefficient>0).sum()),
        'no_external_validation':True,'model_fits_use_development_only':True}
    (OUT/'modeling_summary.json').write_text(json.dumps(summary,indent=2))
    print(evaluation.to_string(index=False),flush=True)
    print(pd.DataFrame(diffs).to_string(index=False),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=2): main()
