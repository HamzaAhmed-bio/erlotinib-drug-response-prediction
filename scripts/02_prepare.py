"""Audit erlotinib eligibility, align FPKM, and freeze a group-disjoint split."""
from pathlib import Path
import hashlib
import json
import os
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit, GroupKFold

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get('DRUG_DATA_DIR', ROOT/'data'))
OUT = ROOT/'results'
OUT.mkdir(exist_ok=True)

def main():
    response_path = ROOT/'metadata/erlotinib_response.csv'
    # Always derive responses from the source workbook, not a hand-edited subset.
    drug = pd.read_excel(DATA/'GDSC2_fitted_dose_response_27Oct23.xlsx')
    drug = drug[(drug.DRUG_ID == 1168) & (drug.DRUG_NAME == 'Erlotinib')].copy()
    drug.to_csv(response_path, index=False)
    assert len(drug) > 0 and not drug.SANGER_MODEL_ID.duplicated().any()
    meta = pd.read_csv(DATA/'model_list_20241120.csv', low_memory=False)
    assert not meta.model_id.duplicated().any()
    keep_cols = ['model_id','model_name','patient_id','sample_id','parent_id','tissue','cancer_type','model_type']
    # The source has three metadata rows and a fourth gene header before values.
    matrix = pd.read_csv(DATA/'rnaseq_fpkm_20210329.csv', skiprows=[1,2,3], low_memory=False)
    gene_col, symbol_col = matrix.columns[:2]
    assert matrix[gene_col].notna().all() and not matrix[gene_col].duplicated().any()
    genes = pd.DataFrame({'gene_id':matrix[gene_col], 'symbol':matrix[symbol_col]})
    genes.to_csv(OUT/'gene_annotation.csv', index=False)
    expr = matrix.iloc[:,2:].apply(pd.to_numeric, errors='raise').T
    expr.columns = genes.gene_id.astype(str)
    expr.index.name = 'model_id'
    assert expr.index.is_unique and not np.isinf(expr.to_numpy()).any()
    assert np.nanmin(expr.to_numpy()) >= 0
    missing = expr.isna().mean(axis=1)
    audit = drug.merge(meta[keep_cols], left_on='SANGER_MODEL_ID', right_on='model_id', how='left', validate='one_to_one', indicator=True)
    audit['same_window'] = np.isclose(audit.MIN_CONC, 0.002001, atol=1e-9, rtol=0) & (audit.MAX_CONC == 2)
    audit['valid_response'] = np.isfinite(audit.AUC) & audit.AUC.between(0,1) & np.isfinite(audit.RMSE) & (audit.RMSE <= .3)
    audit['has_expression'] = audit.SANGER_MODEL_ID.isin(expr.index)
    audit['expression_missing_fraction'] = audit.SANGER_MODEL_ID.map(missing)
    audit['metadata_matched'] = audit['_merge'].eq('both')
    audit['eligible'] = audit.same_window & audit.valid_response & audit.has_expression & audit.metadata_matched & audit.model_type.eq('Cell Line') & audit.expression_missing_fraction.lt(.1)
    audit.to_csv(OUT/'cohort_audit.csv', index=False)
    cohort = audit[audit.eligible].sort_values('model_id').copy()
    assert len(cohort) >= 100
    # Union-find across all metadata so links through an excluded model are retained.
    parents = {x:x for x in meta.model_id}
    def find(x):
        parents.setdefault(x,x)
        if parents[x] != x: parents[x] = find(parents[x])
        return parents[x]
    def union(a,b):
        a,b = find(a), find(b)
        parents[max(a,b)] = min(a,b)
    for key in ['patient_id','sample_id']:
        for _, frame in meta.dropna(subset=[key]).groupby(key):
            ids = frame.model_id.tolist()
            for x in ids[1:]: union(ids[0],x)
    for row in meta.dropna(subset=['parent_id']).itertuples(): union(row.model_id,row.parent_id)
    cohort['group_id'] = cohort.model_id.map(find)
    cohort['tissue'] = cohort.tissue.fillna('Unknown')
    train, test = next(GroupShuffleSplit(n_splits=1,test_size=.2,random_state=42).split(cohort, groups=cohort.group_id))
    cohort['split'] = 'development'
    cohort.iloc[test,cohort.columns.get_loc('split')] = 'test'
    cohort['cv_fold'] = -1
    dev = cohort.iloc[train]
    for fold, (itr, iva) in enumerate(GroupKFold(n_splits=5,shuffle=True,random_state=42).split(dev,groups=dev.group_id)):
        assert not set(dev.iloc[itr].group_id) & set(dev.iloc[iva].group_id)
        cohort.loc[dev.iloc[iva].index,'cv_fold'] = fold
    assert not set(cohort[cohort.split=='test'].group_id) & set(dev.group_id)
    selected = ['model_id','CELL_LINE_NAME','group_id','tissue','cancer_type','AUC','LN_IC50','MIN_CONC','MAX_CONC','RMSE','split','cv_fold']
    cohort = cohort[selected].reset_index(drop=True)
    manifest_path = OUT/'split_manifest.csv'
    manifest_text = cohort.to_csv(index=False)
    if manifest_path.exists():
        assert manifest_path.read_text() == manifest_text, 'Frozen manifest changed; review source or software drift'
    else: manifest_path.write_text(manifest_text)
    X = np.log2(expr.loc[cohort.model_id].to_numpy(dtype=np.float32)+1)
    np.savez_compressed(DATA/'prepared_expression.npz', X=X, gene_ids=genes.gene_id.astype(str).to_numpy(dtype=str), model_ids=cohort.model_id.to_numpy(dtype=str))
    summary = {'drug':'Erlotinib', 'drug_id':1168, 'source_responses':len(drug),
        'uniform_window_responses':int(audit.same_window.sum()), 'eligible_models':len(cohort),
        'development_models':int((cohort.split=='development').sum()), 'test_models':int((cohort.split=='test').sum()),
        'groups':int(cohort.group_id.nunique()), 'genes':X.shape[1],
        'expression_missing_values':int(np.isnan(X).sum()), 'tissues':int(cohort.tissue.nunique()),
        'split_sha256':hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        'input_exclusions':{'other_concentration_window':int((~audit.same_window).sum()),
            'same_window_without_expression':int((audit.same_window & ~audit.has_expression).sum()),
            'same_window_invalid_response':int((audit.same_window & ~audit.valid_response).sum()),
            'same_window_metadata_unmatched':int((audit.same_window & ~audit.metadata_matched).sum())}}
    (OUT/'cohort_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))

if __name__ == '__main__': main()
