"""Reusable estimators and evaluation helpers; learned transforms are train-only."""
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.ensemble import RandomForestRegressor
from sklearn.dummy import DummyRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from scipy.stats import spearmanr

class ExpressionFeatures(TransformerMixin, BaseEstimator):
    def __init__(self, n_features=1000): self.n_features = n_features
    def fit(self, X, y=None):
        X = np.asarray(X, dtype=np.float64)
        valid = np.any(np.isfinite(X), axis=0)
        indices = np.flatnonzero(valid)
        variance = np.nanvar(X[:,indices], axis=0)
        order = np.argsort(-variance, kind='stable')
        self.indices_ = indices[order[variance[order] > 0]][:self.n_features]
        assert len(self.indices_) > 0
        self.medians_ = np.nanmedian(X[:,self.indices_], axis=0)
        return self
    def transform(self, X):
        result = np.asarray(X, dtype=np.float64)[:,self.indices_].copy()
        return np.where(np.isnan(result), self.medians_, result)

def make_model(name, params=None):
    params = params or {}
    if name == 'Mean baseline': return DummyRegressor(strategy='mean')
    if name == 'Tissue-only ridge':
        return Pipeline([('encode',OneHotEncoder(handle_unknown='ignore',sparse_output=False)),('model',Ridge(**params))])
    estimator = (ElasticNet(max_iter=20000,tol=.0001,random_state=42,**params) if name == 'Elastic Net'
                 else RandomForestRegressor(n_estimators=300,random_state=42,n_jobs=2,**params))
    steps = [('features',ExpressionFeatures())]
    if name == 'Elastic Net': steps.append(('scale',StandardScaler()))
    return Pipeline(steps+[('model',estimator)])

def predict(model, X): return np.clip(model.predict(X),0,1)

def metrics(y,p):
    rho = float(spearmanr(y,p).statistic) if np.std(p)>1e-12 and np.std(y)>1e-12 else None
    return {'RMSE':float(np.sqrt(mean_squared_error(y,p))), 'MAE':float(mean_absolute_error(y,p)),
            'R2':float(r2_score(y,p)), 'Spearman':rho}
