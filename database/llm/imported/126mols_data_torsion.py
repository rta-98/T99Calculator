#!/usr/bin/env python
# coding: utf-8

# ## 0. Load data

# In[1]:


import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


# In[2]:


pfas = pd.read_csv('no_sn_torsions_202.csv')
pfas


# ## 1. Data preparation before fitting models

# In[41]:


X = pfas[pfas.columns[2:34]]
X


# In[4]:


Y = pfas[pfas.columns[34:39]]
Y


# In[6]:


Y1 = pfas[pfas.columns[35:39]]
Y1


# In[7]:


XY = pd.concat([X,Y], axis = 1)
XY


# In[19]:


import matplotlib.pyplot as plt
import seaborn as sns
XY_corr = XY.corr(method='pearson')
plt.figure(figsize=(15,8))
mask = np.triu(np.ones_like(XY_corr, dtype=bool))
sns.heatmap(XY_corr, mask=mask, annot=True, cmap='coolwarm', vmin=-1, vmax=1, fmt='.1f')
plt.title('Correlation Matrix')
plt.show()


# In[39]:


XY_corr[XY_corr.columns[32:37]]


# High correlation with each other (>=0.9)
# - C-C-C-F and C-C-C-C
# - F-C-C-C and C-C-C-C
# - F-C-C-F and C-C-C-F
# - F-C-C-F and F-C-C-C
# - Why are F-C-C-C and C-C-C-F torsions 2 different things and only correlate with each other 0.841904?
# --> (C-C-C-C and/or F-C-C-F) correlate to (F-C-C-C and C-C-C-F) --> Drop C-C-C-C and F-C-C-F.
# - C-C-C-H and F-C-C-H
# - C-C-C-H and H-C-C-H
# --> C-C-C-H correlates to (F-C-C-H and H-C-C-H), but F-C-C-H and H-C-C-H don't strongly correlates to each other.
#   --> Drop C-C-C-H.
# - F-C-O-C and C-C-O-C
# - C-O-C-C and F-C-O-C
# --> Drop F-C-O-C.
# Low correlation with all targets (<= 0.1): 14 descriptors
# H-O-C-C, O-C-C-C, O-C-C-F, F-C-C-O, C-C-O-H, C-C-C-H, C-C-O-C, C-O-C-H, H-C-C-F, C-O-C-O, O-C-O-C

# In[42]:


X.drop(['Torsions H-O-C-C', 'Torsions O-C-C-C', 'Torsions O-C-C-F', 'Torsions F-C-C-O', 'Torsions C-C-O-H', 'Torsions C-C-O-C', 'Torsions C-O-C-H', 'Torsions H-C-C-F', 'Torsions C-O-C-O', 'Torsions O-C-O-C'], axis=1, inplace=True)
X


# In[44]:


X.drop(['Torsions C-C-C-C', 'Torsions F-C-C-F', 'Torsions C-C-C-H'], axis=1, inplace=True)
X


# ## 2. Fitting models

# In[45]:


from sklearn.model_selection import train_test_split, cross_validate, RepeatedKFold
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet, MultiTaskElasticNet, MultiTaskLasso
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, HistGradientBoostingRegressor, ExtraTreesRegressor
from sklearn.svm import SVR
from sklearn.tree import DecisionTreeRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
from sklearn.neural_network import MLPRegressor
from sklearn import metrics
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, make_scorer
from sklearn_evaluation import plot
from sklearn.inspection import permutation_importance


# In[46]:


# Fit model and compute regression scores (R2), including CV:
def report_regression_scores(model, X_train, y_train, X_test, y_test, X, y, cv):
    # fit & predict
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_train_pred = model.predict(X_train)

    # --- holdout (train/test split) metrics ---
    train_r2   = metrics.r2_score(y_train, y_train_pred)
    test_r2    = metrics.r2_score(y_test, y_pred)
    train_MAE = metrics.mean_absolute_error(y_train, y_train_pred)
    test_MAE = metrics.mean_absolute_error(y_test, y_pred)

    # --- cross-validation metrics ---
    scores = cross_validate(
        model, X, y,
        cv=cv,                    
        return_train_score=True,
        scoring='r2',
    )
    scores_table = pd.DataFrame(scores)

    train_cv_scores = np.round(scores_table["train_score"].values, 3)
    test_cv_scores  = np.round(scores_table["test_score"].values, 3)

    # Build & return a summary table
    summary_df = pd.DataFrame({
        "Train R2":   [round(train_r2,3)],
        "Test R2":    [round(test_r2,3)],
        "Train MAE":    [round(train_MAE,3)],
        "Test MAE":    [round(test_MAE,3)],
        "Train CV Scores": [train_cv_scores],
        "Test CV Scores": [test_cv_scores],
    }, index=[model.__class__.__name__])

    return summary_df


# In[47]:


import math
# Plot permutation feature importance
def plot_permutation_importance(models, X_train, y_train, X_test, y_test, feature_names, cols, n_features, n_repeats=30, random_state=15):
    # create subplots
    n = len(models)
    rows = math.ceil(n/cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols*5, rows*4))
    if n == 1:
        axes = [axes]    
    axes = axes.flatten()

    # fit model and plot permutation feature importance 
    for ax, model in zip(axes, models):
        model.fit(X_train, y_train)
        r = permutation_importance(model, X_test, y_test, n_repeats=n_repeats, random_state=random_state)
        idxs = r.importances_mean.argsort()[::-1]

        ax.barh(range(len(idxs)), r.importances_mean[idxs], xerr=r.importances_std[idxs], align='center')
        ax.set_yticks(range(len(idxs)))
        ax.set_yticklabels([feature_names[i] for i in idxs])
        ax.invert_yaxis()
        ax.set_title(f"{model.__class__.__name__}")
        ax.set_xlabel("Permutation Importance")

    for empty_ax in axes[len(models):]:
        empty_ax.set_visible(False)

    plt.tight_layout()
    plt.savefig(f"PFI_{n_features}features.png",        
            dpi=300,
            bbox_inches='tight')
    plt.show()


# In[48]:


# list out all the models I wanna try
models = [LinearRegression(),
          Ridge(alpha=0.5, random_state=42),
          Lasso(alpha=0.01, max_iter=10000, random_state=42),
          MultiTaskLasso(alpha=0.01, max_iter=10000, random_state=42),
          ElasticNet(alpha=0.01, l1_ratio=0.5, max_iter=10000, random_state=42),
          MultiTaskElasticNet(alpha=0.01, l1_ratio=0.5, max_iter=10000, random_state=42),
          RandomForestRegressor(n_estimators=500, max_depth=None, min_samples_leaf=2, random_state=42, n_jobs=-1),
          ExtraTreesRegressor(n_estimators=500, random_state=42, n_jobs=-1),
          DecisionTreeRegressor(random_state=42),
          KNeighborsRegressor(n_neighbors=3, weights="distance", metric="minkowski", p=1),
          MLPRegressor(
              hidden_layer_sizes=(32,16),
              activation="relu",
              solver="adam",
              alpha=1e-3,
              learning_rate_init=1e-3,
              max_iter=5000,
              early_stopping=True,
              validation_fraction=0.2,
              n_iter_no_change=30,
              random_state=1)]


# In[49]:


from sklearn.preprocessing import MinMaxScaler, StandardScaler
X_MMS = MinMaxScaler().fit_transform(X)
Y_SS = StandardScaler().fit_transform(Y)
X_train, X_test, y_train, y_test = train_test_split(X_MMS, Y_SS, test_size=0.25, random_state=42)
rows = []
for m in models:
    df_row = report_regression_scores(m, X_train, y_train, X_test, y_test, X_MMS, Y_SS, 4)
    rows.append(df_row)

# concatenate all metrics into one table
all_summaries = pd.concat(rows)
all_summaries.round(3)
display(all_summaries)


# In[ ]:





# In[50]:


from sklearn.preprocessing import MinMaxScaler, StandardScaler
X_MMS = MinMaxScaler().fit_transform(X)
Y1_SS = StandardScaler().fit_transform(Y1)
X_train, X_test, y_train, y_test = train_test_split(X_MMS, Y1_SS, test_size=0.25, random_state=42)
rows = []
for m in models:
    df_row = report_regression_scores(m, X_train, y_train, X_test, y_test, X_MMS, Y1_SS, 4)
    rows.append(df_row)

# concatenate all metrics into one table
all_summaries = pd.concat(rows)
all_summaries.round(3)
display(all_summaries)


# In[ ]:




