#|%%--%%| <g89zrDTtfA|MMBxSfJ3Iu>
import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
#|%%--%%| <MMBxSfJ3Iu|Lt5iRPfRLu>

# In[2]:


pfas = pd.read_csv('no_sn_202.csv')
pfas


# ## 1. Data preparation before fitting models

# In[3]:


X = pfas[pfas.columns[2:49]]
X


# In[4]:


Y = pfas[pfas.columns[49:54]]
Y


# In[5]:


Y1 = pfas[pfas.columns[50:54]]
Y1


# In[6]:


XY = pd.concat([X,Y], axis = 1)
XY


# In[7]:


import matplotlib.pyplot as plt
import seaborn as sns
XY_corr = XY.corr(method='pearson')
plt.figure(figsize=(15,8))
mask = np.triu(np.ones_like(XY_corr, dtype=bool))
sns.heatmap(XY_corr, mask=mask, annot=True, cmap='coolwarm', vmin=-1, vmax=1, fmt='.1f')
plt.title('Correlation Matrix')
plt.show()


# ## 2. Fitting models

# In[8]:


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


# In[9]:


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


# In[10]:


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


# In[11]:


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


# In[12]:


# Try with data scaled with MinMaxScaler
from sklearn.preprocessing import MinMaxScaler, StandardScaler
X_SS = StandardScaler().fit_transform(X)
X_MMS = MinMaxScaler().fit_transform(X)
X_train, X_test, y_train, y_test = train_test_split(X_MMS, Y, test_size=0.25, random_state=42)
rows = []
for m in models:
    df_row = report_regression_scores(m, X_train, y_train, X_test, y_test, X_MMS, Y, 4)
    rows.append(df_row)

# concatenate all metrics into one table
all_summaries = pd.concat(rows)
display(all_summaries)

#|%%--%%| <Lt5iRPfRLu|pCCHUoGVso>



