import pandas as pd
import os 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import math
from pathlib import Path 
from sklearn.model_selection import train_test_split, cross_validate, RepeatedKFold
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, HistGradientBoostingRegressor
from sklearn.svm import SVR
from sklearn.tree import DecisionTreeRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
from sklearn.neural_network import MLPRegressor
from sklearn import metrics
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
# from sklearn.evaluation import plot
from sklearn.inspection import permutation_importance
#|%%--%%| <cHsimX9I2S|IVyZtyz2hJ>
base = Path.cwd()
sv = base 

pfas = pd.read_csv('PFAS_data_130.csv')
pfas.head()

# %%
X = pfas[pfas.columns[10:]]
X

# %%
Y = pfas[pfas.columns[3]]
Y

# %%
XY = pd.concat([X,Y], axis = 1)
XY

# %% [markdown]
# ## 1. Feature selection with correlation matrix

# %%
XY_corr = XY.corr(method='pearson')
plt.figure(figsize=(15,8))
mask = np.triu(np.ones_like(XY_corr, dtype=bool))
sns.heatmap(XY_corr, mask=mask, annot=True, cmap='coolwarm', vmin=-1, vmax=1, fmt='.2f')
plt.title('Correlation Matrix')
plt.show()

#|%%--%%| <IVyZtyz2hJ|Tul2i4pvpT>
# %%
n_descps = len(X.columns) - 1
high_threshold = 0.9
low_threshold = 0.01

# excluding low correlation with a0
for descr in X.columns:
    res_y = stats.linregress(Y, X[descr])
    rsquared_y = res_y.rvalue**2
    if rsquared_y < low_threshold:
        print(f"{descr} has low correlation with a0, removing it")
        X.drop(columns=[descr], inplace = True)
        n_descps -= 1

# detecting descriptors that are strongly correlated with each other
descp_drop = set()
cols = list(X.columns)

for i, col1 in enumerate(cols):
    if col1 in descp_drop:
        continue  # skip if already removed

    for j in range(i + 1, len(cols)):
        col2 = cols[j]
        if col2 in descp_drop:
            continue  # skip if already removed

        res_x = stats.linregress(X[col2], X[col1])
        rsquared = res_x.rvalue**2

        if rsquared > high_threshold:
            print(f"{col1} has high correlation with {col2}")
            
            # Check which one is more strongly correlated with the target
            res_iy = stats.linregress(Y, X[col2])
            rsquared_iy = res_iy.rvalue**2
            res_jy = stats.linregress(Y, X[col1])
            rsquared_jy = res_jy.rvalue**2

            if rsquared_iy > rsquared_jy:
                print(f"--> Removing {col2}")
                descp_drop.add(col2)
            else:
                print(f"--> Removing {col1}")
                descp_drop.add(col1)
                break  # Stop checking col1 if it's removed

# drop descriptors
print("Dropping descriptors:", descp_drop)
X.drop(columns=descp_drop, inplace=True)
n_descps = len(X.columns) - 1
print(f"Number of descriptors after filtering: {n_descps}")
X.head()

# %% [markdown]
# ## 2. A glimpse of chemical space with UMAP and 4 different ways of scaling

# %%
from sklearn.preprocessing import StandardScaler, MinMaxScaler, MaxAbsScaler
import umap

# %%
# try 3 different ways of scaling data
X_SS = StandardScaler(with_mean=False).fit_transform(X)
X_MMS = MinMaxScaler().fit_transform(X)
X_MAS = MaxAbsScaler().fit_transform(X)

# %%
# UMAP for no scaling data
reducer = umap.UMAP(n_components = 2)
embedding = reducer.fit_transform(X)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1])
plt.gca().set_aspect('equal', 'datalim')

# %%
# UMAP for standard scaling data
reducer1 = umap.UMAP(n_components = 2)
embedding_SS = reducer1.fit_transform(X_SS)

plt.scatter(
    embedding_SS[:, 0],
    embedding_SS[:, 1])
plt.gca().set_aspect('equal', 'datalim')

# %%
# UMAP for min max scaling data
reducer2 = umap.UMAP(n_components = 2)
embedding_MMS = reducer2.fit_transform(X_MMS)

plt.scatter(
    embedding_MMS[:, 0],
    embedding_MMS[:, 1])
plt.gca().set_aspect('equal', 'datalim')

# %%
# UMAP for max abs scaling data
reducer3 = umap.UMAP(n_components = 2)
embedding_MAS = reducer3.fit_transform(X_MAS)

plt.scatter(
    embedding_MAS[:, 0],
    embedding_MAS[:, 1])
plt.gca().set_aspect('equal', 'datalim')

# %% [markdown]
# **We can see that data forms different clusters, so random split might not be a good idea, but I haven't figured out how to access to each cluster to split them yet.**

# %% [markdown]
# ## 3. Predicting a0 with Standard Scaling data, random split and some simple models

#|%%--%%| <Tul2i4pvpT|nzW9G2XKVa>



# %%
# Fit model and compute regression scores (R2), including CV:
def report_regression_scores(model, X_train, y_train, X_test, y_test, X, y, cv):
    # fit & predict
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_train_pred = model.predict(X_train)

    # --- holdout (train/test split) metrics ---
    train_r2   = metrics.r2_score(y_train, y_train_pred)
    test_r2    = metrics.r2_score(y_test, y_pred)

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
        "Train CV Scores": [train_cv_scores],
        "Test CV Scores": [test_cv_scores],
    }, index=[model.__class__.__name__])

    return summary_df


# %%
# Print permutation feature importance
def print_permutation_importance(model, X_train, y_train, X_test, y_test, feature_names, n_repeats=30, random_state=15):
    
    model.fit(X_train, y_train)
    r = permutation_importance(model, X_test, y_test, n_repeats=n_repeats, random_state=random_state)

    print(f"\nPermutation importances for {model.__class__.__name__}:")
    print(f"{'Feature':<20} Mean ± Std")
    for idx in r.importances_mean.argsort()[::-1]:
        mean, std = r.importances_mean[idx], r.importances_std[idx]
        if mean - 2*std > 0:
            print(f"{feature_names[idx]:<20} {mean:.3f} ± {std:.3f}")
    return r


# %%
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


# %%
# list out all the models I wanna try
kernel = ConstantKernel(1.0, (1e-3, 1e3)) * RBF(1.0, (1e-2, 1e2)) + WhiteKernel(noise_level=1e-3, noise_level_bounds=(1e-8, 1e1))

models = [Ridge(alpha=0.5, random_state=42), 
          Lasso(alpha=0.01, random_state=42),
          ElasticNet(alpha=0.01, l1_ratio=0.5, random_state=42),
          SVR(kernel='rbf', C=10, epsilon=0.1, gamma='scale'),
          KNeighborsRegressor(n_neighbors=11, weights="distance", metric="minkowski",p=1),
          GradientBoostingRegressor(n_estimators=300, learning_rate=0.05, max_depth=3, random_state=42),
          RandomForestRegressor(n_estimators=500, max_depth=None, min_samples_leaf=2, random_state=42),
          DecisionTreeRegressor(random_state=42),
          GaussianProcessRegressor(kernel=kernel, normalize_y=True, n_restarts_optimizer=5, alpha=1e-8, random_state=42),
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

# %% [markdown]
# ### Try different train/test split ratio

# %% [markdown]
# **Train/Test = 90:10**

# %%
X_train, X_test, y_train, y_test = train_test_split(X_SS, Y, test_size=0.1, random_state=42)
# Report metrics per model
rows = []
for m in models:
    df_row = report_regression_scores(m, X_train, y_train, X_test, y_test, X_SS, Y, cv=10)
    rows.append(df_row)
    
# concatenate all metrics into one table
all_summaries = pd.concat(rows)
display(all_summaries)

# %% [markdown]
# Clearly overfit: Test R2 < Train R2, and cross-validation (CV) R2 values vary by a lot.

# %% [markdown]
# **Train/Test = 83.3:16.7**

# %%
X_train, X_test, y_train, y_test = train_test_split(X_SS, Y, test_size=0.167, random_state=42)
# Report metrics per model
rows = []
for m in models:
    df_row = report_regression_scores(m, X_train, y_train, X_test, y_test, X_SS, Y, cv=6)
    rows.append(df_row)
    
# concatenate all metrics into one table
all_summaries = pd.concat(rows)
display(all_summaries)

# %% [markdown]
# Less overfit except GBR, RFR, DTR and kNN but the cross-validation (CV) values still vary a lot. --> Should be because of random split.

# %% [markdown]
# **Train/Test = 80/20**

# %%
X_train, X_test, y_train, y_test = train_test_split(X_SS, Y, test_size=0.2, random_state=42)
# Report metrics per model
rows = []
for m in models:
    df_row = report_regression_scores(m, X_train, y_train, X_test, y_test, X_SS, Y, cv=5)
    rows.append(df_row)
    
# concatenate all metrics into one table
all_summaries = pd.concat(rows)
display(all_summaries)

# big subplot of permutation feature importance
plot_permutation_importance(models, X_train, y_train, X_test, y_test, list(X.columns), cols=2, n_features=11)

# %% [markdown]
# * Much less overfitting except GBR, RFR, DTR, kNN but the cross-validation problem doesn't change.
# * Only 6-7 descriptors are meaningful for most of the models, except SVR --> Need more descriptors.

# %% [markdown]
# **Train/Test = 75:25**

# %%
X_train, X_test, y_train, y_test = train_test_split(X_SS, Y, test_size=0.25, random_state=42)
# Report metrics per model
rows = []
for m in models:
    df_row = report_regression_scores(m, X_train, y_train, X_test, y_test, X_SS, Y, cv=4)
    rows.append(df_row)
    
# concatenate all metrics into one table
all_summaries = pd.concat(rows)
display(all_summaries)

# big subplot of permutation feature importance
plot_permutation_importance(models, X_train, y_train, X_test, y_test, list(X.columns), cols=2, n_features=11)

# %% [markdown]
# * Similarly, the overfitting problem is much less except for GBR, RFR, DTR, and kNN.
# * Although the CV scores of the test set still vary a lot, it is less than all other train/test split ratios.
# * The only neural network model (MLPRegressor) performance decreases drastically, but we may fix that by tuning the hyperparameters of the model.
# * I personally prefer this ratio.
# * Similarly, only 6-7 descriptors matter, except SVR.

# %% [markdown]
# **Train/Test = 66.7:33.3**

# %%
# Try 66.7/33.3 train/test split 
X_train, X_test, y_train, y_test = train_test_split(X_SS, Y, test_size=0.333, random_state=42)
# Report metrics per model
rows = []
for m in models:
    df_row = report_regression_scores(m, X_train, y_train, X_test, y_test, X_SS, Y, cv=3)
    rows.append(df_row)
    
# concatenate all metrics into one table
all_summaries = pd.concat(rows)
display(all_summaries)

# %% [markdown]
# Slightly overfitting.

# %% [markdown]
# **Train/Test = 50:50**

# %%
X_train, X_test, y_train, y_test = train_test_split(X_SS, Y, test_size=0.5, random_state=42)
# Report metrics per model
rows = []
for m in models:
    df_row = report_regression_scores(m, X_train, y_train, X_test, y_test, X_SS, Y, cv=2)
    rows.append(df_row)
    
# concatenate all metrics into one table
all_summaries = pd.concat(rows)
display(all_summaries)

# %% [markdown]
# Clearly overfitting.

# %% [markdown]
# ### Overall conclusion:
# - GBR, RFR, and DTR are clearly overfitting across different train/test splitting ratios. They're also the models that are independent of scaling data. We can exclude them now.
# - 75:25 train/test split ratio gives the best results.
# - The maximum performance for now is 77-78%, which is decent because we haven't optimized everything yet.
# - Cross-validation scores are bad across all cases, showing that the problem lies somewhere else. That could be the random splitting. I'll try to find a way to access these clusters and split them.
# - Permutation feature importance shows that there are only 6-7 descriptors that matter. We might need more descriptors. I wonder if descriptors involving S atoms are missing here.

# %% [markdown]
# ## 4. Try fitting data with Max Abs Scaling

# %%
# list out all the models I wanna try
kernel = ConstantKernel(1.0, (1e-3, 1e3)) * RBF(1.0, (1e-2, 1e2)) + WhiteKernel(noise_level=1e-3, noise_level_bounds=(1e-8, 1e1))

models = [Ridge(alpha=0.5, random_state=42), 
          Lasso(alpha=0.01, random_state=42),
          ElasticNet(alpha=0.01, l1_ratio=0.5, random_state=42),
          SVR(kernel='rbf', C=10, epsilon=0.1, gamma='scale'),
          KNeighborsRegressor(n_neighbors=11, weights="distance", metric="minkowski",p=1),
          GaussianProcessRegressor(kernel=kernel, normalize_y=True, n_restarts_optimizer=5, alpha=1e-8, random_state=42),
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

# %%
X_train, X_test, y_train, y_test = train_test_split(X_MAS, Y, test_size=0.1, random_state=42)
# Report metrics per model
rows = []
for m in models:
    df_row = report_regression_scores(m, X_train, y_train, X_test, y_test, X_MAS, Y, cv=10)
    rows.append(df_row)
    
# concatenate all metrics into one table
all_summaries = pd.concat(rows)
display(all_summaries)

# %%
X_train, X_test, y_train, y_test = train_test_split(X_MAS, Y, test_size=0.167, random_state=42)
# Report metrics per model
rows = []
for m in models:
    df_row = report_regression_scores(m, X_train, y_train, X_test, y_test, X_MAS, Y, cv=6)
    rows.append(df_row)
    
# concatenate all metrics into one table
all_summaries = pd.concat(rows)
display(all_summaries)

# %%
X_train, X_test, y_train, y_test = train_test_split(X_MAS, Y, test_size=0.2, random_state=42)
# Report metrics per model
rows = []
for m in models:
    df_row = report_regression_scores(m, X_train, y_train, X_test, y_test, X_MAS, Y, cv=5)
    rows.append(df_row)
    
# concatenate all metrics into one table
all_summaries = pd.concat(rows)
display(all_summaries)

# %%
X_train, X_test, y_train, y_test = train_test_split(X_MAS, Y, test_size=0.25, random_state=42)
# Report metrics per model
rows = []
for m in models:
    df_row = report_regression_scores(m, X_train, y_train, X_test, y_test, X_MAS, Y, cv=4)
    rows.append(df_row)
    
# concatenate all metrics into one table
all_summaries = pd.concat(rows)
display(all_summaries)

# %%
X_train, X_test, y_train, y_test = train_test_split(X_MAS, Y, test_size=0.333, random_state=42)
# Report metrics per model
rows = []
for m in models:
    df_row = report_regression_scores(m, X_train, y_train, X_test, y_test, X_MAS, Y, cv=3)
    rows.append(df_row)
    
# concatenate all metrics into one table
all_summaries = pd.concat(rows)
display(all_summaries)

# %%
X_train, X_test, y_train, y_test = train_test_split(X_MAS, Y, test_size=0.5, random_state=42)
# Report metrics per model
rows = []
for m in models:
    df_row = report_regression_scores(m, X_train, y_train, X_test, y_test, X_MAS, Y, cv=2)
    rows.append(df_row)
    
# concatenate all metrics into one table
all_summaries = pd.concat(rows)
display(all_summaries)

# %% [markdown]
# ### Conclusion: Doesn't change much of the result:
# - Train/test split ratio: same pattern with 75/25 to be the best with the least overfitting
# - The cross-validation scores of the test sets are still a problem.
# - kNN clearly overfit --> Exclude.

# %% [markdown]
# ## 5. Total conclusion

# %% [markdown]
# ### Problems:
# - R2 might not be ideal yet. (But of course, we're at the screening step)
# - Overfitting: Cross-validation score differs so much!
# ### What can we improve?
# - The way we split data:
#     * Random split might not be ideal. --> Need help for UMAP visualization, then access each cluster to split.
# - The models we use:
#     * Fine-tuning the parameters of models. --> Need to learn.
#     * Use more models --> EasyML: need to learn.
#     * Neural network!!! --> need to learn.
# - We need to generate more descriptors.
# - We can update 200 more molecules (hopefully help with less overfitting).
# ### Fact check: Is it true that only neural networks can predict multiple Y values?

# %%


#|%%--%%| <nzW9G2XKVa|vmzEUsL53i>



