import sys
import os.path
import time
import pandas as pd
import numpy as np
from numpy import arange
import time
import pyreadr
from boruta import BorutaPy
from sklearn.utils import resample
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, RepeatedKFold, GroupShuffleSplit
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import ElasticNetCV
import pickle

start_time = time.time()

# Get the file name from the command line
file = sys.argv[1]
print(os.path.basename(file))
cell = os.path.basename(file).replace('.RDS', '')

# Read in expression RDS file
df = pyreadr.read_r(file)
df = df[None]
print(df.head())

# Remove 'Hispanic or Latin American' ancestry
df = df[df['ancestry'] != 'Hispanic or Latin American']

# Replace class with binary label
df['class'] = df['class'].replace({"control": 0, "disease": 1})
# Replace ancestry with binary label
df['ancestry'] = df['ancestry'].replace({"European": 0, "Asian": 1})

# Save ancestry to add as a feature later
ancestry = df['ancestry']

# # Get the individual IDs for the training and testing sets from the old analysis
# train_ids = pd.read_csv('old_psuedobulk/data.splits/X_train.' + os.path.basename(file).replace('.RDS', '.csv'))
# train_ids = train_ids['rownames']
# test_ids = pd.read_csv('old_psuedobulk/data.splits/X_test.' + os.path.basename(file).replace('.RDS', '.csv'))
# test_ids = test_ids['rownames']

# # Get the training and testing data
# X_train = df[df['individual'].isin(train_ids)].drop(['class', 'individual', 'ancestry'], axis=1)
# X_test = df[df['individual'].isin(test_ids)].drop(['class', 'individual', 'ancestry'], axis=1)
# y_train = df[df['individual'].isin(train_ids)]['class']
# y_test = df[df['individual'].isin(test_ids)]['class']

# Read in data splits
train = pd.read_csv(f'new_pseudobulk/split_{sys.argv[2]}/train_index.csv')
test = pd.read_csv(f'new_pseudobulk/split_{sys.argv[2]}/test_index.csv')

# Get the training and testing data
X_train = df[df['individual'].isin(train['rownames'])].drop(['class', 'individual', 'ancestry'], axis=1)
y_train = df[df['individual'].isin(train['rownames'])]['class']
X_test = df[df['individual'].isin(test['rownames'])].drop(['class', 'individual', 'ancestry'], axis=1)
y_test = df[df['individual'].isin(test['rownames'])]['class']

# Standard scale the data - z-scores
scaler = StandardScaler()
X_train = pd.DataFrame(scaler.fit_transform(X_train), columns=X_train.columns, index=X_train.index)
X_test = pd.DataFrame(scaler.fit_transform(X_test), columns=X_test.columns, index=X_test.index)

# Add ancestry as a feature
X_train['ancestry'] = ancestry[X_train.index]

# Save data splits
X_train.to_csv(f'new_pseudobulk/split_{sys.argv[2]}/data.splits/X_train.'+cell+'.csv', index=True)
y_train.to_csv(f'new_pseudobulk/split_{sys.argv[2]}/data.splits/y_train.'+cell+'.csv', index=True)
X_test.to_csv(f'new_pseudobulk/split_{sys.argv[2]}/data.splits/X_test.'+cell+'.csv', index=True)
y_test.to_csv(f'new_pseudobulk/split_{sys.argv[2]}/data.splits/y_test.'+cell+'.csv', index=True)

### Boruta feature selection ###
X = X_train.values
y = y_train.ravel()
# random forest classifier utilising all cores and sampling in proportion to y labels
param_grid = {'n_estimators': [100, 200, 300, 400],
              'criterion': ['gini', 'entropy'],
              'max_features': ['sqrt', 'log2', 0.3],
                'max_depth': [3, 4, 5, 6, 7],
                'min_samples_split': [2, 5, 8, 10]
}
clf = RandomForestClassifier(n_jobs=-1, class_weight='balanced', random_state=42)
grid_search = GridSearchCV(clf, param_grid, 
                           cv=RepeatedKFold(n_splits=10, n_repeats=3, 
                                            random_state=42), 
                                            n_jobs=8, 
                                            verbose=1,
                                            scoring='accuracy')
# Fit the grid search object to the training data
grid_search.fit(X, y)
# Create an RFECV object with a random forest classifier
rf = RandomForestClassifier(n_estimators=grid_search.best_params_['n_estimators'], 
                            max_features=grid_search.best_params_['max_features'],
                            max_depth=grid_search.best_params_['max_depth'], 
                            min_samples_split=grid_search.best_params_['min_samples_split'],
                            class_weight='balanced',
                            n_jobs=-1)
# define Boruta feature selection method
feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, random_state=1)
# find all relevant features
feat_selector.fit(X, y)

# Save the model
filename = f'new_pseudobulk/split_{sys.argv[2]}/feature.select.model/boruta_'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(feat_selector, open(filename, 'wb'))

# Get feature rankings
feature_ranks = list(feat_selector.ranking_)
# Create a DataFrame with features, their importance, and ranks
feature_df = pd.DataFrame({
    'Feature': X_train.columns,
    'Rank': feature_ranks
})
# Sort the DataFrame based on feature ranks
feature_df.sort_values(by='Rank', ascending=True, inplace=True)
# Save the feature importance to file
feature_df.to_csv(f'new_pseudobulk/split_{sys.argv[2]}/features/boruta_features.'+cell+'.csv', index=False)

### Elastic net feature selection ###
ratios = arange(0, 1.1, 0.1)
alphas = np.logspace(-4, 0, 10)
cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=42)
enet = ElasticNetCV(l1_ratio=ratios, alphas=alphas, cv=cv, n_jobs=8, random_state=42)
enet.fit(X_train, y_train.ravel())
print(enet)

# Create a dataframe of the features and their coefficients
enet_features = pd.DataFrame({'Feature': enet.feature_names_in_, 'coef': enet.coef_})
# Save the features to file
enet_features.to_csv(f'new_pseudobulk/split_{sys.argv[2]}/features/enet_features.'+cell+'.csv', index=False)

# Save the model
filename = f'new_pseudobulk/split_{sys.argv[2]}/feature.select.model/enet_'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(enet, open(filename, 'wb'))

end_time = time.time()
cpu_time = end_time - start_time

print(f"CPU time used: {cpu_time:.2f} seconds")