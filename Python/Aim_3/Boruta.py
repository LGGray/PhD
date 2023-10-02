import sys
import os.path
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

# Get the file name from the command line
file = sys.argv[1]
print(os.path.basename(file))
cell = os.path.basename(file).replace('.RDS', '')

# Read in expression RDS file
df = pyreadr.read_r(file)
df = df[None]
print(df.head())

# Replace classes with binary label
if sum(df['class'] == 'control') > 0:
  df['class'] = df['class'].replace({"control": 0, "disease": 1})
else:
  df['class'] = df['class'].replace({"managed": 0, "flare": 1})

# # Downsample the majority class - replace=False to ensure no duplicates
# df_majority = df[df['class'] == 1]
# df_minority = df[df['class'] == 0]
# df_majority_downsampled = resample(df_majority,
#                                    replace=False,
#                                    n_samples=len(df_minority),
#                                    random_state=0)
# df = pd.concat([df_majority_downsampled, df_minority])

# # Write downsamples to .RDS file
# pyreadr.write_rds('exp.matrix/downsampled.'+cell+'.RDS', df)

### Split the data into train, tune and test sets ###

# Collect individual IDs
individuals = df['individual'].unique()
n_individuals = len(individuals)

# Get the number of individuals in each condition
individual_class = df['individual'].astype(str) + '_' + df['class'].astype(str)
n_control = len(individual_class[individual_class.str.endswith('_0')].unique())
n_disease = len(individual_class[individual_class.str.endswith('_1')].unique())

# Determine number of controls and disease samples to include in each dataset
n_test_control = int(n_control * 0.2)
n_train_control = n_control - n_test_control

n_test_disease = int(n_disease * 0.2)
n_train_disease = n_disease - n_test_disease

# Randomly assign controls to each dataset
test_control_individuals = np.random.choice(
    df[df['class'] == 0]['individual'].unique(),
    size=n_test_control,
    replace=False
)
train_control_individuals = np.setdiff1d(
    df[df['class'] == 0]['individual'].unique(),
    np.concatenate([test_control_individuals])
)

# Randomly assign disease samples to each dataset
test_disease_individuals = np.random.choice(
    df[df['class'] == 1]['individual'].unique(),
    size=n_test_disease,
    replace=False
)
train_disease_individuals = np.setdiff1d(
    df[df['class'] == 1]['individual'].unique(),
    np.concatenate([test_disease_individuals])
)

# Get the corresponding cells for each dataset
test_index = df['individual'].isin(np.concatenate([test_control_individuals, test_disease_individuals]))
train_index = df['individual'].isin(np.concatenate([train_control_individuals, train_disease_individuals]))

# Split data into training, tuning, and testing sets
X_train, X_test = df.loc[train_index,].drop(['class', 'individual'], axis=1), df.loc[test_index,].drop(['class', 'individual'], axis=1)
y_train, y_test = df.loc[train_index, 'class'], df.loc[test_index, 'class']

# Standard scale the data - z-scores
scaler = StandardScaler()
X_train = pd.DataFrame(scaler.fit_transform(X_train), columns=X_train.columns, index=X_train.index)
X_test = pd.DataFrame(scaler.fit_transform(X_test), columns=X_test.columns, index=X_test.index)

# Save the data to temporary files
X_train.to_csv('psuedobulk/X_train.'+cell+'.csv', index=True)
y_train.to_csv('psuedobulk/y_train.'+cell+'.csv', index=True)
X_test.to_csv('psuedobulk/X_test.'+cell+'.csv', index=True)
y_test.to_csv('psuedobulk/y_test.'+cell+'.csv', index=True)

### Boruta feature selection ###
X = X_train.values
y = y_train.ravel()
# random forest classifier utilising all cores and sampling in proportion to y labels
param_grid = {'n_estimators': [100, 200, 300, 400],
              'criterion': ['gini', 'entropy'],
              'max_features': ['sqrt', 'log2', 0.3],
                'max_depth': [3, 5, 7],
                'min_samples_split': [2, 5, 8, 10]
}
clf = RandomForestClassifier(n_jobs=-1, class_weight='balanced')
grid_search = GridSearchCV(clf, param_grid, cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0), n_jobs=-1, verbose=1)
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
# find all relevant features - 5 features should be selected
feat_selector.fit(X, y)
# Return features
boruta_features = X_train.columns[feat_selector.support_].tolist()
# Save the features to file
pd.DataFrame(boruta_features).to_csv('psuedobulk/features/boruta_features.'+cell+'.csv', index=False)

### Elastic net feature selection ###
ratios = arange(0, 1, 0.1)
alphas = np.logspace(-4, 0, 10)
cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0)
enet = ElasticNetCV(l1_ratio=ratios, alphas=alphas, cv=cv, n_jobs=-1, random_state=0)
enet.fit(X_train, y_train.ravel())

enet_features = X_train.columns[enet.coef_ > 0].tolist()
# Save the features to file
pd.DataFrame(enet_features).to_csv('psuedobulk/features/enet_features.'+cell+'.csv', index=False)

# Save the model
import pickle
filename = 'psuedobulk/feature.select.model/enet_'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(enet, open(filename, 'wb'))

filename = 'psuedobulk/feature.select.model/boruta_'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(feat_selector, open(filename, 'wb'))

