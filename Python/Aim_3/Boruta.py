import sys
import os.path
import pandas as pd
import numpy as np
import time
import pyreadr
from boruta import BorutaPy
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, RepeatedKFold, GroupShuffleSplit
from sklearn.ensemble import RandomForestClassifier

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
n_tune_control = int(n_control * 0.2)
n_train_control = n_control - n_test_control - n_tune_control

n_test_disease = int(n_disease * 0.2)
n_tune_disease = int(n_disease * 0.2)
n_train_disease = n_disease - n_test_disease - n_tune_disease

# Randomly assign controls to each dataset
test_control_individuals = np.random.choice(
    df[df['class'] == 0]['individual'].unique(),
    size=n_test_control,
    replace=False
)
tune_control_individuals = np.random.choice(
    np.setdiff1d(
        df[df['class'] == 0]['individual'].unique(),
        test_control_individuals
    ),
    size=n_tune_control,
    replace=False
)
train_control_individuals = np.setdiff1d(
    df[df['class'] == 0]['individual'].unique(),
    np.concatenate([test_control_individuals, tune_control_individuals])
)

# Randomly assign disease samples to each dataset
test_disease_individuals = np.random.choice(
    df[df['class'] == 1]['individual'].unique(),
    size=n_test_disease,
    replace=False
)
tune_disease_individuals = np.random.choice(
    np.setdiff1d(
        df[df['class'] == 1]['individual'].unique(),
        test_disease_individuals
    ),
    size=n_tune_disease,
    replace=False
)
train_disease_individuals = np.setdiff1d(
    df[df['class'] == 1]['individual'].unique(),
    np.concatenate([test_disease_individuals, tune_disease_individuals])
)

# Get the corresponding cells for each dataset
test_index = df['individual'].isin(np.concatenate([test_control_individuals, test_disease_individuals]))
tune_index = df['individual'].isin(np.concatenate([tune_control_individuals, tune_disease_individuals]))
train_index = df['individual'].isin(np.concatenate([train_control_individuals, train_disease_individuals]))

# Split data into training, tuning, and testing sets
X_train, X_test, X_tune = df.loc[train_index,].drop(['class', 'individual'], axis=1), df.loc[test_index,].drop(['class', 'individual'], axis=1), df.loc[tune_index,].drop(['class', 'individual'], axis=1)
y_train, y_test, y_tune = df.loc[train_index, 'class'], df.loc[test_index, 'class'], df.loc[tune_index, 'class']

# Standard scale the data - z-scores
scaler = StandardScaler()
X_train = pd.DataFrame(scaler.fit_transform(X_train), columns=X_train.columns, index=X_train.index)
X_tune = pd.DataFrame(scaler.fit_transform(X_tune), columns=X_tune.columns, index=X_tune.index)
X_test = pd.DataFrame(scaler.fit_transform(X_test), columns=X_test.columns, index=X_test.index)

# Save the data to temporary files
X_train.to_csv('data.splits/X_train.'+cell+'.csv', index=False)
y_train.to_csv('data.splits/y_train.'+cell+'.csv', index=False)
X_tune.to_csv('data.splits/X_tune.'+cell+'.csv', index=False)
y_tune.to_csv('data.splits/y_tune.'+cell+'.csv', index=False)
X_test.to_csv('data.splits/X_test.'+cell+'.csv', index=False)
y_test.to_csv('data.splits/y_test.'+cell+'.csv', index=False)

### Boruta feature selection ###
X = X_tune.values
y = y_tune.ravel()
# random forest classifier utilising all cores and sampling in proportion to y labels
param_grid = {'n_estimators': [100, 200, 300, 400],
              'max_features': ['sqrt', 'log2', 0.3],
                'max_depth': [3, 5, 7],
                'min_samples_split': [2, 5, 8, 10],
                'ccp_alpha': [0, 0.01, 0.1, 0.1]
}
clf = RandomForestClassifier(n_jobs=-1)
grid_search = GridSearchCV(clf, param_grid, cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0), n_jobs=-1, verbose=1)
# Fit the grid search object to the training data
grid_search.fit(X, y)
# Create an RFECV object with a random forest classifier
rf = RandomForestClassifier(n_estimators=grid_search.best_params_['n_estimators'], 
                            max_features=grid_search.best_params_['max_features'],
                            max_depth=grid_search.best_params_['max_depth'], 
                            min_samples_split=grid_search.best_params_['min_samples_split'], 
                            ccp_alpha=grid_search.best_params_['ccp_alpha'],
                            n_jobs=-1)
# define Boruta feature selection method
feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, random_state=1)
# find all relevant features - 5 features should be selected
feat_selector.fit(X, y)
# Return features
features = X_tune.columns[feat_selector.support_].tolist()

# Save the features to a temporary file
features = pd.DataFrame(features)
features.to_csv(sys.args[2]+'/features.'+sys.args[3]+'.csv', index=False)

cv = StratifiedKFold(5)
rfecv = RFECV(
    estimator=rf,
    step=1,
    cv=cv,
    scoring="accuracy",
    min_features_to_select=1,
    n_jobs=-1,
)
rfecv.fit(X_tune, y_tune)

features = rfecv.get_feature_names_out()

