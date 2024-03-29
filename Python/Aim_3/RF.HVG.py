# Import required libraries
import sys
import glob
import os.path
import time
import pandas as pd
import numpy as np
import pyreadr
from boruta import BorutaPy
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RepeatedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score
from sklearn.feature_selection import RFECV
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.inspection import permutation_importance
import matplotlib.pyplot as plt
import pickle

start_time = time.process_time()

# Get the file name from the command line
file = sys.argv[1]
print(os.path.basename(file))
test = sys.argv[2]
print(test)

# Read in expression RDS file
df = pyreadr.read_r(file)
df = df[None]
print(df.head())

cell = file.replace('exp.matrix/', '').replace('.RDS', '')

# Replace classes with binary label
if sum(df['class'] == 'control') > 0:
  df['class'] = df['class'].replace({"control": 0, "disease": 1})
else:
  df['class'] = df['class'].replace({"managed": 0, "flare": 1})

# Collect individual IDs
individuals = df['individual'].unique()
n_individuals = len(individuals)

# Get the number of individuals in each condition
individual_class = df['individual'].astype(str) + '_' + df['class'].astype(str)
n_control = len(individual_class[individual_class.str.endswith('_0')].unique())
n_disease = len(individual_class[individual_class.str.endswith('_1')].unique())
    
# Calculate the number of individuals to assign to each dataset
n_test = int(n_individuals * 0.2)
n_tune = int(n_individuals * 0.2)
n_train = int(n_individuals * 0.6)

# Randomly assign individuals to each dataset
np.random.seed(42)
test_individuals = np.random.choice(individuals, size=n_test, replace=False)
individuals = np.setdiff1d(individuals, test_individuals)
tune_individuals = np.random.choice(individuals, size=n_tune, replace=False)
individuals = np.setdiff1d(individuals, tune_individuals)
train_individuals = individuals

# Get the corresponding cells for each dataset
test_index = df['individual'].isin(test_individuals)
tune_index = df['individual'].isin(tune_individuals)
train_index = df['individual'].isin(train_individuals)

# Split the data into training, tuning, and testing sets
X_train, X_test, X_tune = df.loc[train_index,].drop(['class','individual'], axis=1), df.loc[test_index,].drop(['class','individual'], axis=1), df.loc[tune_index,].drop(['class','individual'], axis=1)
y_train, y_test, y_tune = df.loc[train_index,'class'], df.loc[test_index,'class'], df.loc[tune_index,'class']

# # Boruta feature selection
# X = X_tune.values
# y = y_tune.ravel()
# # random forest classifier utilising all cores and sampling in proportion to y labels
# rf = RandomForestClassifier(n_jobs=-1, class_weight='balanced', max_depth=5)
# # define Boruta feature selection method
# feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, random_state=1)
# # find all relevant features - 5 features should be selected
# feat_selector.fit(X, y)
# # Return features
# features = X_tune.columns[feat_selector.support_].tolist()

# load the model from disk
RF = pickle.load(open('ML.models/RF_model_'+cell+'.sav', 'rb'))
# Read in chrX genes
chrX = pd.read_csv('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.csv', header=None)[1].values
features = RF.feature_names_in_

# Add condition to remove chrX genes or randomly sampled non-chrX genes
if test == '-X':
    chrX_features = features[np.isin(features, chrX)]
    features = [f for f in features if f not in chrX_features]
elif test == '-random':
    chrX_features = features[np.isin(features, chrX)]
    non_chrX_features = features[~np.isin(features, chrX)]
    non_chrX_sample = np.random.choice(non_chrX_features, size=len(chrX_features), replace=False)
    features = features[~np.isin(features, non_chrX_sample)]

# Perform a grid search to find the best parameters
# Create the parameter grid
param_grid = {'n_estimators': [100, 200, 300, 400],
              'max_features': ['sqrt', 'log2', 0.3],
                'max_depth': [5, 10, 15, 30],
                'min_samples_split': [2, 5, 8, 10]
}
clf = RandomForestClassifier(n_jobs=-1)
grid_search = GridSearchCV(clf, param_grid, cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0), n_jobs=-1, verbose=1)

# Fit the grid search object to the training data
grid_search.fit(X_tune.loc[:,features], y_tune)

# Create an RFECV object with a random forest classifier
clf = RandomForestClassifier(n_estimators=grid_search.best_params_['n_estimators'], 
                            max_depth=grid_search.best_params_['max_depth'], 
                            min_samples_split=grid_search.best_params_['min_samples_split'], n_jobs=-1)
# Fit the model
clf.fit(X_train.loc[:, features], y_train)
# Predict the test set
y_pred = clf.predict(X_test.loc[:, features])
# Get the predicted probabilities
y_pred_proba = clf.predict_proba(X_test.loc[:, features])[:, 1]

# Calculate the metrics
accuracy = accuracy_score(y_test, y_pred)
precision = precision_score(y_test, y_pred)
recall = recall_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred)
auc = roc_auc_score(y_test, y_pred)

# Create dataframe of metrics and save to file
metrics = pd.DataFrame({'Accuracy': [accuracy], 
                        'Precision': [precision], 
                        'Recall': [recall], 
                        'F1': [f1], 
                        'AUC': [auc]})
metrics.to_csv('exp.matrix/metrics/RF_metrics_'+os.path.basename(file).replace('.RDS', test)+'.csv', index=False)

# Save confusion matrix to file
confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
confusion.to_csv('exp.matrix/metrics/RF_confusion_'+os.path.basename(file).replace('.RDS', test)+'.csv', index=False)

# Print the AUROC curve
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
plt.plot(fpr, tpr, label='AUC-ROC (area = %0.2f)' % auc)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('RF: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.legend(loc="lower right")
plt.savefig('exp.matrix/AUROC/RF_'+os.path.basename(file).replace('.RDS', test)+'.pdf', bbox_inches='tight')

# Print the PR curve
precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
average_precision = average_precision_score(y_test, y_pred_proba)
disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
disp.plot()
disp.ax_.set_title('RF: ' + os.path.basename(file).replace('.RDS', test).replace('.', ' '))
plt.savefig('exp.matrix/PRC/RF_'+os.path.basename(file).replace('.RDS', test)+'.pdf', bbox_inches='tight')

# Save the model
import pickle
filename = 'ML.models/RF_model_'+os.path.basename(file).replace('.RDS', test)+'.sav'
pickle.dump(clf, open(filename, 'wb'))

end_time = time.process_time()
cpu_time = end_time - start_time

print(f"CPU time used: {cpu_time:.2f} seconds")