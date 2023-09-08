# Import required libraries
import sys
import os.path
import pandas as pd
import numpy as np
import time
import pyreadr
from boruta import BorutaPy
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import GridSearchCV, RepeatedKFold, GroupShuffleSplit
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score, cohen_kappa_score
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import roc_curve, auc, roc_auc_score
import matplotlib.pyplot as plt

start_time = time.process_time()

# Get the file name from the command line
file = sys.argv[1]
print(os.path.basename(file))

# Read in expression RDS file
df = pyreadr.read_r(file)
df = df[None]
print(df.head())

# Replace classes with binary label
if sum(df['class'] == 'control') > 0:
  df['class'] = df['class'].replace({"control": 0, "disease": 1})
else:
  df['class'] = df['class'].replace({"managed": 0, "flare": 1})

# ### Split the data into train, tune and test sets ###

# # Collect individual IDs
# individuals = df['individual'].unique()
# n_individuals = len(individuals)

# # Get the number of individuals in each condition
# individual_class = df['individual'].astype(str) + '_' + df['class'].astype(str)
# n_control = len(individual_class[individual_class.str.endswith('_0')].unique())
# n_disease = len(individual_class[individual_class.str.endswith('_1')].unique())

# # Determine number of controls and disease samples to include in each dataset
# n_test_control = int(n_control * 0.2)
# n_tune_control = int(n_control * 0.2)
# n_train_control = n_control - n_test_control - n_tune_control

# n_test_disease = int(n_disease * 0.2)
# n_tune_disease = int(n_disease * 0.2)
# n_train_disease = n_disease - n_test_disease - n_tune_disease

# # Randomly assign controls to each dataset
# test_control_individuals = np.random.choice(
#     df[df['class'] == 0]['individual'].unique(),
#     size=n_test_control,
#     replace=False
# )
# tune_control_individuals = np.random.choice(
#     np.setdiff1d(
#         df[df['class'] == 0]['individual'].unique(),
#         test_control_individuals
#     ),
#     size=n_tune_control,
#     replace=False
# )
# train_control_individuals = np.setdiff1d(
#     df[df['class'] == 0]['individual'].unique(),
#     np.concatenate([test_control_individuals, tune_control_individuals])
# )

# # Randomly assign disease samples to each dataset
# test_disease_individuals = np.random.choice(
#     df[df['class'] == 1]['individual'].unique(),
#     size=n_test_disease,
#     replace=False
# )
# tune_disease_individuals = np.random.choice(
#     np.setdiff1d(
#         df[df['class'] == 1]['individual'].unique(),
#         test_disease_individuals
#     ),
#     size=n_tune_disease,
#     replace=False
# )
# train_disease_individuals = np.setdiff1d(
#     df[df['class'] == 1]['individual'].unique(),
#     np.concatenate([test_disease_individuals, tune_disease_individuals])
# )

# # Get the corresponding cells for each dataset
# test_index = df['individual'].isin(np.concatenate([test_control_individuals, test_disease_individuals]))
# tune_index = df['individual'].isin(np.concatenate([tune_control_individuals, tune_disease_individuals]))
# train_index = df['individual'].isin(np.concatenate([train_control_individuals, train_disease_individuals]))

# # Split data into training, tuning, and testing sets
# X_train, X_test, X_tune = df.loc[train_index,].drop(['class', 'individual'], axis=1), df.loc[test_index,].drop(['class', 'individual'], axis=1), df.loc[tune_index,].drop(['class', 'individual'], axis=1)
# y_train, y_test, y_tune = df.loc[train_index, 'class'], df.loc[test_index, 'class'], df.loc[tune_index, 'class']

# ### Boruta feature selection ###
# X = X_tune.values
# y = y_tune.ravel()
# # random forest classifier utilising all cores and sampling in proportion to y labels
# param_grid = {'n_estimators': [100, 200, 300, 400],
#               'max_features': ['sqrt', 'log2', 0.3],
#                 'max_depth': [5, 10, 15, 30],
#                 'min_samples_split': [2, 5, 8, 10]
# }
# clf = RandomForestClassifier(n_jobs=-1)
# grid_search = GridSearchCV(clf, param_grid, cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0), n_jobs=-1, verbose=1)
# # Fit the grid search object to the training data
# grid_search.fit(X, y)
# # Create an RFECV object with a random forest classifier
# rf = RandomForestClassifier(n_estimators=grid_search.best_params_['n_estimators'], 
#                             max_depth=grid_search.best_params_['max_depth'], 
#                             min_samples_split=grid_search.best_params_['min_samples_split'], n_jobs=-1)
# # define Boruta feature selection method
# feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, random_state=1)
# # find all relevant features - 5 features should be selected
# feat_selector.fit(X, y)
# # Return features
# features = X_tune.columns[feat_selector.support_].tolist()

# # Scale data to have min of 0 and max of 1. Required for SVM
# scaler = MinMaxScaler()
# X_train = pd.DataFrame(scaler.fit_transform(X_train), columns=X_train.columns)
# X_tune = pd.DataFrame(scaler.fit_transform(X_tune), columns=X_tune.columns)
# X_test = pd.DataFrame(scaler.fit_transform(X_test), columns=X_test.columns)

# Read in tune, train, test and features
X_tune = pd.read_csv('data.splits/X_tune.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_tune = pd.read_csv('data.splits/y_tune.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
X_train = pd.read_csv('data.splits/X_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_train = pd.read_csv('data.splits/y_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
X_test = pd.read_csv('data.splits/X_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_test = pd.read_csv('data.splits/y_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
features = pd.read_csv('data.splits/'+os.path.basename(file).replace('.RDS', '')+'.csv')

# Perform a grid search to find the best parameters
# Create the parameter grid
param_grid = {'kernel': ['linear', 'rbf', 'poly', 'sigmoid'],
                'C': [0.1, 1, 10, 100, 1000]
}
clf = SVC(probability=True, max_iter=20000, random_state=42)
grid_search = GridSearchCV(clf, param_grid, cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0), n_jobs=-1, verbose=1)
grid_search.fit(X_tune.loc[:, features.iloc[:,0]], y_tune['class'])
# Get the best estimator with the optimal hyperparameters
clf = SVC(probability=True, max_iter=20000, random_state=42,
          kernel=grid_search.best_params_['kernel'],
          C=grid_search.best_params_['C'])

# Fit the model
clf.fit(X_train.loc[:, features.iloc[:,0]], y_train['class'])

# Bootstrap to aggregate multiple models
bootstrapped_models = []
bootstrapped_f1 = []
for i in range(0, 5):
   gss = GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=i)
   groups = df.loc[X_train.index, 'individual']
   train, test = next(gss.split(X_train.loc[:, features.iloc[:,0]], y_train['class'], groups=groups))
   X_train = X_train.loc[:, features.iloc[:,0]]
   X_tr, X_te, y_tr, y_te = X_train.iloc[train], X_train.iloc[test], y_train.iloc[train], y_train.iloc[test]
   clf.fit(X_tr, y_tr['class'])
   # Store the trained model
   bootstrapped_models.append(clf)
   # Store the F1 score
   y_pred = clf.predict(X_te.loc[:, features.iloc[:,0]])
   f1 = f1_score(y_te, y_pred)
   bootstrapped_f1.append(f1)

# Ensemble classifier
from sklearn.ensemble import VotingClassifier
eclf = VotingClassifier(estimators=[('clf1', bootstrapped_models[0]), 
                                    ('clf2', bootstrapped_models[1]), 
                                    ('clf3', bootstrapped_models[2]), 
                                    ('clf4', bootstrapped_models[3]), 
                                    ('clf5', bootstrapped_models[4])], 
                                    voting='soft')
eclf.fit(X_train.loc[:, features.iloc[:,0]], y_train['class'])
y_pred = eclf.predict(X_test.loc[:, features.iloc[:,0]])
y_pred_proba = eclf.predict_proba(X_test.loc[:, features.iloc[:,0]])[:, 1]

# Calculate the metrics
accuracy = accuracy_score(y_test, y_pred)
precision = precision_score(y_test, y_pred)
recall = recall_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred)
auc = roc_auc_score(y_test, y_pred)
kappa = cohen_kappa_score(y_test, y_pred)

# Define bootstrap parameters
n_bootstraps = 1000
confidence_level = 0.9
# Initialize an empty list to store bootstrap scores
bootstrapped_scores = []

from sklearn.utils import resample
# Loop over bootstrap samples
for i in range(n_bootstraps):
    # Resample with replacement
    y_test_resampled, y_pred_resampled = resample(y_test, y_pred, stratify=y_test)
    # Calculate F1 score
    score = f1_score(y_test_resampled, y_pred_resampled)
    # Append score to list
    bootstrapped_scores.append(score)
# Sort the scores
sorted_scores = np.array(bootstrapped_scores)

# Calculate lower and upper bounds of confidence interval
alpha = (1 - confidence_level) / 2
lower_bound = sorted_scores[int(alpha * len(sorted_scores))]
upper_bound = sorted_scores[int((1 - alpha) * len(sorted_scores))]

# Create dataframe of metrics and save to file
metrics = pd.DataFrame({'Accuracy': [accuracy], 
                        'Precision': [precision], 
                        'Recall': [recall], 
                        'F1': [f1],
                        'F1_lower': [lower_bound],
                        'F1_upper': [upper_bound],
                        'AUC': [auc]})
metrics.to_csv('exp.matrix/metrics/SVM_metrics_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Save confusion matrix to file
confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
confusion.to_csv('exp.matrix/metrics/SVM_confusion_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Print the AUROC curve
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
plt.figure()
plt.plot(fpr, tpr, label='AUC-ROC (area = %0.2f)' % auc)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('SVM: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.legend(loc="lower right")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])
plt.savefig('exp.matrix/AUROC/SVM_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Print the PR curve
precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
average_precision = average_precision_score(y_test, y_pred_proba)
disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
disp.plot()
disp.ax_.set_title('SVM: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.savefig('exp.matrix/PRC/SVM_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Save the model
import pickle
filename = 'ML.models/SVM_model_'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(eclf, open(filename, 'wb'))

end_time = time.process_time()
cpu_time = end_time - start_time

print(f"CPU time used: {cpu_time:.2f} seconds")