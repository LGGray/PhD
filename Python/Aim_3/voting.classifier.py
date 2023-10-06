import sys
import pickle
import os.path
import pandas as pd
import numpy as np
import pyreadr
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score, cohen_kappa_score
import matplotlib.pyplot as plt
from sklearn.ensemble import VotingClassifier
from sklearn.inspection import permutation_importance

file = sys.argv[1]

# Read in expression RDS file
df = pyreadr.read_r(file)
df = df[None]
print(df.head())

cell = file.replace('psuedobulk/', '').replace('.RDS', '')

# load the model from disk
logit = pickle.load(open('psuedobulk/ML.models/logit_model_'+cell+'.sav', 'rb'))
RF = pickle.load(open('psuedobulk/ML.models/RF_model_'+cell+'.sav', 'rb'))
SVM = pickle.load(open('psuedobulk/ML.models/SVM_model_'+cell+'.sav', 'rb'))
GBM = pickle.load(open('psuedobulk/ML.models/GBM_model_'+cell+'.sav', 'rb'))
MLP = pickle.load(open('psuedobulk/ML.models/MLP_model_'+cell+'.sav', 'rb'))

# Read in tune, train, test and features
X_train = pd.read_csv('psuedobulk/data.splits/X_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_train = pd.read_csv('psuedobulk/data.splits/y_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
X_test = pd.read_csv('psuedobulk/data.splits/X_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_test = pd.read_csv('psuedobulk/data.splits/y_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
features = pd.read_csv('psuedobulk/features/enet_features.'+os.path.basename(file).replace('.RDS', '')+'.csv')

# Create a voting classifier
voting_clf = VotingClassifier(estimators=[('logit', logit), ('RF', RF), ('SVM', SVM), ('GBM', GBM), ('MLP', MLP)], voting='soft')

# Fit the voting classifier
voting_clf.fit(X_train.loc[:, features.iloc[:,0]], y_train['class'])

# Predict the test set
y_pred = voting_clf.predict(X_test.loc[:, features.iloc[:,0]])
# Get the predicted probabilities
y_pred_proba = voting_clf.predict_proba(X_test.loc[:, features.iloc[:,0]])[:, 1]

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
                        'AUC': [auc],
                        'Kappa': [kappa]})
metrics.to_csv('psuedobulk/ML.models/ensemble/metrics_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Save confusion matrix to file
confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
confusion.to_csv('psuedobulk/ML.models/ensemble/confusion_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Print the PR curve
precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
average_precision = average_precision_score(y_test, y_pred_proba)
disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
disp.plot()
disp.ax_.set_title('ensemble: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.savefig('psuedobulk/ML.models/ensemble/PRcurve_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Save the model
import pickle
filename = 'psuedobulk/ML.models/ensemble/'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(voting_clf, open(filename, 'wb'))


### Feature permutation importance ###
result = permutation_importance(voting_clf, X_test.loc[:, features.iloc[:,0]], y_test['class'], n_repeats=5, 
                           random_state=42, n_jobs=8, scoring=['f1', 'roc_auc'])
# Print the feature importance for F1
F1_features = []
for i in result['f1']['importances_mean'].argsort()[::-1]:
    if result['f1']['importances_mean'][i] - 2 * result['f1']['importances_std'][i] > 0:
        F1_features.append(X_test.loc[:, features.iloc[:,0]].columns[i])
        print(f"{X_test.loc[:, features.iloc[:,0]].columns[i]:<8}"
              f"{result['f1']['importances_mean'][i]:.3f}"
              f" +/- {result['f1']['importances_std'][i]:.3f}")
        
# Print the feature importance for AUC
auc_features = []
for i in result['roc_auc']['importances_mean'].argsort()[::-1]:
    if result['roc_auc']['importances_mean'][i] - 2 * result['roc_auc']['importances_std'][i] > 0:
        auc_features.append(X_test.loc[:, features.iloc[:,0]].columns[i])
        print(f"{X_test.loc[:, features.iloc[:,0]].columns[i]:<8}"
              f"{result['roc_auc']['importances_mean'][i]:.3f}"
              f" +/- {result['roc_auc']['importances_std'][i]:.3f}")

# Intersection of F1 and AUC features
important_features = set(F1_features + auc_features)

# Retrain model with important features
voting_clf.fit(X_train.loc[:, important_features], y_train['class'])

# Predict the test set
y_pred = voting_clf.predict(X_test.loc[:, important_features])
# Get the predicted probabilities
y_pred_proba = voting_clf.predict_proba(X_test.loc[:, important_features])[:, 1]

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
                        'AUC': [auc],
                        'Kappa': [kappa]})
metrics.to_csv('psuedobulk/ML.models/ensemble/perm_metrics_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Save confusion matrix to file
confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
confusion.to_csv('psuedobulk/ML.models/ensemble/perm_confusion_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Print the PR curve
precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
average_precision = average_precision_score(y_test, y_pred_proba)
disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
disp.plot()
disp.ax_.set_title('ensemble permuted features: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.savefig('psuedobulk/ML.models/ensemble/perm_PRcurve_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Save the model
import pickle
filename = 'psuedobulk/ML.models/ensemble/perm.'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(voting_clf, open(filename, 'wb'))