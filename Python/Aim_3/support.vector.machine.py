# Import required libraries
import sys
import os.path
import pandas as pd
import numpy as np
import time
import pyreadr
from boruta import BorutaPy
from sklearn.model_selection import RepeatedKFold
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score
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

# Boruta feature selection
X = X_tune.values
y = y_tune.ravel()
# random forest classifier utilising all cores and sampling in proportion to y labels
rf = RandomForestClassifier(n_jobs=-1, class_weight='balanced', max_depth=5)
# define Boruta feature selection method
feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, random_state=1)
# find all relevant features - 5 features should be selected
feat_selector.fit(X, y)
# Return features
features = X_tune.columns[feat_selector.support_].tolist()

# Scale data to have min of 0 and max of 1. Required for SVM
scaler = MinMaxScaler()
X_train = pd.DataFrame(scaler.fit_transform(X_train), columns=X_train.columns)
X_tune = pd.DataFrame(scaler.fit_transform(X_tune), columns=X_tune.columns)
X_test = pd.DataFrame(scaler.fit_transform(X_test), columns=X_test.columns)

# Perform a grid search to find the best parameters
# Create the parameter grid
param_grid = {'kernel': ['linear', 'rbf', 'poly', 'sigmoid'],
                'C': [0.1, 1, 10, 100, 1000]
}
clf = SVC(probability=True, max_iter=20000, random_state=42)
grid_search = GridSearchCV(clf, param_grid, cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0), n_jobs=-1, verbose=1)
grid_search.fit(X_tune.loc[:, features], y_tune)
# Get the best estimator with the optimal hyperparameters
clf = SVC(probability=True, max_iter=20000, random_state=42,
          kernel=grid_search.best_params_['kernel'],
          C=grid_search.best_params_['C'])
# Fit model
clf.fit(X_train.loc[:, features], y_train)
# Predict on test data
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
pickle.dump(clf, open(filename, 'wb'))

end_time = time.process_time()
cpu_time = end_time - start_time

print(f"CPU time used: {cpu_time:.2f} seconds")