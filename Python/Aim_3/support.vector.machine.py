# Import required libraries
import sys
import os.path
import pandas as pd
import time
import pyreadr
from sklearn.model_selection import RepeatedKFold
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score
from sklearn.feature_selection import RFECV
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
df['class'] = df['class'].replace({"control": 0, "disease": 1})

# Split the data into features (X) and target (y)
X = df.iloc[:, 1:]
y = df.iloc[:, 0]

# Split the original dataset into a training set and a test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)

# Use StratifiedShuffleSplit() to split the training set into two additional subsets: 
# a subset for parameter tuning and a subset for final testing
cv = StratifiedShuffleSplit(n_splits=10, test_size=0.25, random_state=42)
train_index, tune_index = next(cv.split(X_train, y_train))

# Get the training and parameter tuning sets
X_train_final, y_train_final = X_train.iloc[train_index,], y_train.iloc[train_index]
X_tune, y_tune = X_train.iloc[tune_index], y_train.iloc[tune_index]

# Scale data to have min of 0 and max of 1. Required for SVM
scaler = MinMaxScaler()
X_train_final = pd.DataFrame(scaler.fit_transform(X_train_final), columns=X_train_final.columns)
X_tune = pd.DataFrame(scaler.fit_transform(X_tune), columns=X_tune.columns)
X_test = pd.DataFrame(scaler.transform(X_test), columns=X_test.columns)

# Build linear SVM for feature selection
param_grid = {'C': [0.1, 1, 10, 100, 1000]}
clf = SVC(kernel='linear', probability=True, max_iter=10000, random_state=42)
grid_search = GridSearchCV(clf, param_grid, cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0), n_jobs=-1, verbose=1)
grid_search.fit(X_tune, y_tune)
rfecv = RFECV(grid_search.best_estimator_, cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0), scoring='accuracy', n_jobs=-1)
rfecv = rfecv.fit(X_train_final, y_train_final)

# Perform a grid search to find the best parameters
# Create the parameter grid
param_grid = {'kernel': ['linear', 'rbf', 'poly', 'sigmoid'],
                'C': [0.1, 1, 10, 100, 1000]
}
clf = SVC(probability=True, max_iter=10000, random_state=42)
grid_search = GridSearchCV(clf, param_grid, cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0), n_jobs=-1, verbose=1)
grid_search.fit(X_tune.loc[:, rfecv.support_], y_tune)
# Get the best estimator with the optimal hyperparameters
clf = SVC(probability=True, max_iter=10000, random_state=42,
          kernel=grid_search.best_params_['kernel'],
          C=grid_search.best_params_['C'])
# Fit model
clf.fit(X_train_final.loc[:, rfecv.support_], y_train_final)
# Predict on test data
y_pred = clf.predict(X_test.iloc[:, rfecv.support_])

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
fpr, tpr, thresholds = roc_curve(y_test, y_pred)
plt.figure()
plt.plot(fpr, tpr, label='AUC-ROC (area = %0.2f)' % auc)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('logit: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.legend(loc="lower right")
plt.savefig('exp.matrix/AUROC/SVM_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Print the PR curve
precision, recall, thresholds = precision_recall_curve(y_test, y_pred)
average_precision = average_precision_score(y_test, y_pred)
disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
disp.plot()
disp.ax_.set_title('logit: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.savefig('exp.matrix/PRC/SVM_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Save the model
import pickle
filename = 'ML.models/SVM_model_'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(clf, open(filename, 'wb'))

end_time = time.process_time()
cpu_time = end_time - start_time

print(f"CPU time used: {cpu_time:.2f} seconds")