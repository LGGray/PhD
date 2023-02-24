# Import required libraries
import sys
import os.path
import time
import pandas as pd
import pyreadr
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import LogisticRegression
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

X_train.iloc[train_index]

# Get the training and parameter tuning sets
X_train_final, y_train_final = X_train.iloc[train_index,], y_train.iloc[train_index]
X_tune, y_tune = X_train.iloc[tune_index], y_train.iloc[tune_index]

# Scale data. Required to speed up analysis with 'saga' solver
scaler = StandardScaler()
X_train_final = pd.DataFrame(scaler.fit_transform(X_train_final), columns=X_train_final.columns)
X_tune = pd.DataFrame(scaler.fit_transform(X_tune), columns=X_tune.columns)
X_test = pd.DataFrame(scaler.transform(X_test), columns=X_test.columns)

# Tune the model to find the optimal C parameter
param_grid = {'C': [0.001, 0.01, 0.1, 1, 10]}
clf = LogisticRegression(solver='saga', penalty='elasticnet', l1_ratio=0.5, max_iter=10000, random_state=42, n_jobs=-1)
grid_search = GridSearchCV(clf, param_grid, cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0), scoring='accuracy', n_jobs=-1)
grid_search.fit(X_tune, y_tune)

# Build the logistical regression model using the saga solver and elasticnet penalty
# Create the recursive feature eliminator that scores features by mean squared errors
clf = LogisticRegression(solver='saga', C=grid_search.best_params_['C'], penalty='elasticnet', l1_ratio=0.5, max_iter=10000, random_state=42, n_jobs=-1)
rfecv = RFECV(clf, cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0), scoring='accuracy', n_jobs=-1)
# Fit the RFECV object to the training data
rfecv = rfecv.fit(X_train_final, y_train_final)
print('Model training complete')
print('Optimal number of features: ', rfecv.n_features_)
print('Best features: ', rfecv.get_feature_names_out().tolist())

# Fit the model
clf.fit(X_train_final.loc[:, rfecv.support_], y_train_final)
# Predict the test set
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
metrics.to_csv('exp.matrix/metrics/logit_metrics_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Save confusion matrix to file
confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
confusion.to_csv('exp.matrix/metrics/logit_confusion_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Print the AUROC curve
fpr, tpr, thresholds = roc_curve(y_test, y_pred)
plt.plot(fpr, tpr, label='AUC-ROC (area = %0.2f)' % auc)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc="lower right")
plt.savefig('exp.matrix/AUROC/logit_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Print the PR curve
precision, recall, thresholds = precision_recall_curve(y_test, y_pred)
average_precision = average_precision_score(y_test, y_pred)
disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
disp.plot()
disp.ax_.set_title('logit: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.savefig('exp.matrix/PRC/logit_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Save the model
import pickle
filename = 'ML.models/logit_model_'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(rfecv, open(filename, 'wb'))

end_time = time.process_time()
cpu_time = end_time - start_time

print(f"CPU time used: {cpu_time:.2f} seconds")