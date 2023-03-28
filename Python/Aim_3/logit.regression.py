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
from sklearn.inspection import permutation_importance
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

# Scale data. Required to speed up analysis with 'saga' solver
scaler = StandardScaler()
X_train_final = pd.DataFrame(scaler.fit_transform(X_train_final), columns=X_train_final.columns)
X_tune = pd.DataFrame(scaler.fit_transform(X_tune), columns=X_tune.columns)
X_test = pd.DataFrame(scaler.transform(X_test), columns=X_test.columns)

# Create the recursive feature eliminator that scores features by mean squared errors
clf = LogisticRegression(solver='saga', penalty='elasticnet', l1_ratio=0.5, max_iter=10000, random_state=42, n_jobs=-1)
rfecv = RFECV(clf, cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=42), scoring='accuracy', n_jobs=-1)
# Fit the RFECV object to the training data
rfecv = rfecv.fit(X_tune, y_tune)
print('Model training complete')
print('Optimal number of features: ', rfecv.n_features_)
print('Best features: ', rfecv.get_feature_names_out().tolist())

features = rfecv.get_feature_names_out().tolist()

# Tune the model to find the optimal C parameter
param_grid = {'C': [0.001, 0.01, 0.1, 1, 10]}
clf = LogisticRegression(solver='saga', penalty='elasticnet', l1_ratio=0.5, max_iter=10000, random_state=42, n_jobs=-1)
grid_search = GridSearchCV(clf, param_grid, cv=RepeatedKFold(n_splits=len(X_tune.index), n_repeats=3, random_state=0), scoring='accuracy', n_jobs=-1)
grid_search.fit(X_tune.loc[:, features], y_tune)

clf = grid_search.best_estimator_

# # Permute features and calculate feature importance
# if rfecv.n_features_ == 1:
#     # Fit the model
#     clf.fit(X_train_final.loc[:, rfecv.support_], y_train_final)
#     # Predict the test set
#     y_pred = clf.predict(X_test.iloc[:, rfecv.support_])
# else:
#     # Fit the model
#     clf.fit(X_train_final.loc[:, rfecv.support_], y_train_final)
#     # Predict the test set
#     y_pred = clf.predict(X_test.iloc[:, rfecv.support_])
#     r = permutation_importance(clf, X_train_final.loc[:, rfecv.support_], y_train_final,
#                         n_repeats=30,
#                         random_state=42,
#                         n_jobs=-1)

clf.fit(X_train_final.loc[:, features], y_train_final)
# Predict the test set
y_pred = clf.predict(X_test.loc[:, features])

# # Identify which features improve the model  
# selected_features = []
# # Loop over each feature in order of importance
# for i in r.importances_mean.argsort():
#     # Remove the i-th feature from the data
#     X_train_new = np.delete(X_train, i, axis=1)
#     X_val_new = np.delete(X_val, i, axis=1)

#     # Train a new model using the remaining features
#     model_new = LogisticRegression()
#     model_new.fit(X_train_new, y_train)

#     # Predict on the validation set and calculate the F1 score
#     y_pred_new = model_new.predict(X_val_new)
#     f1_new = f1_score(y_val, y_pred_new)

#     # Compare the F1 score to the baseline and add the feature to the list if it improves the score
#     if f1_new > f1_full:
#         selected_features.append(i)

# clf.fit(X_train_final.loc[:, new_features], y_train_final)
# y_pred = clf.predict(X_test.loc[:, new_features])

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
plt.title('logit: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
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
pickle.dump(clf, open(filename, 'wb'))

end_time = time.process_time()
cpu_time = end_time - start_time

print(f"CPU time used: {cpu_time:.2f} seconds")