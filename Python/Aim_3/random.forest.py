# Import required libraries
import sys
import os.path
import pandas as pd
import pyreadr
from sklearn.model_selection import RepeatedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score
from sklearn.feature_selection import RFECV
from sklearn.metrics import roc_curve, auc, roc_auc_score
import matplotlib.pyplot as plt

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

# Create the stratified sampling object
sss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)

# Loop through the splits
for train_index, test_index in sss.split(X, y):
    # Get the training and testing data
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = y.iloc[train_index], y.iloc[test_index]

# Perform a grid search to find the best parameters
# Create the parameter grid
param_grid = {'n_estimators': [100, 200, 300, 400],
                'max_depth': [5, 10, 15, 30],
                'min_samples_split': [2, 5, 8, 10]
}
clf = RandomForestClassifier()
grid_search = GridSearchCV(clf, param_grid, cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0), n_jobs=-1, verbose=1)

# Fit the grid search object to the training data
grid_search.fit(X_train, y_train)

# Create an RFECV object with a random forest classifier
clf = RandomForestClassifier(n_estimators=grid_search.best_params_['n_estimators'], 
                            max_depth=grid_search.best_params_['max_depth'], 
                            min_samples_split=grid_search.best_params_['min_samples_split'], n_jobs=-1)
rfecv = RFECV(clf, cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0), scoring='accuracy', n_jobs=-1)

# Fit the RFECV object to the training data
rfecv = rfecv.fit(X_train, y_train)
print('Model training complete')
print('Optimal number of features: ', rfecv.n_features_)

y_pred = rfecv.predict(X_test)

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
metrics.to_csv('exp.matrix/metrics/RF_metrics_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Save confusion matrix to file
confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
confusion.to_csv('exp.matrix/metrics/RF_confusion_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Print the AUROC curve
fpr, tpr, thresholds = roc_curve(y_test, y_pred)
plt.figure()
plt.plot(fpr, tpr, label='AUC-ROC (area = %0.2f)' % auc)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc="lower right")
plt.savefig('exp.matrix/AUROC/RF_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Save the model
import pickle
filename = 'ML.models/RF_model_'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(rfecv, open(filename, 'wb'))