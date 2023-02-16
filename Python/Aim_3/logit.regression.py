# Import required libraries
import sys
import os.path
import pandas as pd
import pyreadr
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import LogisticRegression
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

# Scale data. Required to speed up analysis with 'saga' solver
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

# Build the logistical regression model using the saga solver and elasticnet penalty
# Create the recursive feature eliminator that scores features by mean squared errors
clf = LogisticRegression(solver='saga', penalty='elasticnet', l1_ratio=0.5, max_iter=10000, n_jobs=-1)
rfecv = RFECV(clf, cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0), scoring='accuracy', n_jobs=-1)
# Fit the RFECV object to the training data
rfecv = rfecv.fit(X_train, y_train)
print('Model training complete')
print('Optimal number of features: ', rfecv.n_features_)
# Fit the model
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

# Save the model
import pickle
filename = 'ML.models/logit_model_'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(rfecv, open(filename, 'wb'))