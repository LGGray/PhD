# Import required libraries
import sys
import pandas as pd
import pyreadr
import multiprocessing
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import confusion_matrix
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc, roc_auc_score
import matplotlib.pyplot as plt

# Get the file name from the command line
file = sys.argv[0]

# Get cell name
cell = file.replace('.RDS', '')

# Read in differentially expressed genes
deg = pd.read_csv('psuedobulk/Fibroblasts.edgeR-LRT.txt', delimiter='\t')
# filter deg for FDR < 0.05 and |logFC| > 0.5
deg = deg[(deg['FDR'] < 0.05) & (abs(deg['logFC']) > 0.5)]

# Read in expression RDS file
df = pyreadr.read_r(file)
df = df[None]

# Replace classes with binary label
df['class'] = df['class'].replace({"control": 0, "disease": 1})

# Filter df column by deg['gene'] and 'class'
cols_to_keep = ['class'] + deg['gene'].tolist()
df = df.loc[:,df.columns.isin(cols_to_keep)]

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

# Build the logistical regression model
# Create the recursive feature eliminator that scores features by mean squared errors
clf = LogisticRegression(solver='sag', penalty='l2', C=0.25, max_iter=10000)
rfecv = RFECV(clf, cv=5, scoring='accuracy', n_jobs=4)

# Fit the RFECV object to the training data
X_train_selected = rfecv.fit_transform(X_train, y_train)
# Transform the test data
X_test_selected = rfecv.transform(X_test)
# Fit the model
y_pred = clf.predict(X_test_selected)

accuracy = accuracy_score(y_test, y_pred)
precision = precision_score(y_test, y_pred)
recall = recall_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred)
auc = roc_auc_score(y_test, y_pred)

# Print the results
print('Accuracy: ', accuracy)
print('Precision: ', precision)
print('Recall: ', recall)
print('F1: ', f1)
print('AUC: ', auc)

# Print the confusion matrix
print(confusion_matrix(y_test, y_pred))

# Print the ROC curve
fpr, tpr, thresholds = roc_curve(y_test, y_pred)
roc_auc = auc(fpr, tpr)
plt.figure()
plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.show()
# save the figure
plt.savefig('ROC_curve_Fibroblasts.png')


