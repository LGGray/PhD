import pickle
import pyreadr
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist
import glob
import os.path
import sys
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score, cohen_kappa_score
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler, MinMaxScaler


model_path = sys.argv[1]
model = os.path.basename(model_path).replace('.sav', '')

# Load in model
eclf = pickle.load(open(model_path, 'rb'))
features = eclf.feature_names_in_

# Read in expression RDS file
exp_path = sys.argv[2]
exp = pyreadr.read_r(exp_path)
exp = exp[None]
print(exp.head())
scRNA = os.path.basename(exp_path).replace('.RDS', '')

# Replace class labels with 0 and 1
exp['class'] = exp['class'].replace({"cSLE": 1, "aSLE": 1, "cHD": 0, "aHD": 0})

# Check which features are missing in exp
missing = np.setdiff1d(features, exp.columns)
exp[missing] = 0

print('Missing features:', missing)

y_test = exp['class']
X_test = exp[features]

scaler = StandardScaler()
X_test = pd.DataFrame(scaler.fit_transform(X_test), columns=X_test.columns)

# Predict
y_pred = eclf.predict(X_test)
y_pred_proba = eclf.predict_proba(X_test)[:,1]

# Calculate the metrics
# Calculate the metrics
accuracy = accuracy_score(y_test, y_pred)
precision = precision_score(y_test, y_pred)
recall = recall_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred, average='weighted')
auroc = roc_auc_score(y_test, y_pred)
auprc = average_precision_score(y_test, y_pred)
kappa = cohen_kappa_score(y_test, y_pred)

# Print the metrics
metrics = pd.DataFrame({'Accuracy': [accuracy],
                        'Precision': [precision],
                        'Recall': [recall],
                        'F1': [f1],
                        'AUC': [auroc],
                        'AUPRC': [auprc],
                        'Kappa': [kappa]})
print(metrics)

metrics.to_csv('psuedobulk/scRNA/'+'metrics_'+model+'.csv', index=False)