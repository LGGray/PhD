import pickle
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

model_path = 'psuedobulk/ML.models/ensemble/Non.classical.monocytes.chrX.sav'
# model_path = sys.argv[1]
model = os.path.basename(model_path).replace('.sav', '')

# Load in model
eclf = pickle.load(open(model_path, 'rb'))
features = eclf.feature_names_in_

# Read in bulk RNA
exp_path = '/directflow/SCCGGroupShare/projects/lacgra/bulk.data/E-GEOD-72509/E-GEOD-72509.female.csv'
# exp_path = sys.argv[2]
exp = pd.read_csv(exp_path, sep=',')

# Replace class labels with 0 and 1
exp['class'] = exp['class'].replace({"control": 0, "disease": 1})

# Check which features are missing in exp
missing = np.setdiff1d(features, exp.columns)
exp[missing] = 0

y_test = exp['class']
X_test = exp[features]

scaler = StandardScaler()
X_test = pd.DataFrame(scaler.fit_transform(X_test), columns=X_test.columns)

# Predict
y_pred = eclf.predict(X_test)
y_pred_proba = eclf.predict_proba(X_test)[:,1]

# Calculate the metrics
accuracy = accuracy_score(y_test, y_pred)
precision = precision_score(y_test, y_pred)
recall = recall_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred)
auc = roc_auc_score(y_test, y_pred)
kappa = cohen_kappa_score(y_test, y_pred)

# Print the metrics
metrics = pd.DataFrame({'accuracy': [accuracy],
                        'precision': [precision],
                        'recall': [recall],
                        'f1': [f1],
                        'auc': [auc],
                        'kappa': [kappa]})
print(metrics)

Memory B cells:
accuracy  precision    recall        f1       auc     kappa
0       0.5   0.666667  0.571429  0.615385  0.452381 -0.086957
Treg:
   accuracy  precision    recall        f1      auc     kappa
0       0.6   0.714286  0.714286  0.714286  0.52381  0.047619
Th:
   accuracy  precision    recall        f1       auc     kappa
0       0.7   0.785714  0.785714  0.785714  0.642857  0.285714
nCM:
   accuracy  precision    recall        f1       auc     kappa
0       0.7   0.785714  0.785714  0.785714  0.642857  0.285714


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

metrics.to_csv('psuedobulk/bulk_RNA/'+'metrics_'+model+'_'+bulk+'.csv', index=False)

# Save confusion matrix to file
confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
confusion.to_csv('psuedobulk/bulk_RNA/'+'confusion_'+model+'_'+bulk+'.csv', index=False)

# Print the AUROC curve
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
plt.figure()
plt.plot(fpr, tpr, label='AUC-ROC (area = %0.2f)' % auc)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title(model+'_'+bulk.replace('_', ' ').replace('.', ' '))
plt.legend(loc="lower right")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])
plt.savefig('psuedobulk/bulk_RNA/AUROC/AUROC_'+model+'_'+bulk+'.pdf', bbox_inches='tight')

# Print the PR curve
precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
average_precision = average_precision_score(y_test, y_pred_proba)
disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
disp.plot()
disp.ax_.set_title(model+'_'+bulk.replace('_', ' ').replace('.', ' '))
plt.savefig('psuedobulk/bulk_RNA/PRC/PRC_'+model+'_'+bulk+'.pdf', bbox_inches='tight')