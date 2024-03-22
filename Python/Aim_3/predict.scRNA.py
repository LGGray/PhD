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
from sklearn.utils import resample


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

if sys.argv[3] == 'adult':
    # Subset class to include aHD and aSLE
    exp = exp[exp['class'].isin(['aHD', 'aSLE'])]
if sys.argv[3] == 'child':
    # Subset class to include cHD and cSLE
    exp = exp[exp['class'].isin(['cHD', 'cSLE'])]
if sys.argv[3] == 'all':
    # Subset class to include all
    exp = exp[exp['class'].isin(['cHD', 'cSLE', 'aHD', 'aSLE'])]

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

# Define bootstrap parameters
n_bootstraps = 1000
confidence_level = 0.9
# Initialize an empty list to store bootstrap scores
bootstrapped_f1 = []
bootstrapped_AUC = []
bootstrapped_AUPRC = []

# Loop over bootstrap samples
for i in range(n_bootstraps):
    # Resample with replacement
    y_test_resampled, y_pred_resampled = resample(y_test, y_pred, stratify=y_test)
    # Calculate F1 score
    f1 = f1_score(y_test_resampled, y_pred_resampled, average='weighted')
    # Calculate AUROC
    auroc = roc_auc_score(y_test_resampled, y_pred_resampled)
    # Calculate AUPRC
    auprc = average_precision_score(y_test_resampled, y_pred_resampled)
    # Append score to list
    bootstrapped_f1.append(f1)
    bootstrapped_AUC.append(auroc)
    bootstrapped_AUPRC.append(auprc)

# Calculate percentile for confidence intervals
lower_percentile = (1 - confidence_level) / 2 * 100
upper_percentile = (1 + confidence_level) / 2 * 100

f1_lower_bound = np.percentile(bootstrapped_f1, lower_percentile)
f1_upper_bound = np.percentile(bootstrapped_f1, upper_percentile)

auroc_lower_bound = np.percentile(bootstrapped_AUC, lower_percentile)
auroc_upper_bound = np.percentile(bootstrapped_AUC, upper_percentile)

auprc_lower_bound = np.percentile(bootstrapped_AUPRC, lower_percentile)
auprc_upper_bound = np.percentile(bootstrapped_AUPRC, upper_percentile)

# Create dataframe of metrics and save to file
metrics = pd.DataFrame({'Accuracy': [accuracy], 
                        'Precision': [precision], 
                        'Recall': [recall], 
                        'F1': [f1],
                        'F1_lower': [f1_lower_bound],
                        'F1_upper': [f1_upper_bound],
                        'AUC': [auroc],
                        'AUC_lower': [auroc_lower_bound],
                        'AUC_upper': [auroc_upper_bound],
                        'AUPRC': [auprc],
                        'AUPRC_lower': [auprc_lower_bound],
                        'AUPRC_upper': [auprc_upper_bound],
                        'Kappa': [kappa]})

if sys.argv[3] == 'adult':
    metrics.to_csv('psuedobulk/scRNA/'+'metrics_'+model+'_adult.csv', index=False)
if sys.argv[3] == 'child':
    metrics.to_csv('psuedobulk/scRNA/'+'metrics_'+model+'_child.csv', index=False)
if sys.argv[3] == 'all':
    metrics.to_csv('psuedobulk/scRNA/'+'metrics_'+model+'_all.csv', index=False)

# Create confusion matrix
confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
# Define class names
classes = ['Control', 'Disease']
fig, ax = plt.subplots()
# Set the color map to 'coolwarm'
cmap = plt.cm.coolwarm
# Create the heatmap for the confusion matrix
cax = ax.matshow(confusion, cmap=cmap)
# Add color bar
plt.colorbar(cax)
# Add counts to the confusion matrix cells
confusion_values = confusion.values
for (i, j), val in np.ndenumerate(confusion_values):
    ax.text(j, i, f'{val}', ha='center', va='center', color='white', fontsize=12)
# Set axis labels
ax.set_xlabel('Predicted labels')
ax.set_ylabel('True labels')
ax.set_xticks(range(len(classes)))
ax.set_yticks(range(len(classes)))
ax.set_xticklabels(classes)
ax.set_yticklabels(classes)
# Set the title
ax.set_title(os.path.basename(sys.argv[2]).replace('.RDS', '').replace('.', ' '))
# Annotate with F1 score
plt.annotate(f'F1 Score: {f1:.2f}', xy=(0.5, -0.1), xycoords='axes fraction', 
             ha='center', va='center', fontsize=12, color='black')
# Adjust layout for visibility
plt.tight_layout()
# Save the figure
plt.savefig('ML.plots/confusion_'+os.path.basename(sys.argv[2]).replace('.RDS', '_')+sys.argv[3]+'.pdf', bbox_inches='tight')
plt.close()
