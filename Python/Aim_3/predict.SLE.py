import sys
import pickle
import os.path
import pandas as pd
import numpy as np
import pyreadr
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score, cohen_kappa_score
import matplotlib.pyplot as plt
from sklearn.ensemble import VotingClassifier
from sklearn.inspection import permutation_importance
from sklearn.utils import resample

# Load in ensemble model
model_path = f'new_pseudobulk/split_{sys.argv[1]}/intersection/ensemble/{sys.argv[2]}'
eclf = pickle.load(open(model_path, 'rb'))
features = eclf.feature_names_in_

cell = {sys.argv[2]}.replace('.sav', '')

### Read in independent test set
test = pyreadr.read_r(f'/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/SLE_GSE135779/psuedobulk/{cell}.RDS')
test = test[None]

# Subset for adult, child or all
if sys.argv[3] == 'adult':
    # Subset class to include aHD and aSLE
    test = test[test['class'].isin(['aHD', 'aSLE'])]
if sys.argv[3] == 'child':
    # Subset class to include cHD and cSLE
    test = test[test['class'].isin(['cHD', 'cSLE'])]
if sys.argv[3] == 'all':
    # Subset class to include all
    test = test[test['class'].isin(['cHD', 'cSLE', 'aHD', 'aSLE'])]

# Replace class labels with 0 and 1
test['class'] = test['class'].replace({"cSLE": 1, "aSLE": 1, "cHD": 0, "aHD": 0})

# Set class as target variable
y_test = test[['class']]

# Check which features are missing and fill with 0
missing = np.setdiff1d(features, test.columns)
test[missing] = 0
X_test = test[features]

# Scale the data
scaler = StandardScaler()
X_test = pd.DataFrame(scaler.fit_transform(X_test), columns=X_test.columns, index=X_test.index)

# Predict the test set
y_pred = eclf.predict(X_test)
y_pred_proba = eclf.predict_proba(X_test)[:, 1]

# Calculate Youden's J statistic to find the optimal threshold
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
j_scores = tpr - fpr
# Find the optimal threshold
optimal_idx = np.argmax(j_scores)
optimal_threshold = thresholds[optimal_idx]
# Convert probabilities to binary predictions based on optimal threshold
y_pred = (y_pred_proba >= optimal_threshold).astype(int)

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
                        'Kappa': [kappa],
                        'n_features': [len(features)]})

print(metrics)

# metrics.to_csv(f'new_pseudobulk/split_1/intersection/ensemble/metrics_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# # Save confusion matrix to file
# confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
# confusion.to_csv(f'new_pseudobulk/split_1/intersection/ensemble/confusion_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# # Define class names
# classes = ['Control', 'Disease']
# fig, ax = plt.subplots()
# # Set the color map to 'coolwarm'
# cmap = plt.cm.coolwarm
# # Create the heatmap for the confusion matrix
# cax = ax.matshow(confusion, cmap=cmap)
# # Add color bar
# plt.colorbar(cax)
# # Add counts to the confusion matrix cells
# confusion_values = confusion.values
# for (i, j), val in np.ndenumerate(confusion_values):
#     ax.text(j, i, f'{val}', ha='center', va='center', color='white')
# # Set axis labels
# ax.set_xlabel('Predicted labels')
# ax.set_ylabel('True labels')
# ax.set_xticks(range(len(classes)))
# ax.set_yticks(range(len(classes)))
# ax.set_xticklabels(classes)
# ax.set_yticklabels(classes)
# # Set the title
# ax.set_title('Ensemble: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
# # Annotate with F1 score
# plt.annotate(f'F1 Score: {f1:.2f}', xy=(0.5, -0.1), xycoords='axes fraction', 
#              ha='center', va='center', fontsize=12, color='black')
# # Adjust layout for visibility
# plt.tight_layout()
# # Save the figure
# plt.savefig(f'new_pseudobulk/split_1/intersection/ensemble/confusion_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')
# plt.close()

# # Print the PR curve
# precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
# average_precision = average_precision_score(y_test, y_pred_proba)
# disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
# disp.plot()
# disp.ax_.set_title('Ensemble: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
# plt.savefig(f'new_pseudobulk/split_1/intersection/ensemble/PRcurve_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# # Save the model
# import pickle
# filename = f'new_pseudobulk/split_1/intersection/ensemble/'+os.path.basename(file).replace('.RDS', '')+'.sav'
# pickle.dump(voting_clf, open(filename, 'wb'))

# # Save features to file
# if sys.argv[3] == 'intersection':
#     features.to_csv(f'new_pseudobulk/split_1/features/intersection_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)
# elif sys.argv[3] == 'combined':
#     features.to_csv(f'new_pseudobulk/split_1/features/combined_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)