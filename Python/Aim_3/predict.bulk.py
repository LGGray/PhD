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

model_path = 'psuedobulk/ML.models/ensemble/perm.all.cells.chrX.sav'
# model_path = sys.argv[1]
model = os.path.basename(model_path).replace('.sav', '')

# Load in model
eclf = pickle.load(open(model_path, 'rb'))
features = eclf.feature_names_in_

# Read in bulk RNA and metadata
exp_path = '/directflow/SCCGGroupShare/projects/lacgra/bulk.data/GSE108497.tsv'
meta_path = '/directflow/SCCGGroupShare/projects/lacgra/bulk.data/metadata.txt'
# exp_path = sys.argv[2]
exp = pd.read_csv(exp_path, sep='\t', index_col = 0)
bulk = os.path.basename(exp_path).replace('.tsv', '')
# meta_path = sys.argv[3]
meta = pd.read_csv(meta_path, sep='\t')
meta['Condition'] = meta['Condition'].replace({"Healthy": 0, "SLE": 1})
meta = meta.loc[meta['GSE'] == bulk,]
meta.index = meta['Sample']
# Match samples in metadata to samples in expression data
exp = exp.reindex(meta.index)
exp['class'] = meta['Condition'].values

# Subset exp rows by meta['Gender'] == 'Female'
if('Female' in meta['Gender']):
    exp = exp.loc[meta['Gender'] == 'Female',]
    # Predict sex based on XIST and RPS4Y1 then filter out males

# # Identify females in the data
# # Scale the data
# exp_scaled = (exp[['XIST', 'RPS4Y1']] - exp[['XIST', 'RPS4Y1']].mean()) / exp[['XIST', 'RPS4Y1']].std()
# # Calculate the dissimilarity matrix
# dissimilarity = pdist(exp_scaled.values, metric='euclidean')
# # Perform hierarchical clustering
# Z = linkage(dissimilarity, method='median')
# # Cut the dendrogram to obtain two clusters
# cluster_result = fcluster(Z, 2, criterion='maxclust')
# # Calculate the mean XIST expression for each cluster
# xist_1 = exp['XIST'][cluster_result == 1].mean()
# xist_2 = exp['XIST'][cluster_result == 2].mean()
# # Assign sex based on XIST expression
# sex_list = []
# if xist_1 > xist_2:
#     sex_list = np.where(cluster_result == 1, 'F', 'M')
# else:
#     sex_list = np.where(cluster_result == 1, 'M', 'F')
# sex_df = pd.DataFrame({'sample':exp.index, 'sex': sex_list})

# exp = exp.loc[exp.index.isin(sex_df.loc[sex_df['sex'] == 'F', 'sample'])]

# Convert features to a pandas Series
features_series = pd.Series(features)
# Check which features are in exp.columns[:-1]
exp[features[~features_series.isin(exp.columns[:-1])]] = 0

# Create datasets
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