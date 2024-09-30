import sys
import pickle
import os.path
import re
import pandas as pd
import numpy as np
import pyreadr
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score, cohen_kappa_score, matthews_corrcoef
import matplotlib.pyplot as plt
from sklearn.ensemble import VotingClassifier
from sklearn.inspection import permutation_importance
from sklearn.utils import resample

# Define the substrings to remove
substrings_to_remove = ['.chrX', '.autosome', '.HVG.autosome', '.HVG', '.SLE', '.sav']
# Create a regex pattern to match any of the substrings
pattern = '|'.join(map(re.escape, substrings_to_remove))
# Replace the substrings with an empty string
cell = re.sub(pattern, '', sys.argv[1])

# Read in expression RDS file
df = pyreadr.read_r(f'new_pseudobulk/flare/'+cell+'.RDS')
df = df[None]
print(df.head())

# Remove 'Hispanic or Latin American' and 'African American' from ancestry
df = df[~df['ancestry'].isin(['Hispanic or Latin American', 'African American'])]

# Replace class with binary label
df['class'] = df['class'].replace({"control": 0, "disease": 1})
# Replace ancestry with binary label
df['ancestry'] = df['ancestry'].replace({"European": 0, "Asian": 1})

flare = pd.read_csv(f'new_pseudobulk/split_{sys.argv[2]}/flare_test_index.csv')

test = df[df.index.isin(flare['rownames'])]

# load the model from disk
model = pickle.load(open(f'new_pseudobulk/split_{sys.argv[2]}/intersection/ensemble/'+sys.argv[1], 'rb'))

# Get the features
features = model.feature_names_in_.tolist()

# Get the test set
X_test = test.loc[:, features]
y_test = test['class']

# scale data
scaler = StandardScaler()
X_test = pd.DataFrame(scaler.fit_transform(X_test), columns=X_test.columns, index=X_test.index)

# Predict the test set
y_pred = model.predict(X_test)
y_pred_proba = model.predict_proba(X_test)[:, 1]

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
mcc = matthews_corrcoef(y_test, y_pred)

# Create dataframe of metrics and save to file
metrics = pd.DataFrame({'Accuracy': [accuracy], 
                        'Precision': [precision], 
                        'Recall': [recall], 
                        'F1': [f1],
                        'AUC': [auroc],
                        'AUPRC': [auprc],
                        'Kappa': [kappa],
                        'MCC': [mcc],
                        'n_features': [len(features)]})
metrics.to_csv(f'new_pseudobulk/split_{sys.argv[2]}/flare/metrics_'+sys.argv[2].replace('.sav', '.csv'), index=False)

# Save confusion matrix to file
confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
confusion.to_csv(f'new_pseudobulk/split_{sys.argv[2]}/flare/confusion_'+sys.argv[2].replace('.sav', '.csv'), index=False)

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
    ax.text(j, i, f'{val}', ha='center', va='center', color='white')
# Set axis labels
ax.set_xlabel('Predicted labels')
ax.set_ylabel('True labels')
ax.set_xticks(range(len(classes)))
ax.set_yticks(range(len(classes)))
ax.set_xticklabels(classes)
ax.set_yticklabels(classes)
# Set the title
ax.set_title(sys.argv[1].replace('.sav', '').replace('.', ' '))
# Annotate with F1 score
plt.annotate(f'MCC: {mcc:.2f}', xy=(0.5, -0.1), xycoords='axes fraction', 
             ha='center', va='center', fontsize=12, color='black')
# Adjust layout for visibility
plt.tight_layout()
# Save the figure
plt.savefig(f'new_pseudobulk/split_{sys.argv[2]}/flare/confusion_'+sys.argv[2].replace('.sav', '.pdf'), bbox_inches='tight')
plt.close()

# Print the PR curve
precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
average_precision = average_precision_score(y_test, y_pred_proba)
disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
disp.plot()
disp.ax_.set_title(sys.argv[1].replace('.sav', '').replace('.', ' '))
plt.savefig(f'new_pseudobulk/split_{sys.argv[2]}/flare/PRcurve_'+sys.argv[2].replace('.sav', '.pdf'), bbox_inches='tight')
plt.close()