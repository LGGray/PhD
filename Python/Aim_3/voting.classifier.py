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

file = sys.argv[1]

# Read in expression RDS file
df = pyreadr.read_r(file)
df = df[None]
print(df.head())

cell = file.replace('new_pseudobulk/', '').replace('.RDS', '')

# load the model from disk
logit = pickle.load(open(f'new_pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ML.models/logit_model_'+cell+'.sav', 'rb'))
RF = pickle.load(open(f'new_pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ML.models/RF_model_'+cell+'.sav', 'rb'))
SVM = pickle.load(open(f'new_pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ML.models/SVM_model_'+cell+'.sav', 'rb'))
GBM = pickle.load(open(f'new_pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ML.models/GBM_model_'+cell+'.sav', 'rb'))
MLP = pickle.load(open(f'new_pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ML.models/MLP_model_'+cell+'.sav', 'rb'))

# Read in tune, train, test and features
X_train = pd.read_csv(f'new_pseudobulk/split_{sys.argv[2]}/data.splits/X_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_train = pd.read_csv(f'new_pseudobulk/split_{sys.argv[2]}/data.splits/y_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
X_test = pd.read_csv(f'new_pseudobulk/split_{sys.argv[2]}/data.splits/X_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_test = pd.read_csv(f'new_pseudobulk/split_{sys.argv[2]}/data.splits/y_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
enet_features = pd.read_csv(f'new_pseudobulk/split_{sys.argv[2]}/features/enet_features.'+os.path.basename(file).replace('.RDS', '')+'.csv')
boruta_features = pd.read_csv(f'new_pseudobulk/split_{sys.argv[2]}/features/boruta_features.'+os.path.basename(file).replace('.RDS', '')+'.csv')

# Subset for selected and tentitive features from boruta
boruta_features = boruta_features[boruta_features['Rank'] == 1]
# Subset elastic net features to those with absolute value of coefficients in 80th percentile
threshold = np.percentile(np.abs(enet_features['coef']), 90)
enet_features = enet_features[np.abs(enet_features['coef']) >= threshold]

#### Condition for command-line argument indicating feature type ###
if sys.argv[3] == 'intersection':
    # Intersection of features selected by Boruta and Elastic Net
    features = pd.merge(enet_features, boruta_features, on='Feature', how='inner')['Feature']
    if(len(features) == 0):
        print("No common features between Boruta and Elastic Net")
        sys.exit()
elif sys.argv[3] == 'combined':
    # Features selected by Boruta and Elastic Net
    features = pd.merge(enet_features, boruta_features, on='Feature', how='outer')['Feature']
elif sys.argv[3] == 'boruta':
    features = boruta_features['Feature']
elif sys.argv[3] == 'enet':
    features = enet_features['Feature']

# Create a voting classifier
voting_clf = VotingClassifier(estimators=[('logit', logit), ('RF', RF), ('SVM', SVM), ('GBM', GBM), ('MLP', MLP)], voting='soft')

# Fit the voting classifier
voting_clf.fit(X_train.loc[:, features], y_train['class'])

# Predict the test set
y_pred = voting_clf.predict(X_test.loc[:, features])
y_pred_proba = voting_clf.predict_proba(X_test.loc[:, features])[:, 1]

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
                        'Kappa': [kappa]})
metrics.to_csv(f'new_pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ensemble/metrics_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Save confusion matrix to file
confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
confusion.to_csv(f'new_pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ensemble/confusion_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

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
    ax.text(j, i, f'{val}', ha='center', va='center', color='black')
# Set axis labels
ax.set_xlabel('Predicted labels')
ax.set_ylabel('True labels')
ax.set_xticks(range(len(classes)))
ax.set_yticks(range(len(classes)))
ax.set_xticklabels(classes)
ax.set_yticklabels(classes)
# Set the title
ax.set_title('Ensemble: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
# Annotate with F1 score
plt.annotate(f'F1 Score: {f1:.2f}', xy=(0.5, -0.1), xycoords='axes fraction', 
             ha='center', va='center', fontsize=12, color='white')
# Adjust layout for visibility
plt.tight_layout()
# Save the figure
plt.savefig(f'new_pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ensemble/confusion_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')
plt.close()

# Print the PR curve
precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
average_precision = average_precision_score(y_test, y_pred_proba)
disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
disp.plot()
disp.ax_.set_title('Ensemble: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.savefig(f'new_pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ensemble/PRcurve_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Save the model
import pickle
filename = f'new_pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ensemble/'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(voting_clf, open(filename, 'wb'))

# ### Feature permutation importance ###
# result = permutation_importance(voting_clf, X_test.loc[:, features], y_test['class'], n_repeats=30, 
#                            random_state=42, n_jobs=8, scoring=['f1_weighted', 'average_precision'])
# # Print the feature importance for F1
# F1_features = []
# for i in result['f1_weighted']['importances_mean'].argsort()[::-1]:
#     if result['f1_weighted']['importances_mean'][i] - 2 * result['f1_weighted']['importances_std'][i] > 0:
#         F1_features.append(X_test.loc[:, features].columns[i])
#         print(f"{X_test.loc[:, features].columns[i]:<8}"
#               f"{result['f1_weighted']['importances_mean'][i]:.3f}"
#               f" +/- {result['f1_weighted']['importances_std'][i]:.3f}")
        
# # Print the feature importance for AUC
# auprc_features = []
# for i in result['average_precision']['importances_mean'].argsort()[::-1]:
#     if result['average_precision']['importances_mean'][i] - 2 * result['average_precision']['importances_std'][i] > 0:
#         auprc_features.append(X_test.loc[:, features].columns[i])
#         print(f"{X_test.loc[:, features].columns[i]:<8}"
#               f"{result['average_precision']['importances_mean'][i]:.3f}"
#               f" +/- {result['average_precision']['importances_std'][i]:.3f}")

# # Intersection of F1 and AUC features
# important_features = set(F1_features + auprc_features)

# # Retrain model with important features
# voting_clf.fit(X_train.loc[:, important_features], y_train['class'])

# # Predict the test set
# y_pred = voting_clf.predict(X_test.loc[:, important_features])
# # Get the predicted probabilities
# y_pred_proba = voting_clf.predict_proba(X_test.loc[:, important_features])[:, 1]

# # Calculate Youden's J statistic to find the optimal threshold
# fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
# j_scores = tpr - fpr
# # Find the optimal threshold
# optimal_idx = np.argmax(j_scores)
# optimal_threshold = thresholds[optimal_idx]
# # Convert probabilities to binary predictions based on optimal threshold
# y_pred = (y_pred_proba >= optimal_threshold).astype(int)

# # Calculate the metrics
# accuracy = accuracy_score(y_test, y_pred)
# precision = precision_score(y_test, y_pred)
# recall = recall_score(y_test, y_pred)
# f1 = f1_score(y_test, y_pred, average='weighted')
# auroc = roc_auc_score(y_test, y_pred)
# auprc = average_precision_score(y_test, y_pred)
# kappa = cohen_kappa_score(y_test, y_pred)

# # Define bootstrap parameters
# n_bootstraps = 1000
# confidence_level = 0.9
# # Initialize an empty list to store bootstrap scores
# bootstrapped_f1 = []
# bootstrapped_AUC = []
# bootstrapped_AUPRC = []

# # Loop over bootstrap samples
# for i in range(n_bootstraps):
#     # Resample with replacement
#     y_test_resampled, y_pred_resampled = resample(y_test, y_pred, stratify=y_test)
#     # Calculate F1 score
#     f1 = f1_score(y_test_resampled, y_pred_resampled, average='weighted')
#     # Calculate AUROC
#     auroc = roc_auc_score(y_test_resampled, y_pred_resampled)
#     # Calculate AUPRC
#     auprc = average_precision_score(y_test_resampled, y_pred_resampled)
#     # Append score to list
#     bootstrapped_f1.append(f1)
#     bootstrapped_AUC.append(auroc)
#     bootstrapped_AUPRC.append(auprc)

# # Calculate percentile for confidence intervals
# lower_percentile = (1 - confidence_level) / 2 * 100
# upper_percentile = (1 + confidence_level) / 2 * 100

# f1_lower_bound = np.percentile(bootstrapped_f1, lower_percentile)
# f1_upper_bound = np.percentile(bootstrapped_f1, upper_percentile)

# auroc_lower_bound = np.percentile(bootstrapped_AUC, lower_percentile)
# auroc_upper_bound = np.percentile(bootstrapped_AUC, upper_percentile)

# auprc_lower_bound = np.percentile(bootstrapped_AUPRC, lower_percentile)
# auprc_upper_bound = np.percentile(bootstrapped_AUPRC, upper_percentile)

# # Create dataframe of metrics and save to file
# metrics = pd.DataFrame({'Accuracy': [accuracy], 
#                         'Precision': [precision], 
#                         'Recall': [recall], 
#                         'F1': [f1],
#                         'F1_lower': [f1_lower_bound],
#                         'F1_upper': [f1_upper_bound],
#                         'AUC': [auroc],
#                         'AUC_lower': [auroc_lower_bound],
#                         'AUC_upper': [auroc_upper_bound],
#                         'AUPRC': [auprc],
#                         'AUPRC_lower': [auprc_lower_bound],
#                         'AUPRC_upper': [auprc_upper_bound],
#                         'Kappa': [kappa]})
# metrics.to_csv('psuedobulk/ML.models/ensemble/perm_metrics_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# # Save confusion matrix to file
# confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
# confusion.to_csv('psuedobulk/ML.models/ensemble/perm_confusion_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

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
#     ax.text(j, i, f'{val}', ha='center', va='center', color='black')
# # Set axis labels
# ax.set_xlabel('Predicted labels')
# ax.set_ylabel('True labels')
# ax.set_xticks(range(len(classes)))
# ax.set_yticks(range(len(classes)))
# ax.set_xticklabels(classes)
# ax.set_yticklabels(classes)
# # Set the title
# ax.set_title('Ensemble permuted features: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
# # Annotate with F1 score
# plt.annotate(f'F1 Score: {f1:.2f}', xy=(0.5, -0.1), xycoords='axes fraction', 
#              ha='center', va='center', fontsize=12, color='black')
# # Adjust layout for visibility
# plt.tight_layout()
# # Save the figure
# plt.savefig('psuedobulk/ML.models/ensemble/perm_confusion_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')
# plt.close()

# # Print the PR curve
# precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
# average_precision = average_precision_score(y_test, y_pred_proba)
# disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
# disp.plot()
# disp.ax_.set_title('Ensemble permuted features: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
# plt.savefig('psuedobulk/ML.models/ensemble/perm_PRcurve_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# # Save the model
# import pickle
# filename = 'psuedobulk/ML.models/ensemble/perm.'+os.path.basename(file).replace('.RDS', '')+'.sav'
# pickle.dump(voting_clf, open(filename, 'wb'))