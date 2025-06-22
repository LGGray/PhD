import pickle
import pandas as pd
import sys
import os
import numpy as np
from collections import Counter
from sklearn.ensemble import VotingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score, cohen_kappa_score, matthews_corrcoef
import matplotlib.pyplot as plt
import pyreadr
from sklearn.preprocessing import MinMaxScaler, StandardScaler

file = sys.argv[1]

cell = file.replace('pseudobulk/', '').replace('.RDS', '')

### Read in independent test set
test = pyreadr.read_r(f'/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/SLE_GSE135779/pseudobulk_update/{cell}.RDS')
test = test[None]

### Read in all features file
all_features = pd.read_csv('figures.chrX_vs_SLE/top_celltypes_all_features.csv')
features = all_features[all_features['celltype'] == cell]['feature'].tolist()

# Which features are in test.columns
features = [feature for feature in features if feature in test.columns]

# Read in ensembe models
ensemble_models = []
for i in range(1, 11):
    model = pickle.load(open(f'pseudobulk/split_{i}/combined/ensemble/{cell}.sav', 'rb'))
    ensemble_models.append(model)

# Function to calculate the average or mode of hyperparameters
def average_hyperparams(models, param_path):
    values = [eval(f'model.{param_path}') for model in models]
    if isinstance(values[0], (int, float)):
        return int(np.mean(values))
    else:
        return Counter(values).most_common(1)[0][0]
    
# Averaging or taking the mode of hyperparameters
if average_hyperparams(ensemble_models, 'estimators[0][1].C') == 0:
    average_logit_C = 0.1
else:
    average_logit_C = average_hyperparams(ensemble_models, 'estimators[0][1].C')

average_logit_penalty = average_hyperparams(ensemble_models, 'estimators[0][1].penalty')
average_logit_solver = average_hyperparams(ensemble_models, 'estimators[0][1].solver')
average_logit_l1_ratio = average_hyperparams(ensemble_models, 'estimators[0][1].l1_ratio')

average_rf_max_depth = average_hyperparams(ensemble_models, 'estimators[1][1].max_depth')
average_rf_min_samples_split = average_hyperparams(ensemble_models, 'estimators[1][1].min_samples_split')
average_rf_n_estimators = average_hyperparams(ensemble_models, 'estimators[1][1].n_estimators')

if average_hyperparams(ensemble_models, 'estimators[2][1].C') == 0:
    average_svm_C = 0.1
else:
    average_svm_C = average_hyperparams(ensemble_models, 'estimators[2][1].C')

#average_gbm_max_features = average_hyperparams(ensemble_models, 'estimators[3][1].max_features')

average_mlp_alpha = average_hyperparams(ensemble_models, 'estimators[4][1].alpha')
average_mlp_max_iter = average_hyperparams(ensemble_models, 'estimators[4][1].max_iter')

# Creating the new model with the average/mode hyperparameters
ensemble = VotingClassifier(
    estimators=[
        ('logit', LogisticRegression(C=0.1, penalty=average_logit_penalty, solver=average_logit_solver, l1_ratio=average_logit_l1_ratio, class_weight='balanced', max_iter=10000, n_jobs=8, random_state=42)),
        ('RF', RandomForestClassifier(max_depth=average_rf_max_depth, min_samples_split=average_rf_min_samples_split, n_estimators=average_rf_n_estimators, class_weight='balanced', n_jobs=8)),
        ('SVM', SVC(C=0.1, class_weight='balanced', probability=True, random_state=42)),
        ('GBM', GradientBoostingClassifier(max_features='sqrt', random_state=42, subsample=1)),
        ('MLP', MLPClassifier(alpha=average_mlp_alpha, early_stopping=True, max_iter=average_mlp_max_iter, random_state=42))
    ],
    voting='soft'
)

# Retrain the ensemble model on the training set
X_train = pd.read_csv(f'pseudobulk/split_1/data.splits/X_train.'+cell+'.csv', index_col=0)
y_train = pd.read_csv(f'pseudobulk/split_1/data.splits/y_train.'+cell+'.csv', index_col=0)
X_test = pd.read_csv(f'pseudobulk/split_1/data.splits/X_test.'+cell+'.csv', index_col=0)
y_test = pd.read_csv(f'pseudobulk/split_1/data.splits/y_test.'+cell+'.csv', index_col=0)

# Fit the ensemble model
ensemble.fit(X_train.loc[:, features], y_train['class'])

# Predict the test set
y_pred = ensemble.predict(X_test.loc[:, features])
y_pred_proba = ensemble.predict_proba(X_test.loc[:, features])[:, 1]

# Calculate Youden's J statistic to find the optimal threshold
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
j_scores = tpr - fpr
# Find the optimal threshold
optimal_idx = np.argmax(j_scores)
optimal_threshold = thresholds[optimal_idx]
# Convert probabilities to binary predictions based on optimal threshold
y_pred = (y_pred_proba >= optimal_threshold).astype(int)

# calculate the MCC
mcc = matthews_corrcoef(y_test['class'], y_pred)

# Subset for indepdnent data for adult, child or all
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

test = test.drop(columns=['class', 'individual'])

# Scale the data
scaler = StandardScaler()
X_test = pd.DataFrame(scaler.fit_transform(test), columns=test.columns, index=test.index)

# Predict the test set
y_pred = ensemble.predict(X_test.loc[:, features])
y_pred_proba = ensemble.predict_proba(X_test.loc[:, features])[:, 1]

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
metrics.to_csv(f'/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/SLE_GSE135779/ML.plots_update/metrics_{cell}.chrX.csv', index=False)
