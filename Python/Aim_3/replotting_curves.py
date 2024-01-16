import sys
import os
import pickle
import pandas as pd
import numpy as np
import pyreadr
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score
import matplotlib.pyplot as plt
import pickle

file = sys.argv[1]

cell = file.replace('psuedobulk/', '').replace('.RDS', '')

# load the model from disk
logit = pickle.load(open('psuedobulk/ML.models/logit_model_'+cell+'.sav', 'rb'))
RF = pickle.load(open('psuedobulk/ML.models/RF_model_'+cell+'.sav', 'rb'))
SVM = pickle.load(open('psuedobulk/ML.models/SVM_model_'+cell+'.sav', 'rb'))
GBM = pickle.load(open('psuedobulk/ML.models/GBM_model_'+cell+'.sav', 'rb'))
MLP = pickle.load(open('psuedobulk/ML.models/MLP_model_'+cell+'.sav', 'rb'))

# Read in tune, train, test and features
X_train = pd.read_csv('psuedobulk/data.splits/X_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_train = pd.read_csv('psuedobulk/data.splits/y_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
X_test = pd.read_csv('psuedobulk/data.splits/X_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_test = pd.read_csv('psuedobulk/data.splits/y_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
enet_features = pd.read_csv('psuedobulk/features/enet_features.'+os.path.basename(file).replace('.RDS', '')+'.csv')
boruta_features = pd.read_csv('psuedobulk/features/boruta_features.'+os.path.basename(file).replace('.RDS', '')+'.csv')

# Subset for best and tentitive features selected by boruta
boruta_features = boruta_features[boruta_features['Rank'] <= 2]
# Subset elastic net features to those with absolute value of coefficients in 90th percentile
enet_features = enet_features[enet_features['coef'].abs() >= enet_features['coef'].abs().quantile(0.9)]

# Intersection of features selected by Boruta and Elastic Net
features = pd.merge(enet_features, boruta_features, on='Feature', how='inner')['Feature']

# Write features to file
features.to_csv('psuedobulk/features/combined_features.'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

name_mapping = {
    'LogisticRegression': 'logit',
    'RandomForestClassifier': 'RF',
    'SVC': 'SVM',
    'GradientBoostingClassifier': 'GBM',
    'MLPClassifier': 'MLP'
}

for clf in [logit, RF, SVM, GBM, MLP]:
    # Predict the test set
    y_pred = clf.predict(X_test.loc[:, features])
    y_pred_proba = clf.predict_proba(X_test.loc[:, features])[:, 1]

    # Calculate Youden's J statistic to find the optimal threshold
    fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
    j_scores = tpr - fpr
    # Find the optimal threshold
    optimal_idx = np.argmax(j_scores)
    optimal_threshold = thresholds[optimal_idx]
    # Convert probabilities to binary predictions based on optimal threshold
    y_pred = (y_pred_proba >= optimal_threshold).astype(int)

    # Calculate auroc
    auroc = roc_auc_score(y_test, y_pred)

    # Print the AUROC curve
    fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
    plt.plot(fpr, tpr, label='AUC-ROC (area = %0.2f)' % auroc)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(name_mapping[type(clf).__name__] + ': ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
    plt.legend(loc="lower right")
    plt.savefig('psuedobulk/AUROC/'+name_mapping[type(clf).__name__]+'_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')
    plt.close()

    # Print the PR curve
    precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
    average_precision = average_precision_score(y_test, y_pred_proba)
    disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
    disp.plot()
    disp.ax_.set_title(name_mapping[type(clf).__name__] + ': ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
    plt.savefig('psuedobulk/PRC/'+name_mapping[type(clf).__name__]+'_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')
    plt.close()
