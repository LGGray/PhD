import sys
import os
import pickle
import pandas as pd
import numpy as np
import pyreadr
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score
import matplotlib.pyplot as plt

file = sys.argv[1]

# Read in expression RDS file
df = pyreadr.read_r(file)
df = df[None]
print(df.head())

cell = file.replace('psuedobulk/', '').replace('.RDS', '')

# load the model from disk
logit = pickle.load(open('psuedobulk/ML.models/logit_model_'+cell+'.sav', 'rb'))
RF = pickle.load(open('psuedobulk/ML.models/RF_model_'+cell+'.sav', 'rb'))
SVM = pickle.load(open('psuedobulk/ML.models/SVM_model_'+cell+'.sav', 'rb'))
GBM = pickle.load(open('psuedobulk/ML.models/GBM_model_'+cell+'.sav', 'rb'))
MLP = pickle.load(open('psuedobulk/ML.models/MLP_model_'+cell+'.sav', 'rb'))

# Replace classes with binary label
df['class'] = df['class'].replace({"control": 0, "disease": 1})

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

# Get the predicted probabilities
logit_pred_proba = logit.predict_proba(X_test.loc[:, features])[:, 1]
RF_pred_proba = RF.predict_proba(X_test.loc[:, features])[:, 1]
SVM_pred_proba = SVM.predict_proba(X_test.loc[:, features])[:, 1]
GBM_pred_proba = GBM.predict_proba(X_test.loc[:, features])[:, 1]
MLP_pred_proba = MLP.predict_proba(X_test.loc[:, features])[:, 1]

# Create a dictionary of model names and predicted probabilities
models = {'logit': logit_pred_proba,
          'RF': RF_pred_proba,
          'SVM': SVM_pred_proba,
          'GBM': GBM_pred_proba,
          'MLP': MLP_pred_proba}
# Loop through the dictionary and plot models
for name, proba in models.items():
    # Calculate the precision and recall for each model
    precision, recall, _ = precision_recall_curve(y_test['class'], proba)
    # Plot the curve and add the label
    plt.plot(recall, precision, label=name)
# Add labels and legend
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve: ' + cell.replace('.', ' '))
plt.legend()
plt.savefig('psuedobulk/ML.plots/PRCurve_'+ cell +'.pdf', dpi=300)
plt.show()
plt.close()