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

cell = file.replace('exp.matrix/', '').replace('.RDS', '')

# load the model from disk
logit = pickle.load(open('ML.models/logit_model_'+cell+'.sav', 'rb'))
RF = pickle.load(open('ML.models/RF_model_'+cell+'.sav', 'rb'))
SVM = pickle.load(open('ML.models/SVM_model_'+cell+'.sav', 'rb'))
GBM = pickle.load(open('ML.models/GBM_model_'+cell+'.sav', 'rb'))
MLP = pickle.load(open('ML.models/MLP_model_'+cell+'.sav', 'rb'))

# Replace classes with binary label
df['class'] = df['class'].replace({"control": 0, "disease": 1})

# Read in tune, train, test and features
X_tune = pd.read_csv('data.splits/X_tune.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_tune = pd.read_csv('data.splits/y_tune.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
X_train = pd.read_csv('data.splits/X_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_train = pd.read_csv('data.splits/y_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
X_test = pd.read_csv('data.splits/X_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_test = pd.read_csv('data.splits/y_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
features = pd.read_csv('data.splits/'+os.path.basename(file).replace('.RDS', '')+'.csv')

# Get the predicted probabilities
logit_pred_proba = logit.predict_proba(X_test.loc[:, features.iloc[:,0]])[:, 1]
RF_pred_proba = RF.predict_proba(X_test.loc[:, features.iloc[:,0]])[:, 1]
SVM_pred_proba = SVM.predict_proba(X_test.loc[:, features.iloc[:,0]])[:, 1]
GBM_pred_proba = GBM.predict_proba(X_test.loc[:, features.iloc[:,0]])[:, 1]
MLP_pred_proba = MLP.predict_proba(X_test.loc[:, features.iloc[:,0]])[:, 1]

# # Calculate F1 scores for each model
# logit_f1 = f1_score(y_test['class'], logit.predict(X_test.loc[:, features.iloc[:,0]]))
# RF_f1 = f1_score(y_test['class'], RF.predict(X_test.loc[:, features.iloc[:,0]]))
# SVM_f1 = f1_score(y_test['class'], SVM.predict(X_test.loc[:, features.iloc[:,0]]))
# GBM_f1 = f1_score(y_test['class'], GBM.predict(X_test.loc[:, features.iloc[:,0]]))
# MLP_f1 = f1_score(y_test['class'], MLP.predict(X_test.loc[:, features.iloc[:,0]]))

# f1_scores = {'Logistic':logit_f1, 'RF':RF_f1, 'SVM':SVM_f1, 'GBM':GBM_f1, 'MLP':MLP_f1}

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
plt.savefig('ML.plots/PRCurve_'+ cell +'.png', dpi=300)
plt.show()
plt.close()