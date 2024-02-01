import sys
import os
import pickle
import pandas as pd
import numpy as np
import pyreadr
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score
import matplotlib.pyplot as plt
import shap

file = sys.argv[1]

cell = file.replace('psuedobulk/', '').replace('.RDS', '')

# Read in tune, train and test
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

# load the model from disk
# RF = pickle.load(open('psuedobulk/ML.models/RF_model_'+cell+'.sav', 'rb'))
GBM = pickle.load(open('psuedobulk/ML.models/GBM_model_'+cell+'.sav', 'rb'))

# GBM SHAP values
explainer = shap.Explainer(GBM)
explanation = explainer(X_test.loc[:, features])

shap.plots.beeswarm(explanation, max_display=len(features))
plt.savefig('psuedobulk/SHAP/GBM_'+cell+'.beeswarm.pdf', bbox_inches='tight')
plt.close()

shap.plots.heatmap(explanation, max_display=len(features), instance_order=explanation.sum(1))
plt.savefig('psuedobulk/SHAP/GBM_'+cell+'.heatmap.pdf', bbox_inches='tight')
plt.close()

shap.plots.bar(explanation, max_display=len(features))
plt.savefig('psuedobulk/SHAP/GBM_'+cell+'.barplot.pdf', bbox_inches='tight')
plt.close()

# # RF SHAP values
# explainer = shap.Explainer(RF)
# explanation = explainer(X_test.loc[:, features])

# shap.plots.beeswarm(explanation[:,:,1], max_display=len(features))
# plt.savefig('psuedobulk/SHAP/RF_'+cell+'.beeswarm.pdf', bbox_inches='tight')
# plt.close()

# shap.plots.heatmap(explanation[:,:,1], max_display=len(features))
# plt.savefig('psuedobulk/SHAP/RF_'+cell+'.heatmap.pdf', bbox_inches='tight')
# plt.close()

# shap.plots.bar(explanation[:,:,1], max_display=len(features))
# plt.savefig('psuedobulk/SHAP/RF_'+cell+'.barplot.pdf', bbox_inches='tight')
# plt.close()