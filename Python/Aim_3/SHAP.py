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

cell = file.replace('pseudobulk/', '').replace('.RDS', '')

# Read in tune, train, test and features
X_train = pd.read_csv('pseudobulk/data.splits/X_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_train = pd.read_csv('pseudobulk/data.splits/y_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
X_test = pd.read_csv('pseudobulk/data.splits/X_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_test = pd.read_csv('pseudobulk/data.splits/y_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
enet_features = pd.read_csv('pseudobulk/features/enet_features.'+os.path.basename(file).replace('.RDS', '')+'.csv')
boruta_features = pd.read_csv('pseudobulk/features/boruta_features.'+os.path.basename(file).replace('.RDS', '')+'.csv')

# Subset for best and tentitive features selected by boruta
boruta_features = boruta_features[boruta_features['Rank'] == 1]
# Subset elastic net features to those with absolute value of coefficients in 90th percentile
threshold = np.percentile(np.abs(enet_features['coef']), 90)
enet_features = enet_features[enet_features['coef'] > threshold]

#### Condition for command-line argument indicating feature type ###
if sys.argv[2] == 'intersection':
    # Intersection of features selected by Boruta and Elastic Net
    features = pd.merge(enet_features, boruta_features, on='Feature', how='inner')['Feature']
    if(len(features) == 0):
        print("No common features between Boruta and Elastic Net")
        sys.exit()
elif sys.argv[2] == 'combined':
    # Features selected by Boruta and Elastic Net
    features = pd.merge(enet_features, boruta_features, on='Feature', how='outer')['Feature']
elif sys.argv[2] == 'boruta':
    features = boruta_features['Feature']
elif sys.argv[2] == 'enet':
    features = enet_features['Feature']

# load the model from disk
ensemble = pickle.load(open('pseudobulk/'+sys.argv[2]+'/ML.models/ensemble/'+cell+'.sav', 'rb'))

# Create a wrapper function for the predict_proba method of VotingClassifier
def voting_classifier_proba(data):
    return ensemble.predict_proba(data)

explainer = shap.KernelExplainer(voting_classifier_proba, X_train.loc[:, features])

# shap_values = explainer.shap_values(X_test.loc[:, features])

explanation = explainer(X_test.loc[:, features])

shap_values_single_class = explanation[..., 1]  # Adjust index based on the class you are interested in
shap.plots.beeswarm(shap_values_single_class, max_display=len(features))
plt.savefig('pseudobulk/'+sys.argv[2]+'/ML.models/ensemble/'+cell+'.beeswarm.pdf', bbox_inches='tight')
plt.close()

# shap.plots.heatmap(explanation, max_display=len(features), instance_order=explanation.sum(1))
# plt.savefig('psuedobulk/SHAP/GBM_'+cell+'.heatmap.pdf', bbox_inches='tight')
# plt.close()

# shap.plots.bar(explanation, max_display=len(features))
# plt.savefig('psuedobulk/SHAP/GBM_'+cell+'.barplot.pdf', bbox_inches='tight')
# plt.close()

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