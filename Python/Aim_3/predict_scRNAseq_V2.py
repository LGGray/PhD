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

voting_clf = pickle.load(open('psuedobulk/ML.models/ensemble/'+sys.argv[1], 'rb'))

celltype = os.path.basename(sys.argv[1]).replace('perm_', '').replace('.RDS', '')

# save model features in to an vector
features = pd.Series(voting_clf.feature_names_in_)

infile = sys.argv[2]

# Read in expression RDS file
df = pyreadr.read_r(infile)
df = df[None]
print(df.head())

# Replace class labels with binary label
df['class'] = df['class'].replace({"control": 0, "disease": 1})

# Create dummy variable for missing features
missing_features = features[~features.isin(df.columns[1:])]
missing_features = pd.DataFrame(np.zeros((df.shape[0], missing_features.shape[0])), columns=missing_features)
missing_features.index = df.index
df = pd.concat([df, missing_features], axis=1)
# Create X and y
X = df.loc[:, features]
y = df['class']

# Scale data
scaler = StandardScaler()
X = pd.DataFrame(scaler.fit_transform(X), columns=X.columns, index=X.index)

# Predict the data
y_pred = voting_clf.predict(X)
# Get the predicted probabilities
y_pred_proba = voting_clf.predict_proba(X)[:, 1]

# Calculate the metrics
accuracy = accuracy_score(y, y_pred)
precision = precision_score(y, y_pred)
recall = recall_score(y, y_pred)
f1 = f1_score(y, y_pred)
auc = roc_auc_score(y, y_pred)
kappa = cohen_kappa_score(y, y_pred)

metrics = pd.DataFrame({'Accuracy': [accuracy], 
                        'Precision': [precision], 
                        'Recall': [recall], 
                        'F1': [f1],
                        'AUC': [auc],
                        'Kappa': [kappa]})
metrics.to_csv(os.path.dirname(infile)+'/'+celltype+'_metrics.csv', index=False)

print(metrics)

confusion = pd.DataFrame(confusion_matrix(y, y_pred))
confusion.to_csv(os.path.dirname(infile)+'/'+celltype+'_confusion.csv', index=False)

print(confusion)


# #### OneK1K data ncM ####
# >>> shape(df)
# 15341, 21

# >>> metrics
#    Accuracy  Precision    Recall        F1       AUC     Kappa
# 0  0.206122   0.062437  0.843023  0.116264  0.503544  0.001036

# >>> confusion
#       0     1
# 0  1282  6532
# 1    81   435

# >>> np.unique(y_pred, return_counts=True)
# (array([0, 1]), array([1363, 6967]))

# #### pSS data ncM ####
# >>> metrics
#    Accuracy  Precision   Recall        F1       AUC     Kappa
# 0  0.612654   0.646125  0.84908  0.733828  0.530569  0.068622
# >>> confusion = pd.DataFrame(confusion_matrix(y, y_pred))
# >>> confusion
#      0    1
# 0  102  379
# 1  123  692

# #### GSE108497 ncM ####
#    Accuracy  Precision    Recall        F1       AUC     Kappa
# 0  0.675781   0.685315  0.904615  0.779841  0.591345  0.207149
# # confusion
#     0    1
# 0  52  135
# 1  31  294
# ### GSE108497 ABC ####
#    Accuracy  Precision    Recall        F1       AUC     Kappa
# 0  0.623047   0.673684  0.787692  0.726241  0.562295  0.132889
# # confusion
#     0    1
# 0  63  124
# 1  69  256

# ### GSE108497 ncM ####
#    Accuracy  Precision    Recall        F1       AUC     Kappa
# 0  0.710145   0.754098  0.901961  0.821429  0.534314  0.083665
# # confusion
#     0   1
# 0   6  30
# 1  10  92

# ### GSE108497 ABC ####
#    Accuracy  Precision    Recall        F1       AUC     Kappa
# 0  0.746377   0.813084  0.852941  0.832536  0.648693  0.311377
# # confusion
#     0   1
# 0  16  20
# 1  15  87