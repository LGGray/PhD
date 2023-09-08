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

# Read in tune, train, test and features
X_tune = pd.read_csv('data.splits/X_tune.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_tune = pd.read_csv('data.splits/y_tune.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
X_train = pd.read_csv('data.splits/X_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_train = pd.read_csv('data.splits/y_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
X_test = pd.read_csv('data.splits/X_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_test = pd.read_csv('data.splits/y_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
features = pd.read_csv('data.splits/'+os.path.basename(file).replace('.RDS', '')+'.csv')

# Create a voting classifier
voting_clf = VotingClassifier(estimators=[('logit', logit), ('RF', RF), ('SVM', SVM), ('GBM', GBM), ('MLP', MLP)], voting='soft')

# Fit the voting classifier
voting_clf.fit(X_train.loc[:, features.iloc[:,0]], y_train['class'])

# Predict the test set
y_pred = voting_clf.predict(X_test.loc[:, features.iloc[:,0]])
# Get the predicted probabilities
y_pred_proba = voting_clf.predict_proba(X_test.loc[:, features.iloc[:,0]])[:, 1]

# Calculate the metrics
accuracy = accuracy_score(y_test, y_pred)
precision = precision_score(y_test, y_pred)
recall = recall_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred)
auc = roc_auc_score(y_test, y_pred)
kappa = cohen_kappa_score(y_test, y_pred)

# Create dataframe of metrics and save to file
metrics = pd.DataFrame({'Accuracy': [accuracy], 
                        'Precision': [precision], 
                        'Recall': [recall], 
                        'F1': [f1], 
                        'AUC': [auc],
                        'Kappa': [kappa]})

metrics.to_csv('ML.models/ensemble/metrics_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Save confusion matrix to file
confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
confusion.to_csv('ML.models/ensemble/confusion_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Print the PR curve
precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
average_precision = average_precision_score(y_test, y_pred_proba)
disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
disp.plot()
disp.ax_.set_title('ensemble: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.savefig('ML.models/ensemble/PRcurve_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Save the model
import pickle
filename = 'ML.models/ensemble/'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(voting_clf, open(filename, 'wb'))
