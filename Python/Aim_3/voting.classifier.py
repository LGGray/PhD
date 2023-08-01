import sys
import pickle
import os.path
import pandas as pd
import numpy as np
import pyreadr
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score
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

# create a list of all features and remove duplicates
features = []
for model in models:
    features.extend(model.feature_names_in_)
features = list(set(features))

# Replace classes with binary label
df['class'] = df['class'].replace({"control": 0, "disease": 1})

# Stratified split the data into train (80%) and test (20%) sets
X = df.drop(['class','individual'], axis=1)
y = df['class']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=df['class'])

# Scale the data
scaler = StandardScaler()
X_train = pd.DataFrame(scaler.fit_transform(X_train), columns=X_train.columns)
X_test = pd.DataFrame(scaler.fit_transform(X_test), columns=X_test.columns)

# Create a voting classifier
voting_clf = VotingClassifier(estimators=[('logit', logit), ('RF', RF), ('SVM', SVM), ('GBM', GBM), ('MLP', MLP)], voting='soft')

# Fit the voting classifier
voting_clf.fit(X_train.loc[:, features], y_train)

# Predict the test set
y_pred = voting_clf.predict(X_test.loc[:, features])
# Get the predicted probabilities
y_pred_proba = voting_clf.predict_proba(X_test.loc[:, features])[:, 1]

# Calculate the metrics
accuracy = accuracy_score(y_test, y_pred)
precision = precision_score(y_test, y_pred)
recall = recall_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred)
auc = roc_auc_score(y_test, y_pred)

# Create dataframe of metrics and save to file
metrics = pd.DataFrame({'Accuracy': [accuracy], 
                        'Precision': [precision], 
                        'Recall': [recall], 
                        'F1': [f1], 
                        'AUC': [auc]})

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
