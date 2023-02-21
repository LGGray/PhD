import pickle
import pyreadr
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score

# Load in new data
test = pyreadr.read_r('exp.matrix/Plasma.cells.deg.RDS')
test = test[None]
print(test.head())
test['class'] = test['class'].replace({"control": 0, "disease": 1})

# Load the model
with open('../UC_GSE125527/ML.models/RF_model_Plasma.cells.deg.sav', 'rb') as f:
    trained_model = pickle.load(f)

features = trained_model.get_feature_names_out().tolist()
# Subset the data to the features used in the model even if missing
Y = test['class']
test = test.loc[:, test.columns.isin(features)]
features = test.columns.tolist()
# Scale data to match the model
scaler = StandardScaler()
test = pd.DataFrame(scaler.fit_transform(test))
test.columns = features

clf = RandomForestClassifier(n_estimators=100, 
                            max_depth=30, 
                            min_samples_split=2, n_jobs=-1)
clf = clf.fit(test, Y)
y_pred = clf.predict(test)

accuracy = accuracy_score(Y, y_pred)
precision = precision_score(Y, y_pred)
recall = recall_score(Y, y_pred)
f1 = f1_score(Y, y_pred)
auc = roc_auc_score(Y, y_pred)
confusion_matrix(Y, y_pred)

# tabulate Y values
Y.value_counts()
