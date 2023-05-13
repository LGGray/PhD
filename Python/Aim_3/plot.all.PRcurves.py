import sys
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
# # get basename of file
# cell=os.path.basename(file).replace('.RDS','')

# load the model from disk
logit = pickle.load(open('ML.models/logit_model_'+cell+'.sav', 'rb'))
RF = pickle.load(open('ML.models/RF_model_'+cell+'.sav', 'rb'))
SVM = pickle.load(open('ML.models/SVM_model_'+cell+'.sav', 'rb'))
GBM = pickle.load(open('ML.models/GBM_model_'+cell+'.sav', 'rb'))
MLP = pickle.load(open('ML.models/MLP_model_'+cell+'.sav', 'rb'))

# # Identify model features missing from data
# logit_features = logit.feature_names_in_
# df.columns = df.columns.astype(str)
# missing = np.setdiff1d(logit_features, df.columns)
# # Add missing features to data
# df[missing] = 0
# # fill missing columns by k-nearest neighbors imputation
# from sklearn.impute import KNNImputer
# imputer = KNNImputer(n_neighbors=5)
# df.iloc[:,2::] = pd.DataFrame(imputer.fit_transform(df.iloc[:,2::]), columns = df.iloc[:,2::].columns)

# Replace classes with binary label
df['class'] = df['class'].replace({"control": 0, "disease": 1})

# Collect individual IDs
individuals = df['individual'].unique()
n_individuals = len(individuals)

# Get the number of individuals in each condition
individual_class = df['individual'].astype(str) + '_' + df['class'].astype(str)
n_control = len(individual_class[individual_class.str.endswith('_0')].unique())
n_disease = len(individual_class[individual_class.str.endswith('_1')].unique())

# Calculate the number of individuals to assign to each dataset
n_test = int(n_individuals * 0.2)
n_tune = int(n_individuals * 0.2)
n_train = int(n_individuals * 0.6)

# Randomly assign individuals to each dataset
np.random.seed(42)
test_individuals = np.random.choice(individuals, size=n_test, replace=False)
individuals = np.setdiff1d(individuals, test_individuals)
tune_individuals = np.random.choice(individuals, size=n_tune, replace=False)
individuals = np.setdiff1d(individuals, tune_individuals)
train_individuals = individuals

# Get the corresponding cells for each dataset
test_index = df['individual'].isin(test_individuals)
tune_index = df['individual'].isin(tune_individuals)
train_index = df['individual'].isin(train_individuals)

# Split the data into training, tuning, and testing sets
X_train, X_test, X_tune = df.loc[train_index,].drop(['class','individual'], axis=1), df.loc[test_index,].drop(['class','individual'], axis=1), df.loc[tune_index,].drop(['class','individual'], axis=1)
y_train, y_test, y_tune = df.loc[train_index,'class'], df.loc[test_index,'class'], df.loc[tune_index,'class']

# Get the predicted probabilities
RF_pred_proba = RF.predict_proba(X_test.loc[:, RF.feature_names_in_])[:, 1]
GBM_pred_proba = GBM.predict_proba(X_test.loc[:, GBM.feature_names_in_])[:, 1]
MLP_pred_proba = MLP.predict_proba(X_test.loc[:, MLP.feature_names_in_])[:, 1]
# Scale data. Required for logit to speed up analysis with 'saga' solver
scaler = StandardScaler()
X_test_logit = pd.DataFrame(scaler.fit_transform(X_test), columns=X_test.columns)
logit_pred_proba = logit.predict_proba(X_test_logit.loc[:, logit.feature_names_in_])[:, 1]
# Scale data to have min of 0 and max of 1. Required for SVM
scaler = MinMaxScaler()
X_test_SVM = pd.DataFrame(scaler.fit_transform(X_test), columns=X_test.columns)
SVM_pred_proba = SVM.predict_proba(X_test_SVM.loc[:, SVM.feature_names_in_])[:, 1]

# Calculate F1 scores for each model
logit_f1 = f1_score(y_test, logit.predict(X_test_logit.loc[:, logit.feature_names_in_]))
RF_f1 = f1_score(y_test, RF.predict(X_test.loc[:, RF.feature_names_in_]))
SVM_f1 = f1_score(y_test, SVM.predict(X_test_SVM.loc[:, SVM.feature_names_in_]))
GBM_f1 = f1_score(y_test, GBM.predict(X_test.loc[:, GBM.feature_names_in_]))
MLP_f1 = f1_score(y_test, MLP.predict(X_test.loc[:, MLP.feature_names_in_]))

# f1_scores = {'Logistic':logit_f1, 'RF':RF_f1, 'SVM':SVM_f1, 'GBM':GBM_f1, 'MLP':MLP_f1}

# Calulate precicion and recall values for each model
logit_precision, logit_recall, _ = precision_recall_curve(y_test, logit_pred_proba)
RF_precision, RF_recall, _ = precision_recall_curve(y_test, RF_pred_proba)
SVM_precision, SVM_recall, _ = precision_recall_curve(y_test, SVM_pred_proba)
GBM_precision, GBM_recall, _ = precision_recall_curve(y_test, GBM_pred_proba)
MLP_precision, MLP_recall, _ = precision_recall_curve(y_test, MLP_pred_proba)

# Plot precision-recall curves for all models
# Plot precision-recall curves for all models with F1 score > 0.8
if logit_f1 > 0.7:
    plt.plot(logit_recall, logit_precision, label='Logistic')
if RF_f1 > 0.7:
    plt.plot(RF_recall, RF_precision, label='Random Forest')
if SVM_f1 > 0.7:
    plt.plot(SVM_recall, SVM_precision, label='Support Vector Machine')
if GBM_f1 > 0.7:
    plt.plot(GBM_recall, GBM_precision, label='Gradient Boosting Machine')
if MLP_f1 > 0.7:
    plt.plot(MLP_recall, MLP_precision, label='Multi-Layer Perceptron')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve: ' + cell.replace('.', ' '))
plt.legend()
plt.savefig('precision_recall_curve_'+ cell + '.png', dpi=300)
# plt.savefig('../precision_recall_curve_SLE-pSS_'+ cell +'.png', dpi=300)
plt.show()
plt.close()