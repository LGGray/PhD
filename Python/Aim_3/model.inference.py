import pickle as pk
import pyreadr
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score

# Read in best model file
file = sys.argv[1]
best_model = pd.read_csv(file, sep='\t')

study1 = basename(file).split('.')[0]
study2 = file.split('.')[1]

model_dict = {'LR': LogisticRegression(), 'RF': RandomForestClassifier(), 'GBM': GradientBoostingClassifier(), 'SVM': SVC()}

# Iterate through each model in best_model to evaluate new data
for i in range(0, len(best_model)):
    ML_model = best_model.iloc[i,0].split('_')[0] + '_model_'
    cell_type = best_model.iloc[i,0].split('_')[1]
    # Load data
    df = pyreadr.read_r(study2 + '/exp.matrix/' + cell_type + '.RDS')
    df = df[None]
    # print(df.head())
    df['class'] = df['class'].replace({"control": 0, "disease": 1})

    model = pk.load(open(study1 + '/ML.models/' + ML_model + cell_type + '.sav', 'rb'))

    features = model.feature_names_in_.tolist()
    # # remove second element from list
    # features.pop(1)
    # check if all features in model.columns
    if all(elem in df.columns.tolist() for elem in features):
        print('All features in model')
        Y = df['class']
        X = df[features]
        y_pred = model.predict(X)

        f1 = f1_score(Y, y_pred)
        auroc = roc_auc_score(Y, y_pred)

        # print('F1 score: ', f1)
        # print('AUC: ', auroc)
    else:
        print('Not all features in model')
        

for i in range(0, len(best_model)):
    ML_model = best_model.iloc[i,0].split('_')[0] + '_model_'
    cell_type = best_model.iloc[i,0].split('_')[1]
    # Load data
    df = pyreadr.read_r(study2 + '/exp.matrix/' + cell_type + '.RDS')
    df = df[None]
    # print(df.head())
    df['class'] = df['class'].replace({"control": 0, "disease": 1})

    model = pk.load(open(study1 + '/ML.models/' + ML_model + cell_type + '.sav', 'rb'))

    features = model.feature_names_in_.tolist()
    Y = df['class']
    X = df[features]
    y_pred = model.predict(X)

    f1 = f1_score(Y, y_pred)
    auroc = roc_auc_score(Y, y_pred)

    print('F1 score: ', f1)
    print('AUC: ', auroc)


else:
    