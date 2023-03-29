import os
import pandas as pd
import pyreadr
import pickle as pk
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler


studies = ['AD_GSE147424', 'MS_GSE193770', 'pSS_GSE157278', 'SLE_SDY997', 'UC_GSE125527', 'UC_GSE182270']
celltype = ['Regulatory.T.cells', 'Tem.Trm.cytotoxic.T.cells', 'Tcm.Naive.helper.T.cells', 'Tem.Effector.helper.T.cells']
models = ['GBM_model_', 'RF_model_', 'SVM_model_', 'logit_model_']

results = []
for i, study in enumerate(studies):
    for cell in celltype:
        for j, model_name in enumerate(models):
            model_path = study+'/ML.models/' + model_name + cell + '.common.sav'
            if os.path.exists(model_path):
                model = pk.load(open(model_path, 'rb'))
                for k, test_study in enumerate(studies):
                    if i != k:
                        df_path = test_study+'/exp.matrix/'+cell+'.common.RDS'
                        if os.path.exists(df_path):
                            df = pyreadr.read_r(df_path)
                            df = df[None]
                            df['class'] = df['class'].replace({"control": 0, "disease": 1})
                            X = df.iloc[:, 1:]
                            y = df.iloc[:, 0]
                            try:
                                X = X[model.feature_names_in_]
                                if model_name == 'logit_model_':
                                    scaler = StandardScaler()
                                    X = pd.DataFrame(scaler.fit_transform(X), columns=X.columns)
                                elif model_name == 'SVM_model_':
                                    scaler = MinMaxScaler()
                                    X = pd.DataFrame(scaler.fit_transform(X), columns=X.columns)
                                else:
                                    X = X
                                y_pred = model.predict(X)
                                f1 = f1_score(y, y_pred)
                                auc = roc_auc_score(y, y_pred)
                                if f1 > 0.7 and auc > 0.6:
                                    result = {'ID':study+'_'+model_name+'_'+cell+'_'+test_study, 'f1': f1, 'auc': auc, 'confusion_matrix': confusion_matrix(y, y_pred).tolist()}
                                    results.append(result)
                            except KeyError as e:
                                print(f"Skipping file {df_path} due to KeyError: {e}")



for result in results:
    print(result)


model = pk.load(open('pSS_GSE157278/ML.models/logit_model_Tem.Trm.cytotoxic.T.cells.common.sav', 'rb'))
df = pyreadr.read_r('SLE_SDY997/exp.matrix/Tem.Trm.cytotoxic.T.cells.common.RDS')
df = df[None]
df['class'] = df['class'].replace({"control": 0, "disease": 1})
X = df.iloc[:, 1:]
y = df.iloc[:, 0]
X = X[model.feature_names_in_]
scaler = StandardScaler()
X = pd.DataFrame(scaler.fit_transform(X), columns=X.columns)
y_pred = model.predict(X)
f1 = f1_score(y, y_pred)
auc = roc_auc_score(y, y_pred)
print(f1, auc)
print(confusion_matrix(y, y_pred))

# Print the PR curve
precision, recall, thresholds = precision_recall_curve(y, y_pred)
average_precision = average_precision_score(y, y_pred)
disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
disp.plot()
disp.ax_.set_title('Trm-CTL.logit.pSS-SLE')
plt.savefig('Trm-CTL.logit.pSS-SLE.pdf', bbox_inches='tight')