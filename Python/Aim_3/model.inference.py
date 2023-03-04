import pickle as pk
import pyreadr
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score

# Read in best model file
best_model = pd.read_csv(sys.argv[1], sep='\t')

study1 = sys.argv[1].split('.')[0]
study2 = sys.argv[1].split('.')[1]

model_dict = {'LR': LogisticRegression(), 'RF': RandomForestClassifier(), 'GBM': GradientBoostingClassifier(), 'SVM': SVC()}

# Iterate through each model in best_model to evaluate new data
ML_model = best_model.iloc[0,0].split('_')[0] + '_model_'
cell_type = best_model.iloc[0,0].split('_')[1]
# Load data
df = pyreadr.read_r(study1 + '/exp.matrix/' + cell_type + '.RDS')
df = df[None]
print(df.head())
df['class'] = df['class'].replace({"control": 0, "disease": 1})

model = pk.load(open(study2 + '/ML.models/' + ML_model + cell_type + '.sav', 'rb'))

features = model.get_feature_names_out().tolist()

Y = df['class']
X = df[features]

y_pred = model.predict(X)
