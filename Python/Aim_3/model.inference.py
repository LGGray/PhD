import pickle
import pyreadr
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score

# Get cell type from input
cell_type = sys.argv[1]
# Get model from input
model = sys.argv[2]

# Load data
df = pyreadr.read_r(file)
df = df[None]
print(df.head())