# Import required libraries
import sys
import os.path
import time
import pandas as pd
import numpy as np
import pyreadr
from boruta import BorutaPy
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RepeatedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score
from sklearn.feature_selection import RFECV
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.inspection import permutation_importance
import matplotlib.pyplot as plt

start_time = time.process_time()

# Get the file name from the command line
file = sys.argv[1]
print(os.path.basename(file))

# Read in expression RDS file
df = pyreadr.read_r(file)
df = df[None]
print(df.head())

# Replace classes with binary label
df['class'] = df['class'].replace({"control": 0, "disease": 1})

# Collect individual IDs
individuals = df['individual'].unique()
n_individuals = len(individuals)

# Get the number of individuals in each condition
individual_class = df['individual'].astype(str) + '_' + df['class'].astype(str)
n_control = len(individual_class[individual_class.str.endswith('_0')].unique())
n_disease = len(individual_class[individual_class.str.endswith('_1')].unique())

# Condition to factor in studies where classes are imbalanced
if n_control > 2:
    
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
else:
    # Randomly choose control individual
    np.random.seed(42)
    control_individuals = df['individual'][df['class'] == 0].unique()
    control_individual = np.random.choice(control_individuals, size=1, replace=False)
    disease_individuals = df['individual'][df['class'] == 1].unique()
    disease_individual = np.random.choice(disease_individuals, size=round(len(disease_individuals)*0.8), replace=False)

    # Get the corresponding cells for each dataset
    train_index = df['individual'].isin(control_individual.tolist() + disease_individual.tolist())
    test_index = ~df['individual'].isin(control_individual.tolist() + disease_individual.tolist())

    # Split the data into training, tuning, and testing sets
    X_train, X_test = df.loc[train_index,].drop(['class','individual'], axis=1), df.loc[test_index,].drop(['class','individual'], axis=1)
    y_train, y_test = df.loc[train_index,'class'], df.loc[test_index,'class']

    # Further split training into tuning and training sets
    X_train, X_tune, y_train, y_tune = train_test_split(X_train, y_train, test_size=0.25, random_state=42)

# Boruta feature selection
X = X_tune.values
y = y_tune.ravel()
# random forest classifier utilising all cores and sampling in proportion to y labels
rf = RandomForestClassifier(n_jobs=-1, class_weight='balanced', max_depth=5)
# define Boruta feature selection method
feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, random_state=1)
# find all relevant features - 5 features should be selected
feat_selector.fit(X, y)
# Return features
features = X_tune.columns[feat_selector.support_].tolist()

# Scale data. Required to speed up analysis with 'saga' solver
scaler = StandardScaler()
X_train = pd.DataFrame(scaler.fit_transform(X_train), columns=X_train.columns)
X_tune = pd.DataFrame(scaler.fit_transform(X_tune), columns=X_tune.columns)
X_test = pd.DataFrame(scaler.fit_transform(X_test), columns=X_test.columns)

# # Fit the RFECV object to the tune data
# rfecv = RFECV(RandomForestClassifier(n_jobs=-1), cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=0), scoring='accuracy', n_jobs=-1)
# rfecv = rfecv.fit(X_tune, y_tune)
# print('Model training complete')
# print('Optimal number of features: ', rfecv.n_features_)
# print('Best features: ', rfecv.get_feature_names_out().tolist())
# features = rfecv.get_feature_names_out().tolist()


# Tune the model to find the optimal C parameter
param_grid = {'C': [0.001, 0.01, 0.1, 1, 10]}
clf = LogisticRegression(solver='saga', penalty='elasticnet', l1_ratio=0.5, max_iter=10000, random_state=42, n_jobs=-1)
grid_search = GridSearchCV(clf, param_grid, cv=RepeatedKFold(n_splits=len(X_tune.index), n_repeats=3, random_state=0), scoring='accuracy', n_jobs=-1)
grid_search.fit(X_tune.loc[:, features], y_tune)

clf = grid_search.best_estimator_

clf.fit(X_train.loc[:, features], y_train)
# Predict the test set
y_pred = clf.predict(X_test.loc[:, features])
# Get the predicted probabilities
y_pred_proba = clf.predict_proba(X_test.loc[:, features])[:, 1]

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
metrics.to_csv('exp.matrix/metrics/logit_metrics_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Save confusion matrix to file
confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
confusion.to_csv('exp.matrix/metrics/logit_confusion_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Print the AUROC curve
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
plt.plot(fpr, tpr, label='AUC-ROC (area = %0.2f)' % auc)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('logit: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.legend(loc="lower right")
plt.savefig('exp.matrix/AUROC/logit_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Print the PR curve
precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
average_precision = average_precision_score(y_test, y_pred_proba)
disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
disp.plot()
disp.ax_.set_title('logit: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.savefig('exp.matrix/PRC/logit_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Save the model
import pickle
filename = 'ML.models/logit_model_'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(clf, open(filename, 'wb'))

end_time = time.process_time()
cpu_time = end_time - start_time

print(f"CPU time used: {cpu_time:.2f} seconds")