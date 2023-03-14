import sys
import os.path
import time
import pandas as pd
import numpy as np
import pyreadr
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit, GridSearchCV
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import RepeatedStratifiedKFold, RepeatedKFold
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.feature_selection import RFECV, VarianceThreshold
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score

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

# Split the data into features (X) and target (y)
X = df.iloc[:, 1:]
y = df.iloc[:, 0]

# Split the original dataset into a training set and a test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)

# Use StratifiedShuffleSplit() to split the training set into two additional subsets: 
# a subset for parameter tuning and a subset for final testing
cv = StratifiedShuffleSplit(n_splits=10, test_size=0.25, random_state=42)
train_index, tune_index = next(cv.split(X_train, y_train))

# Get the training and parameter tuning sets
X_train_final, y_train_final = X_train.iloc[train_index,], y_train.iloc[train_index]
X_tune, y_tune = X_train.iloc[tune_index], y_train.iloc[tune_index]

#scale data
scaler = MinMaxScaler()
X_train_final = pd.DataFrame(scaler.fit_transform(X_train_final), columns=X_train_final.columns)
X_tune = pd.DataFrame(scaler.fit_transform(X_tune), columns=X_tune.columns)
X_test = pd.DataFrame(scaler.transform(X_test), columns=X_test.columns)

# Remove features with low variance
selector = VarianceThreshold(threshold=0)
selector.fit(X_tune)
features = selector.get_feature_names_out(X_tune.columns)
X_tune = pd.DataFrame(selector.fit_transform(X_tune), columns=features)

# # Remove highly correlated features
# corr_matrix = X_tune.corr(method='spearman')
# upper_tri = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
# to_drop = [column for column in upper_tri.columns if any(upper_tri[column].abs() > 0.7)]
# X_tune = X_tune.drop(to_drop, axis=1)

# Perform feature selection and hyperparameter tuning using a pipeline
pipeline = Pipeline([
    ('classification', GridSearchCV(QuadraticDiscriminantAnalysis(), param_grid = {'reg_param': [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]},
                                    cv=RepeatedKFold(n_splits=len(X_tune.index), n_repeats=3, random_state=42), scoring='accuracy', n_jobs=-1))
])
# Fit the pipeline to the tune data
pipeline.fit(X_tune, y_tune)

# Set model with best parameters
clf = QuadraticDiscriminantAnalysis(reg_param=pipeline['classification'].best_params_['reg_param'])
# Match Train and Test features to Tune features
X_train_final = X_train_final.loc[:, X_tune.columns]
X_test = X_test.loc[:, X_tune.columns]
# Fit the model to the training data
clf.fit(X_train_final, y_train_final)
# Predict the test set
y_pred = clf.predict(X_test)

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
metrics.to_csv('exp.matrix/metrics/QDA_metrics_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Save confusion matrix to file
confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
confusion.to_csv('exp.matrix/metrics/QDA_confusion_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Print the AUROC curve
fpr, tpr, thresholds = roc_curve(y_test, y_pred)
plt.plot(fpr, tpr, label='AUC-ROC (area = %0.2f)' % auc)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('QDA: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.legend(loc="lower right")
plt.savefig('exp.matrix/AUROC/QDA_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Print the PR curve
precision, recall, thresholds = precision_recall_curve(y_test, y_pred)
average_precision = average_precision_score(y_test, y_pred)
disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
disp.plot()
disp.ax_.set_title('QDA: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.savefig('exp.matrix/PRC/QDA_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Save the model
import pickle
filename = 'ML.models/QDA_model_'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(clf, open(filename, 'wb'))

end_time = time.process_time()
cpu_time = end_time - start_time

print(f"CPU time used: {cpu_time:.2f} seconds")