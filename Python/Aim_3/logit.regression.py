# Import required libraries
import sys
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt

# Get the file name from the command line
file = sys.argv[0]

# Read in the expression matrix as a pandas dataframe
df = pd.read_csv(file, delimiter='\t', index_col=0)

# Replace classes with binary label
df['class'] = df['class'].replace({"OA": 0, "RA": 1})

# Split the data into features (X) and target (y)
X = df.iloc[:, 1:]
y = df.iloc[:, 0]

# Create the stratified sampling object
sss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)

# Loop through the splits
for train_index, test_index in sss.split(X, y):
    # Get the training and testing data
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = y.iloc[train_index], y.iloc[test_index]

# Perform feature selection with 10-fold cross validation
# Create the recursive feature eliminator that scores features by mean squared errors
rfecv = RFECV(estimator=LogisticRegression(), step=1, cv=StratifiedKFold(10), scoring='accuracy')

# Fit the recursive feature eliminator
rfecv.fit(X_train, y_train)

# Print the optimal number of features
print('Optimal number of features: {}'.format(rfecv.n_features_))

# Plot number of features VS. cross-validation scores
plt.figure()
plt.xlabel("Number of features selected")
plt.ylabel("Cross validation score (nb of correct classifications)")
plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)
plt.savefig('features.vs.cv.pdf')


# Print the features that are selected
print('Selected features: {}'.format(list(X_train.columns[rfecv.support_])))
# Print the features that are not selected
print('Not selected features: {}'.format(list(X_train.columns[~rfecv.support_])))

# Get the selected features
X_train = X_train.iloc[:, rfecv.support_]
X_test = X_test.iloc[:, rfecv.support_]

# Build the logistical regression model
model = LogisticRegression()
model.fit(X_train, y_train)

# Print the accuracy of the model
print(model.score(X_test, y_test))

# Get the predicted probabilities
y_pred = model.predict(X_test)

# Create a confusion matrix
conf_matrix = confusion_matrix(y_test, y_pred)


# Plot AUC curve
# Get the false positive rate, true positive rate, and thresholds
fpr, tpr, thresholds = roc_curve(y_test, y_pred)

# Get the area under the curve
auc = auc(fpr, tpr)

# Plot the ROC curve
plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % auc)
plt.plot([0, 1], [0, 1], 'k--')
plt.savefig('auc.pdf')

