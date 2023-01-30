import pyreadr
import multiprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.feature_selection import RFECV

# Get the file name from the command line
file = sys.argv[0]

df = pyreadr.read_r(file)
df = df[None]

# Check if classes are named 'control' and 'disease' and replace with 0 and 1
if df['class'].unique().tolist() == ['control', 'disease']:
# Replace classes with binary label
    df['class'] = df['class'].replace({"control": 0, "disease": 1})
else:
    # replace class float with int
    df['class'] = df['class'].astype(int)

chrX = pyreadr.read_r("/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata")
chrX = chrX['chrX']

multiprocessing.cpu_count()/2

# Filter df column by chrX rownames and 'class'
cols_to_keep = ['class'] + chrX.index.tolist()
df2 = df.loc[:,df.columns.isin(cols_to_keep)]

# Split the data into features (X) and target (y)
X = df2.iloc[:, 1:]
y = df2.iloc[:, 0]

# Create the stratified sampling object
sss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)

# Loop through the splits
for train_index, test_index in sss.split(X, y):
    # Get the training and testing data
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = y.iloc[train_index], y.iloc[test_index]

# Create an RFECV object with a random forest classifier
clf = RandomForestClassifier()
rfecv = RFECV(clf, cv=5, scoring='accuracy', n_jobs=-1)

# Fit the RFECV object to the training data
X_train_selected = rfecv.fit_transform(X_train, y_train)
clf.fit(X_train_selected, y_train)

# Plot the number of features vs. cross-validation score
plt.figure()
plt.xlabel("Number of features selected")
plt.ylabel("Cross validation score (nb of correct classifications)")
plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)

# Print the optimal number of features
print('Optimal number of features: {}'.format(rfecv.n_features_))

# Print the features that are selected
print('Selected features: {}'.format(list(X_train.columns[rfecv.support_])))

# Predict the labels of the test data
X_test_selected = rfecv.transform(X_test)
y_pred = clf.predict(X_test_selected)

# Print the accuracy
print('Accuracy: {}'.format(accuracy_score(y_test, y_pred)))

# Print the confusion matrix
print('Confusion matrix: {}'.format(confusion_matrix(y_test, y_pred)))

# Print AUROC score
print('AUROC: {}'.format(roc_auc_score(y_test, y_pred)))

# Plot AUROC curce and save to file
fpr, tpr, thresholds = roc_curve(y_test, y_pred)
plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc_score(y_test, y_pred))
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.show()
plt.savefig('auroc.RF.png')







