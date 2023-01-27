import os
import sys
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import confusion_matrix, accuracy_score, classification_report
from sklearn import preprocessing
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")
sns.set(style="whitegrid", color_codes=True)

def randomForest(cell, data_dir, results_dir):
    os.chdir(data_dir)
    expr = pd.read_csv(cell + '.deg.csv', index_col=0)
    expr = expr.dropna(axis=0)
    expr = expr.dropna(axis=1)
    X = expr.drop('condition', axis=1)
    y = expr['condition']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    clf = RandomForestClassifier(n_estimators=100)
    clf.fit(X_train, y_train)
    print(clf.feature_importances_)
    print(clf.n_features_)
    print(clf.n_outputs_)
    print(clf.oob_score_)
    y_pred = clf.predict(X_test)
    print(confusion_matrix(y_test, y_pred))
    print(classification_report(y_test, y_pred))
    print(accuracy_score(y_test, y_pred))
    fpr, tpr, thresholds = roc_curve(y_test, y_pred)
    roc_auc = auc(fpr, tpr)
    plt.title('Receiver Operating Characteristic')
    plt.plot(fpr, tpr, 'b', label='AUC = %0.2f' % roc_auc)
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.savefig(results_dir + cell + '.AUROC.pdf')
    model = {'rf': clf, 'auc': roc_auc}
    model.to_csv(results_dir + cell + '.model.csv')

def main():
    cell = sys.argv[1]
    data_dir = sys.argv[2]
    results_dir = sys.argv[3]
    randomForest(cell, data_dir, results_dir)

if __name__ == '__main__':
    main()