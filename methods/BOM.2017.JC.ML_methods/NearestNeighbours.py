import csv
import numpy as np
from sklearn import preprocessing
from sklearn.neighbors import KNeighborsRegressor
from scipy.stats import pearsonr
from scipy.stats import spearmanr
# load training data
with open('/home/sam/trainingdata.csv', 'rb') as csvfile:
    molecule_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
    row = molecule_reader.next() # discard title row
    feature_names = np.array(row)
    X_train, y_train = [], []
    for row in molecule_reader:
        X_train.append(row)
        y_train.append(row[1]) # target value is "LogS"
    X_train = np.array(X_train)
    y_train = np.array(y_train)
# keep all descriptors
i = []
for j in range(2, 103):
    i.append(j)
X_train = X_train[:, i]
feature_names = feature_names[i]
# load test data
with open('/home/sam/testdata.csv', 'rb') as csvfile:
    molecule_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
    row = molecule_reader.next()
    feature_names = np.array(row)
    X_test, y_test = [], []
    for row in molecule_reader:
        X_test.append(row)
        y_test.append(row[1]) # target value is "LogS"
    X_test = np.array(X_test)
    y_test = np.array(y_test)
# keep all descriptors again
X_test = X_test[:, i]
feature_names = feature_names[i]
# convert to float type and scale the data
X_train = X_train.astype(np.float)
X_test = X_test.astype(np.float)
y_train = y_train.astype(np.float)
y_test = y_test.astype(np.float)
scaler = preprocessing.StandardScaler().fit(X_train)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)
# load and fit regressor
neigh= KNeighborsRegressor(n_neighbors=2)
neigh.fit(X_train, y_train)
preds = neigh.predict(X_test)
# print predictions, RMSE, Pearson correlation and Spearman rank
for x in range(25):
    print(preds[x])
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())
print(rmse(preds, y_test))
print(pearsonr(preds, y_test))
print(spearmanr(preds, y_test))


