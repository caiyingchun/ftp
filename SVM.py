import os
import itertools
import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import VarianceThreshold
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV, cross_val_predict
from sklearn.metrics import precision_recall_fscore_support
import matplotlib.pyplot as plt
from sklearn.externals import joblib

def GetFileList(dir, fileList):
    newDir = dir
    if os.path.isfile(dir):
        fileList.append(dir.decode('gbk'))
    elif os.path.isdir(dir):
        for s in os.listdir(dir):
            newDir=os.path.join(dir,s)
            GetFileList(newDir, fileList)  
    return fileList

def removeLowVariableFeatures(X_train, X_test, X_ext):
    sel = VarianceThreshold().fit(X_train)
    np.savetxt('VT.txt', sel.variances_, fmt='%s', newline='\n')
    sel_X_train = sel.transform(X_train)
    sel_X_test = sel.transform(X_test)
    sel_X_ext = sel.transform(X_ext)
    return sel_X_train, sel_X_test, sel_X_ext

flist = []
for i in GetFileList('./KEGG/', []):
    flist.append(i.encode().replace('./KEGG/', ''))

for i in flist:
    X_data = pd.read_table('./KEGG/' + i, sep=',', header=0)
    X = X_data.ix[:, 3:].values
    y = X_data.ix[:, 2].values
    X_ext_data = pd.read_table('./Rhea/' + i, sep=',', header=0)
    X_ext = X_ext_data.ix[:, 3:].values
    y_ext = X_ext_data.ix[:, 2].values

    #train and test
    _, _, (X_train, X_test, y_train, y_test) = train_test_split(X, y, test_size=0.2, random_state=0)
    X_train, X_test, X_ext = removeLowVariableFeatures(X_train, X_test, X_ext)

    tuned_parameters = {'gamma': [2**x for x in range(3,-15,-2)],
                  'C': [2**x for x in range(-5,15,2)]}
    #tuned_parameters = {}
        
    score_cv = []
    score_test = []
    score_ext = []

    for j in range(5):
        clf0 = GridSearchCV(SVC(class_weight='balanced'), tuned_parameters, cv=10, n_jobs=8)
        clf0.fit(X_train, y_train)
        best_params = clf0.best_params_
        print 'best_params:' + str(best_params)
        clf = SVC(class_weight='balanced',**best_params)
        #####################CV################################
        y_true, y_pred = y_train, cross_val_predict(clf, X_train, y_train, cv=10)
        prec, recall, f1score, _ = precision_recall_fscore_support(y_true, y_pred, average='weighted')
        score_cv.append([prec, recall, f1score])

        clf.fit(X_train, y_train)
        #####################test################################
        y_true, y_pred = y_test, clf.predict(X_test)
        prec, recall, f1score, _ = precision_recall_fscore_support(y_true, y_pred, average='weighted')
        score_test.append([prec, recall, f1score])

        #####################ext################################
        y_true, y_pred = y_ext, clf.predict(X_ext)
        prec, recall, f1score, _ = precision_recall_fscore_support(y_true, y_pred, average='weighted')
        score_ext.append([prec, recall, f1score])
        
    print 'SVM_' + i
    ###########Output_CV##############################
    score_cv_mean = np.array(score_cv).mean(axis=0)
    score_cv_std = np.array(score_cv).std(axis=0)

    for k in [i for j in zip(score_cv_mean, score_cv_std) for i in j]:
        print k,
    print
     ###########Output_test#############################
    score_test_mean = np.array(score_test).mean(axis=0)
    score_test_std = np.array(score_test).std(axis=0)

    for k in [i for j in zip(score_test_mean, score_test_std) for i in j]:
        print k,
    print
    ###########Output_ext##############################
    score_ext_mean = np.array(score_ext).mean(axis=0)
    score_ext_std = np.array(score_ext).std(axis=0)

    for k in [i for j in zip(score_ext_mean, score_ext_std) for i in j]:
        print k,
    print
