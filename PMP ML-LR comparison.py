import sys
import skimage as sk
import numpy as np
import matplotlib;
import matplotlib.pyplot as plt;

plt.ion()
from matplotlib.ticker import MultipleLocator
import time, random
import inspect
import colorama as CLR;

CLR.init(autoreset=True)
import PyQt5
import scipy
import seaborn as sns
import pandas as pd
import itertools
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from skimage.io import imread, imsave
import cv2 as cv
from mpl_interactions import ioff, panhandler, zoom_factory
import copy
import shap
import psutil, os, gc;

PSUprocess = psutil.Process(os.getpid())  # RAM counter

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPalette, QColor, QFont
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QWidget, QScrollArea, QTabWidget, QPushButton, QCheckBox, \
    QHBoxLayout, QVBoxLayout

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import torch
import torch.nn as nn
import torchvision
import torchvision.transforms as transforms
from pytorch_tabnet.tab_model import TabNetClassifier

import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

import sklearn
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve, auc

from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV

import xgboost as xgb
from xgboost import XGBClassifier

import lifelines
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index as lifelines_CI

import xgbse as xgb;
from xgbse import XGBSEStackedWeibull, XGBSEKaplanNeighbors
from xgbse.converters import convert_to_structured
from xgbse.metrics import concordance_index as XGBSE_CI, approx_brier_score, dist_calibration_score

device = torch.device("cuda" if torch.cuda.is_available() else 'cpu')
printweights = False
pd.set_option('display.float_format', '{:.6f}'.format)
np.seterr(divide='ignore', invalid='ignore')

def get_col_type(col):
    if col.map(lambda x: isinstance(x, float) and not pd.isna(x)).any():
        return float
    else:
        return type(col.iloc[0])


na_vals = ['NA', 'null', '-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN', '#N/A N/A', '#N/A',
           'N/A', 'n/a', 'NA', 'null', 'NaN', '-NaN', 'nan', '-nan']

os.chdir("C:/filepath")
PMPraw = pd.read_csv('data_path.csv', na_values=na_vals, keep_default_na=False)
PMPdata = PMPraw[PMPraw['atrisk'] == 1]

FactorColumns = ['Age', 'Sex', 'Site', 'BMI', 'NumPregnancies',
                  'BD', 'Ulceration', 'HistoComType', 'Regression', 'MitoticRate', 'KidneyDisease', 'HeartDisease', 'RespDisease', 'Diabetes', 'Hypertension',
                  'Arthritis', 'Depression', 'HIV', 'Dementia',
                  'OtherMalignancy', 'Chemo5YR', 'Rad5YR', 'PrevMelanoma', 'Hormone5YR',
                  'PrevSkinCancer', 'FHxMelanoma', 'PreNaevus',
                  'StatinsLongTerm', 'NSAIDLongTerm', 'CorticosteroidLongTerm', 'PsychotropicsLongTerm',
                  'Occupation', 'Education', 'Partner', 'Birthplace',
                  'FACTG_Physical', 'FACTG_Social', 'FACTG_Emotional', 'FACTG_Functional', 'EverSmoker',
                  'SunOccup', 'Sunscreen', 'SelfSkinChecks', 'DoctorSkinChecks', 'HealthServiceAccess',
                  'RemotenessCat', 'SEIFADecile']
ContinuousVars = ['FACTG_Physical', 'FACTG_Social', 'FACTG_Emotional', 'FACTG_Functional', 'SEIFADecile'] # Age, BMI, NumPregnancies, BD, MitoticRate are all picked up by subsequent integer check
FactorColumnsRestricted = ['BD', 'Ulceration', 'HistoComType', 'Regression', 'MitoticRate', 'PreNaevus', 'Site']

restricted = False

if restricted:
    X_pre_filter = PMPdata[FactorColumnsRestricted].copy()
else:
    X_pre_filter = PMPdata[FactorColumns].copy()

ColumnType = {col: get_col_type(X_pre_filter[col]) for col in X_pre_filter}
OrdinalCheck = {col: ColumnType[col] not in [np.float64, np.int64, float, int] for col in X_pre_filter}
OrdinalCheck.update({var: False for var in ContinuousVars})

ReferenceCategories = {
    'HistoComType': 'SSM',
    'Site': 'Trunk',
    'MolesAge21': 'None',
    'Occupation': 'Unemployed',
    'Education': 'Less than Grade 12',
    'RemotenessCat': 'Major cities'
}

for col, ref in ReferenceCategories.items():
    if col in X_pre_filter.columns:
        other_categories = [c for c in X_pre_filter[col].unique() if c != ref]
        X_pre_filter[col] = pd.Categorical(
            X_pre_filter[col],
            categories=[ref] + other_categories
        )

X = pd.get_dummies(X_pre_filter, columns=[col for col in X_pre_filter if OrdinalCheck[col]], dtype=int, drop_first=True)
Y = PMPdata['rec7year']
Y_surv = PMPdata['meldeath7year']
Y_dm = PMPdata['dm7year']
feature_names = X.columns.tolist();

n_samples, n_features = X.shape;

SplitYear = PMPdata['CollectionYear'] < 2013
StageIB = np.logical_and(PMPdata['ClinicalStage'] == "IB", PMPdata['CollectionYear'] >= 2013)
StageIBIIA = np.logical_and(np.logical_or(PMPdata['ClinicalStage'] == "IB", PMPdata['ClinicalStage'] == "IIA"),
                            PMPdata['CollectionYear'] >= 2013)
StageIIBIIC = np.logical_and(np.logical_or(PMPdata['ClinicalStage'] == "IIB", PMPdata['ClinicalStage'] == "IIC"),
                             PMPdata['CollectionYear'] >= 2013)
StageIII = np.logical_and(PMPdata['PathologicStage'] == "III", PMPdata['CollectionYear'] >= 2013)
SLNB_Included = np.logical_and(PMPdata['SLNBTotal'] != "Not performed", PMPdata['CollectionYear'] >= 2013)

SLNBNegative = np.logical_and(PMPdata['SLNBpos'] == "Negative", PMPdata['CollectionYear'] >= 2013)
SLNBNotPerformed = np.logical_and(PMPdata['SLNBTotal'] == "Not performed", PMPdata['CollectionYear'] >= 2013)
DelayedRecur = np.logical_and(PMPdata['rec2year'] == 0, PMPdata['CollectionYear'] >= 2013)
DistantRecur = PMPdata['CollectionYear'] >= 2013
MelanomaMortality = PMPdata['CollectionYear'] >= 2013

X_original = X[SplitYear];
Y_original = Y[SplitYear];
PMP_original = PMPdata[SplitYear];
X_original_test = X[~SplitYear];
Y_original_test = Y[~SplitYear];

subset_conditions = {
    "train": SplitYear,
    "test": ~SplitYear,
    "early": StageIBIIA,
    "mid": StageIIBIIC,
    "late": StageIII,
    "SLNBincluded": SLNB_Included,
    "SLNBneg": SLNBNegative,
    "SLNBnone": SLNBNotPerformed,
    "DelayedRecur": DelayedRecur,
    "DistantRecur": DistantRecur,
    "MelanomaMortality": MelanomaMortality
}

for suffix, condition in subset_conditions.items():
    globals()[f"X_{suffix}"] = X[condition]
    globals()[f"Y_{suffix}"] = Y[condition]
    globals()[f"Y_surv_{suffix}"] = Y_surv[condition]

globals()[f"Y_MelanomaMortality"] = Y_surv[MelanomaMortality]
globals()[f"Y_DistantRecur"] = Y_dm[DistantRecur]

Y_train_dm = Y_dm[SplitYear] # n.b. this is just for plotting; all training and thresholding is performed using full recurrence data
Y_test_dm = Y_dm[~SplitYear]

X_dict = {}
X_scaled_dict = {}
Y_dict = {}
Y_surv_dict = {}

X_train = X[SplitYear]

sc = StandardScaler(); sc.fit(X_train)
names = list(subset_conditions.keys())

for name in names:
    X_dict[name] = torch.from_numpy(globals()["X_" + name].values.astype(np.float32))

    X_scaled_dict[name] = torch.from_numpy(sc.transform(X_dict[name]))

    Y_dict[name] = torch.from_numpy(globals()["Y_" + name].values.astype(np.float32))
    Y_dict[name] = Y_dict[name].view(Y_dict[name].shape[0], 1).ravel()

    Y_surv_dict[name] = torch.from_numpy(globals()["Y_surv_" + name].values.astype(np.float32))
    Y_surv_dict[name] = Y_surv_dict[name].view(Y_surv_dict[name].shape[0], 1).ravel()

X_grid = [X_dict["train"], X_dict["test"], X_dict["early"], X_dict["mid"], X_dict["late"], X_dict["SLNBincluded"], X_dict["SLNBneg"], X_dict["SLNBnone"], X_dict["DelayedRecur"], X_dict["DistantRecur"], X_dict["MelanomaMortality"]]
X_train, X_test, X_early, X_mid, X_late, X_SLNBincluded, X_SLNBneg, X_SLNBnone, X_DelayedRecur, X_DistantRecur, X_MelanomaMortality = X_grid

X_grid_scaled = [X_scaled_dict["train"], X_scaled_dict["test"], X_scaled_dict["early"], X_scaled_dict["mid"], X_scaled_dict["late"], X_scaled_dict["SLNBincluded"],
                 X_scaled_dict["SLNBneg"], X_scaled_dict["SLNBnone"], X_scaled_dict["DelayedRecur"], X_scaled_dict["DistantRecur"], X_scaled_dict["MelanomaMortality"]]
X_train_scaled, X_test_scaled, X_early_scaled, X_mid_scaled, X_late_scaled, X_SLNBincluded_scaled, X_SLNBneg_scaled, X_SLNBnone_scaled, X_DelayedRecur_scaled, X_DistantRecur_scaled, X_MelanomaMortality_scaled = X_grid_scaled

Y_surv_grid = [Y_surv_dict["train"], Y_surv_dict["test"], Y_surv_dict["early"], Y_surv_dict["mid"], Y_surv_dict["late"], Y_surv_dict["SLNBincluded"],
               Y_surv_dict["SLNBneg"], Y_surv_dict["SLNBnone"], Y_surv_dict["DelayedRecur"], Y_surv_dict["DistantRecur"], Y_surv_dict["MelanomaMortality"]]
Y_surv_train, Y_surv_test, Y_surv_early, Y_surv_mid, Y_surv_late, Y_surv_SLNBincluded, Y_surv_SLNBneg, Y_surv_SLNBnone, Y_surv_DelayedRecur, Y_surv_DistantRecur, Y_surv_MelanomaMortality = Y_surv_grid

Y_grid = [Y_dict["train"], Y_dict["test"], Y_dict["early"], Y_dict["mid"], Y_dict["late"], Y_dict["SLNBincluded"], Y_dict["SLNBneg"], Y_dict["SLNBnone"], Y_dict["DelayedRecur"], Y_dict["DistantRecur"], Y_dict["MelanomaMortality"]]
Y_train, Y_test, Y_early, Y_mid, Y_late, Y_SLNBincluded, Y_SLNBneg, Y_SLNBnone, Y_DelayedRecur, Y_DistantRecur, Y_MelanomaMortality = Y_grid

def Youden(Y_set, Y_probs, plotting=False, to_plot=False):
    fpr, tpr, thresholds = roc_curve(Y_set, Y_probs)

    optimal_index = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_index]
    auc_calc = sklearn.metrics.auc(fpr, tpr)

    if plotting:
        plt.figure(figsize=(8, 8))
        lw = 2
        plt.plot(fpr, tpr, color='navy',
                 lw=lw, label='AUC (area = %0.2f)' % auc_calc)
        plt.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('Specificity')
        plt.ylabel('Sensitivity')
        plt.title('Receiver Operating Characteristic (ROC) Curve')
        plt.legend(loc="lower right")
        plt.pause(0.001)
        # plt.show()

    if to_plot:
        return np.where(Y_probs >= optimal_threshold, 1, 0), optimal_threshold, fpr, tpr, thresholds
    else:
        return np.where(Y_probs >= optimal_threshold, 1, 0), optimal_threshold

def modelrun(modelname, model, scaling=False, modelarray = names,
             X_grid=X_grid, X_grid_scaled=X_grid_scaled, Y_grid=Y_grid,
             Y_train_override=False):
    X_clone = X_grid.copy();
    X_clone_scaled = X_grid_scaled.copy();
    Y_clone = Y_grid.copy()

    AUC_list = {'train': 0, 'test': 0, 'early': 0, 'mid': 0, 'late': 0, 'SLNBincluded': 0, 'SLNBneg': 0, 'SLNBnone': 0, 'DelayedRecur': 0, 'DistantRecur': 0, 'MelanomaMortality': 0, 'kfold': 0}

    for i in range(len(modelarray)):
        X_clone[i] = X_clone[i].numpy() if isinstance(X_clone[i], torch.Tensor) else X_clone[i]
        X_clone_scaled[i] = X_clone_scaled[i].numpy() if isinstance(X_clone_scaled[i], torch.Tensor) else \
            X_clone_scaled[i]

        if scaling == True:
            X_final = X_clone_scaled[i].copy()
        else:
            X_final = X_clone[i].copy()

        if modelarray[i] == 'train':
            Y_probs = model.predict_proba(X_final)[:, 1]
            Y_predicted, optimal_threshold = Youden(Y_clone[i], Y_probs)

        # we set threshold once, or otherwise it recalculates a new threshold for each validation
        AUC_list[modelarray[i]] = calc_probs(modelname, model, modelarray[i], X_final, Y_clone[i],
                                             direct_threshold=optimal_threshold)

def calc_probs(modelname, model, cohortname, X_final, Y_final, plotting=False, direct_threshold=-1, verbose=True):
    Y_probs = model.predict_proba(X_final)[:, 1]

    if (direct_threshold == -1):
        Y_predicted, final_threshold = Youden(Y_final, Y_probs, plotting)
    else:
        Y_predicted, final_threshold = np.where(Y_probs >= direct_threshold, 1, 0), direct_threshold

    PPV, NPV, Sens, Spec, F1, AUC, TotalPos, TotalPred, TotalN = get_metrics(Y_predicted, Y_probs, Y_final)

    if True:
        print(
            f'{modelname}, {cohortname} cohort (Pred {TotalPred}, Actual {TotalPos}/{TotalN}); Threshold: {final_threshold:.3f}, AUC: {AUC:.3f}, F1: {F1:.2f}, PPV: {PPV:.2f}, NPV: {NPV:.2f}, Sens: {Sens:.2f}, Spec: {Spec:.2f}')

    return AUC

def get_metrics(pred, prob, actual_tensor):
    actual = actual_tensor.numpy() if isinstance(actual_tensor, torch.Tensor) else actual_tensor.clone();
    pred = pred.squeeze();
    actual = actual.squeeze()

    TP = sum(np.logical_and(pred == 1, actual == 1))
    FP = sum(np.logical_and(pred == 1, actual == 0))

    TN = sum(np.logical_and(pred == 0, actual == 0))
    FN = sum(np.logical_and(pred == 0, actual == 1))

    TotalPred = sum(pred == 1)
    TotalPos = sum(actual == 1)
    TotalN = len(pred)

    PPV = TP / (TP + FP)
    NPV = TN / (TN + FN)
    Sens = TP / (TP + FN)
    Spec = TN / (TN + FP)

    AUC = roc_auc_score(actual, prob)

    F1 = (2 * Sens * PPV) / (Sens + PPV)

    return PPV, NPV, Sens, Spec, F1, AUC, TotalPos, TotalPred, TotalN

run_all = True
run_XGB = True

if True:
    modelLR = LogisticRegression(penalty = 'l2', max_iter = 10000, random_state = 40) # L2 regularization
    modelLR.fit(X_train_scaled, Y_train)

    modelrun("Logistic regression", modelLR, scaling=True)

    coefficients = modelLR.coef_
    if (printweights):
        for feature, coef in zip(feature_names, coefficients[0]):
            print(f"{feature}: {coef:.4f}")

    X_train_scaled_np = X_train_scaled.numpy()
    masker = shap.maskers.Independent(X_train_scaled_np)
    explainer_lr = shap.LinearExplainer(modelLR, masker=masker)
    shap_values_lr = explainer_lr(X_train_scaled_np)

if True:
    modelLRIso = LogisticRegression(penalty = None, max_iter = 10000, random_state = 0) # Non-regularized
    modelLRIso.fit(X_train_scaled, Y_train)

    modelrun("Logistic regression, non-regularized", modelLRIso, scaling=True)

    coefficients_iso = modelLRIso.coef_
    if (printweights):
        for feature, coef in zip(feature_names, coefficients_iso[0]):
            print(f"{feature}: {coef:.4f}")

if run_XGB:
    modelXGB = XGBClassifier(learning_rate = 0.1, max_depth = 5, min_child_weight = 1, subsample = 0.5, colsample_bytree = 0.5,
                             n_estimators = 50, alpha = 0, gamma = 1, reg_lambda = 1, random_state = 0)
    modelXGB.fit(X_train_scaled, Y_train)

    importances = modelXGB.feature_importances_ * 1000
    coef_names = list(X.columns)
    annotation = {name: importance for name, importance in zip(coef_names, importances)}
    annotation = dict(sorted(annotation.items(), key=lambda item: item[1], reverse=True))

    correlations = X_original.corrwith(PMP_original['rec7year'])

    for name, importance in annotation.items():
        correlation = correlations[name]
        if correlation <= 0:
            annotation[name] *= -1
        annotation[name] = round(annotation[name], 8)

    annotation_df = pd.DataFrame(list(annotation.items()), columns=['Feature', 'Importance'])
    annotation_df.to_csv('annotations.csv', index=False)

    modelrun("XGBoost classifier", modelXGB, scaling=True)

    explainer_xgb = shap.TreeExplainer(modelXGB)
    shap_values_xgb = explainer_xgb.shap_values(X_train_scaled_np)

    X_train_scaled_df = pd.DataFrame(X_train_scaled_np, columns=X.columns)
    shap_df = pd.DataFrame(shap_values_xgb, columns=X.columns)

    correlations = {}
    for feature in X.columns:
        corr = X_train_scaled_df[feature].corr(shap_df[feature])
        correlations[feature] = corr
    correlation_df = pd.DataFrame({
        'Feature': list(correlations.keys()),
        'Correlation': list(correlations.values())
    })

if True:
    train_probs_lr = modelLR.predict_proba(X_train_scaled)[:, 1]
    test_probs_lr = modelLR.predict_proba(X_test_scaled)[:, 1]

    train_probs_lr_iso = modelLRIso.predict_proba(X_train_scaled)[:, 1]
    test_probs_lr_iso = modelLRIso.predict_proba(X_test_scaled)[:, 1]

    train_probs_xgb = modelXGB.predict_proba(X_train_scaled)[:, 1]
    test_probs_xgb = modelXGB.predict_proba(X_test_scaled)[:, 1]

if True:
    mean_abs_shap_xgb = np.mean(np.abs(shap_values_xgb), axis=0)
    mean_abs_shap_lr = np.mean(np.abs(shap_values_lr.values), axis=0)

    df_shap = pd.DataFrame({
        'Feature': X.columns,
        'LR_Direction': ['Positive' if coef > 0 else 'Negative' for coef in coefficients[0]],
        'MeanAbsSHAP_XGB': mean_abs_shap_xgb,
        'MeanAbsSHAP_LR': mean_abs_shap_lr
    })

    df_shap = df_shap.merge(correlation_df, on='Feature', how='left')
    df_shap['XGB_Direction'] = df_shap['Correlation'].apply(lambda x: 'Positive' if x > 0 else 'Negative')
    df_shap.drop('Correlation', axis=1, inplace=True)

    df_shap.sort_values('MeanAbsSHAP_XGB', ascending=False, inplace=True)
    df_shap.reset_index(drop=True, inplace=True)

    df_shap.to_csv('SHAP.csv', index=False)

PMPdata_train = PMPdata[SplitYear]
PMPdata_test = PMPdata[~SplitYear]

df_train = pd.DataFrame({
    'ID': PMPdata_train["ID"].values,
    'Set': ['Train'] * len(X_train_scaled),
    'LR_Prob': train_probs_lr,
    'LRIso_Prob': train_probs_lr_iso,
    'XGB_Prob': train_probs_xgb,
    'Y_Actual': Y_train,
    'Y_DM': Y_train_dm,
    'Y_MelDeath': Y_surv_train,
    'Stage': PMPdata_train['ClinicalStage'].values,
    'SLNBStatus': PMPdata_train['SLNBTotal'].values,
    'DelayedRec': PMPdata_train['rec2year'].values
})

df_test = pd.DataFrame({
    'ID': PMPdata_test["ID"].values,
    'Set': ['Test'] * len(X_test_scaled),
    'LR_Prob': test_probs_lr,
    'LRIso_Prob': test_probs_lr_iso,
    'XGB_Prob': test_probs_xgb,
    'Y_Actual': Y_test,
    'Y_DM': Y_test_dm,
    'Y_MelDeath': Y_surv_test,
    'Stage': PMPdata_test['ClinicalStage'].values,
    'SLNBStatus': PMPdata_test['SLNBTotal'].values,
    'DelayedRec': PMPdata_test['rec2year'].values
})

output_filename = "XGB_res.csv" if restricted else "XGB.csv"
df_output = pd.concat([df_train, df_test], ignore_index=True)
df_output.to_csv(output_filename, index=False)

