import numpy as np
import pandas as pd
import os
import statsmodels.api as sm
import argparse
from lpa_pipeline import validator, vcf_translate, preprocessing, utils, settings

def load_simple(ETHNICITY):
    genotype = pd.read_csv(settings.CSV_SIMPLE[ETHNICITY], index_col = 0)
    snp_x = pd.read_csv(settings.SNP_X_SIMPLE[ETHNICITY], index_col = 0)
    snp_y = pd.read_csv(settings.SNP_Y_SIMPLE[ETHNICITY], index_col = 0)
    df = pd.concat([genotype, snp_x], axis = 1)
    df, snp_y = df.align(snp_y, axis = 0)
    return df, snp_y

def DEMENTIA_GLM_Binomial(ETHNICITY):
    label = "DEMENTIA"    
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Binomial"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")
    df, snp_y = load_simple(ETHNICITY)
    GLM = sm.GLM(endog = snp_y[label], exog = sm.add_constant(df), family=sm.families.Binomial())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")

def DEMENTIA_GLM_Poisson(ETHNICITY):
    label = "DEMENTIA"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Poisson"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")
    df, snp_y = load_simple(ETHNICITY)
    GLM = sm.GLM(endog = snp_y[label], exog = sm.add_constant(df), family=sm.families.Poisson())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")

def DEMENTIA_MNLogit(ETHNICITY):
    label = "DEMENTIA"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/MNLogit"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")
    df, snp_y = load_simple(ETHNICITY)
    MNL = sm.MNLogit(endog = snp_y[label], exog = sm.add_constant(df))
    # for the matrix is singular, run this
    results = MNL.fit_regularized(method = "l1")
    # and we only have the bfgs param, similar to https://groups.google.com/g/pystatsmodels/c/E3VOUYxms6A?pli=1
    results.params.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    
def DIABETES_GLM_Binomial(ETHNICITY):
    label = "DIABETES"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Binomial"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    GLM = sm.GLM(endog = snp_y[label], exog = sm.add_constant(df), family=sm.families.Binomial())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")

def DIABETES_Logit(ETHNICITY):
    label = "DIABETES"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/Logit"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    logit = sm.Logit(endog = snp_y[label], exog = sm.add_constant(df))
    results = logit.fit_regularized(method='l1')
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def HEART_GLM_Binomial(ETHNICITY):
    label = "HEART"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Binomial"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    GLM = sm.GLM(endog = snp_y[label], exog = sm.add_constant(df), family=sm.families.Binomial())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def HEART_Logit(ETHNICITY):
    label = "HEART"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/Logit"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    logit = sm.Logit(endog = snp_y[label], exog = sm.add_constant(df))
    # for the matrix is singular, run this
    results = logit.fit_regularized(method='l1')
    # and we only have the bfgs param, similar to https://groups.google.com/g/pystatsmodels/c/E3VOUYxms6A?pli=1
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def HYPERTENSION_GLM_Binomial(ETHNICITY):
    label = "HYPERTENSION"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Binomial"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)

    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    GLM = sm.GLM(endog = snp_y[label], exog = sm.add_constant(df), family=sm.families.Binomial())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def HYPERTENSION_Logit(ETHNICITY):
    label =  "HYPERTENSION"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/Logit"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")
    
    df, snp_y = load_simple(ETHNICITY)
    logit = sm.Logit(endog = snp_y[label], exog = sm.add_constant(df))
    # for the matrix is singular, run this
    results = logit.fit_regularized(method='l1')
    # and we only have the bfgs param, similar to https://groups.google.com/g/pystatsmodels/c/E3VOUYxms6A?pli=1
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
    
def LIP01_B_GLM_Gaussian(ETHNICITY):
    label = "LIP01_B"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def LIP01_B_GLM_Gamma(ETHNICITY):
    label = "LIP01_B"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gamma"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    glm = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gamma())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def LIP01_B_GLM_InverseGaussian(ETHNICITY):
    label = "LIP01_B"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_InverseGaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.InverseGaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")


    

    
def LIP02_B_GLM_Gaussian(ETHNICITY):
    label = "LIP02_B"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")

    
    
def LIP02_B_GLM_Gamma(ETHNICITY):
    label = "LIP02_B"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gamma"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    glm = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gamma())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def LIP02_B_GLM_InverseGaussian(ETHNICITY):
    label = "LIP02_B"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_InverseGaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.InverseGaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    

    
def LIP03_B_GLM_Gaussian(ETHNICITY):
    label = "LIP03_B"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def LIP03_B_GLM_Gamma(ETHNICITY):
    label = "LIP03_B"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gamma"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    glm = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gamma())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def LIP03_B_GLM_InverseGaussian(ETHNICITY):
    label = "LIP03_B"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_InverseGaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)

    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.InverseGaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
    

    
def LIP04_2_GLM_Gaussian(ETHNICITY):
    label = "LIP04_2"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)

    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def LIP04_2_GLM_Gamma(ETHNICITY):
    label = "LIP04_2"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gamma"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)

    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    glm = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gamma())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def LIP04_2_GLM_InverseGaussian(ETHNICITY):
    label = "LIP04_2"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_InverseGaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)

    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.InverseGaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")

    
    
def INSL01_2_GLM_Gaussian(ETHNICITY):
    label = "INSL01_2"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)

    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")

    
    
def INSL01_2_GLM_Gamma(ETHNICITY):
    label = "INSL01_2"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gamma"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)

    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    glm = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gamma())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")

    
    
def INSL01_2_GLM_InverseGaussian(ETHNICITY):
    label = "INSL01_2"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_InverseGaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)

    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.InverseGaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")

    
    
    

def INSL02_2_GLM_Gaussian(ETHNICITY):
    label = "INSL02_2"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)

    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def INSL02_2_GLM_Gamma(ETHNICITY):
    label = "INSL02_2"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gamma"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")

    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    glm = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gamma())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def INSL02_2_GLM_InverseGaussian(ETHNICITY):
    label = "INSL02_2"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_InverseGaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")
        
    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.InverseGaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    

def INSL03_2_GLM_Gaussian(ETHNICITY):
    label = "INSL03_2"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")
        
    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def INSL03_2_GLM_Gamma(ETHNICITY):
    label = "INSL03_2"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gamma"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")
        
    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    glm = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gamma())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def INSL03_2_GLM_InverseGaussian(ETHNICITY):
    label = "INSL03_2"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_InverseGaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")
        
    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.InverseGaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def HBA1C_2_GLM_Gaussian(ETHNICITY):
    label = "HBA1C_2"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")
        
    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def HBA1C_2_GLM_Gamma(ETHNICITY):
    label = "HBA1C_2"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_Gamma"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")
        
    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    glm = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.Gamma())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        
    
    
def HBA1C_2_GLM_InverseGaussian(ETHNICITY):
    label = "HBA1C_2"
    rel_path = f"/data_analysis_result/pysm_Lasso/{label}/GLM_InverseGaussian"
    abs_path = settings.LPA_SHARED + rel_path
    os.makedirs(abs_path, exist_ok = True)
    
    print(f"label = {label} \n method = {rel_path.split("/")[-1]}")
        
    df, snp_y = load_simple(ETHNICITY)
    snp_y = snp_y[label]
    snp_y = snp_y.dropna()
    snp_y, df = snp_y.align(df, axis = 0, join = "inner")
    GLM = sm.GLM(endog = snp_y, exog = sm.add_constant(df), family=sm.families.InverseGaussian())
    results = GLM.fit_regularized(method='elastic_net', alpha=1, refit = True)
    output = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    output.columns = ["weight", "standard error", "p-value"]
    output.to_csv(abs_path + f"/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
    print(f"ETHNICITY:{ETHNICITY}: finished, saved to {rel_path}/ETHNICITY={ETHNICITY},n={df.shape[0]}.csv")
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser() 
    parser.add_argument("-i", "--index", default="1", help="the index of the task")
    args = parser.parse_args()
    eth_dict = {0: "complete", 1:1,2:2,3:3}
    task = args.index // 4
    eth = args.index % 4
    if task == 0:
        DEMENTIA_GLM_Binomial(eth_dict[eth])
    elif task == 1:
        DEMENTIA_GLM_Poisson(eth_dict[eth])
    elif task == 2:
        DEMENTIA_MNLogit(eth_dict[eth])
    elif task == 3:
        DIABETES_GLM_Binomial(eth_dict[eth])
    elif task == 4:
        DIABETES_Logit(eth_dict[eth])
    elif task == 5:
        HEART_GLM_Binomial(eth_dict[eth])
    elif task == 6:
        HEART_Logit(eth_dict[eth])
    elif task == 7:
        HYPERTENSION_GLM_Binomial(eth_dict[eth])
    elif task == 8:
        HYPERTENSION_Logit(eth_dict[eth])
    elif task == 9:
        LIP01_B_GLM_Gaussian(eth_dict[eth])
    elif task == 10:
        LIP01_B_GLM_Gamma(eth_dict[eth])
    elif task == 11:
        LIP01_B_GLM_InverseGaussian(eth_dict[eth])
    elif task == 12:
        LIP02_B_GLM_Gaussian(eth_dict[eth])
    elif task == 13:
        LIP02_B_GLM_Gamma(eth_dict[eth])
    elif task == 14:
        LIP02_B_GLM_InverseGaussian(eth_dict[eth])
    elif task == 15:
        LIP03_B_GLM_Gaussian(eth_dict[eth])
    elif task == 16:
        LIP03_B_GLM_Gamma(eth_dict[eth])
    elif task == 17:
        LIP03_B_GLM_InverseGaussian(eth_dict[eth])
    elif task == 18:
        LIP04_2_GLM_Gaussian(eth_dict[eth])
    elif task == 19:
        LIP04_2_GLM_Gamma(eth_dict[eth])
    elif task == 20:
        LIP04_2_GLM_InverseGaussian(eth_dict[eth])
    elif task == 21:
        INSL01_2_GLM_Gaussian(eth_dict[eth])
    elif task == 22:
        INSL01_2_GLM_Gamma(eth_dict[eth])
    elif task == 23:
        INSL01_2_GLM_InverseGaussian(eth_dict[eth])
    elif task == 24:
        INSL02_2_GLM_Gaussian(eth_dict[eth])
    elif task == 25:
        INSL02_2_GLM_Gamma(eth_dict[eth])
    elif task == 26:
        INSL02_2_GLM_InverseGaussian(eth_dict[eth])
    elif task == 27:
        INSL03_2_GLM_Gaussian(eth_dict[eth])
    elif task == 28:
        INSL03_2_GLM_Gamma(eth_dict[eth])
    elif task == 29:
        INSL03_2_GLM_InverseGaussian(eth_dict[eth])
    elif task == 30:
        HBA1C_2_2_GLM_Gaussian(eth_dict[eth])
    elif task == 31:
        HBA1C_2_GLM_Gamma(eth_dict[eth])
    elif task == 32:
        HBA1C_2_GLM_InverseGaussian(eth_dict[eth])
