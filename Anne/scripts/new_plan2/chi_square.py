import numpy as np
from scipy.stats.distributions import chi2
from bioinfokit import analys, visuz

def volcano_plot(df):
    visuz.GeneExpression.volcano(df=df,lfc='log2_fc',pv='p_value_X2',show=True)



def run_chi_square(df):
    df['counts_nonbreast_plus'] = df['counts_nonbreast'] + (1/1885)
    df['counts_breast_plus'] = df['counts_breast'] + (1/285)
    #
    n = np.array([285, 1885])
    S = df[['counts_breast_plus', 'counts_nonbreast_plus']].to_numpy() #+1
    #
    f = n / n.sum()
    C = (n-S)
    #
    E1 = S.sum(axis=1)[:,None]*f
    E2 = C.sum(axis=1)[:,None]*f
    #
    D1 = ((E1-S)**2/E1).sum(axis=1)
    D2 = ((E2-C)/E2).sum(axis=1)
    #
    X2 = D1 + D2
    #
    p_value_X2 = chi2.sf(X2,1)
    #
    constant = np.log(n[0]) - np.log(n[1])
    log2_fc = (np.log(S[:, 0]) - np.log(S[:, 1]) - constant) / np.log(2)
    #1
    #-1 borst kanker is helft voorkomen van nb
    log10_p_value = -np.log10(p_value_X2)
    #
    df['X2'] = X2
    df['p_value_X2'] = p_value_X2
    df['log10_p_value'] = log10_p_value # -log10(p-value)
    df['log2_fc'] = log2_fc # log2(FC)

