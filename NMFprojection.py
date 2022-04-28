import argparse
import warnings

import pandas as pd
import numpy as np

from sklearn.decomposition import non_negative_factorization


seed = 0


def NMFprojection(X, fixed_W, normalized=False):
    # normalize X
    if normalized == False:
        X = (X * 10**4) / X.sum()
        X = np.log1p(X+1)

    if (normalized == False) & (X.max().max() > 20):
        warnings.warn("input X looks not normalized though normalized flag was passed")
    if (normalized == True) & (X.max().max() < 20):
        warnings.warn("input X looks normalized though normalized flag was not passed")

    # intersect genes
    genes = list(set(X.index) & set(fixed_W.index))
    X = X.loc[genes]
    fixed_W = fixed_W.loc[genes]
    fixed_H = fixed_W.T
    
    W, _, n_iter = non_negative_factorization(np.array(X.T, dtype = 'float64'), n_components=fixed_H.shape[0], 
                                          init='custom', random_state=seed, update_H=False, 
                                          H=np.array(fixed_H, dtype = 'float64'))
    H = W.T
    df_H = pd.DataFrame(H, columns=X.columns, index=fixed_W.columns)
    
    return X, df_H, fixed_W
    

def calc_RMSE(X, fixed_W, df_H):
    df_RMSE = pd.DataFrame(np.sqrt(((fixed_W.dot(df_H) - X) ** 2).sum(axis=0) / X.shape[0]), 
                            columns=['Error'], index=X.columns)
    
    print("stats of RSME")
    print(df_RMSE.describe())
    return df_RMSE


def calc_EV(X, fixed_W, df_H):
    l_ev = []
    for c in df_H.index:
        l_ev.append(1 - ((fixed_W[[c]].dot(df_H.loc[[c]]) - X) ** 2).sum().sum() / np.sum(np.array(X) ** 2))

    l_ev.append(1 - ((fixed_W.dot(df_H) - X) ** 2).sum().sum() / np.sum(np.array(X) ** 2))

    df_ev = pd.DataFrame(l_ev, columns=['ExplainedVariance'], index=list(df_H.index) + ['ALL'])
    print("stats of Explained Variance")
    print(df_ev)
    return df_ev


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'NMFprojection')
    parser.add_argument('input', 
                        help='input csv/tsv of gene expressions. UMI, scaledTPM, TPM can be used. Row: genes, Columns: samples.')
    parser.add_argument('fixedW', 
                        help='input csv/tsv of precomputed W. Row: genes, Columns: components.')
    parser.add_argument('--outputprefix', default='NMF',
                        help='output prefix. default=NMF')
    parser.add_argument('--normalized', action='store_true', 
                        help='if normalized and log transformed, specify this flag.')
    args = parser.parse_args()

    X = pd.read_csv(args.input, index_col=0, delim_whitespace=True) # gene x cell or samples
    fixed_W = pd.read_csv(args.fixedW, index_col=0) # gene' x components
    f_outputprefix = args.outputprefix

    X, df_H, fixed_W = NMFprojection(X, fixed_W, normalized=args.normalized)

    df_RMSE = calc_RMSE(X, fixed_W, df_H)
    df_ev = calc_EV(X, fixed_W, df_H)


    df_H.to_csv('{}_projection.csv'.format(f_outputprefix))
    df_RMSE.to_csv('{}_RMSE.csv'.format(f_outputprefix))
    df_ev.to_csv('{}_ExplainedVariance.csv'.format(f_outputprefix))
