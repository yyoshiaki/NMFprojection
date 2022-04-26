import argparse

import pandas as pd
import numpy as np

from sklearn.decomposition import non_negative_factorization


seed = 0


def NMFprojection(f_input, f_fixedW, f_outputprefix, normalized=False)
    X = pd.read_csv(f_input, delim_whitespace=True) # gene x cell or samples
    fixed_W = pd.read_csv(f_fixedW, delim_whitespace=True) # gene' x components
    
    # normalize X
    if normalized == False:
        X = (X * 10**4) / X.sum()
        X = np.log1p(X+1)
    
    # intersect genes
    genes = list(set(X.index) & set(gene.index))
    X = X.loc[gene]
    fixed_W = fixed_W.loc[gene]
    fixed_H = fixed_W.T
    
    W, _, n_iter = non_negative_factorization(X, n_components=2, init='custom', 
                                              random_state=seed, update_H=False, H=fixed_H)
    H = W.T
    df_H = pd.DataFrame(H, columns=X.columns, index=fixed_W.columns)
    df_RSS = pd.DataFrame(((W.dot(H) - X)**2).sum(axis=0), columns=['Error'], index=X.columns)
    
    print("stats of RSS")
    print(df_error.describe())
    
    H.to_csv('{}_projection.csv'.format(f_outputprefix))
    df_RSS.to_csv('{}_RSS.csv'.format(f_outputprefix))
    
    return df_H, df_error
    
    

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
    
    NMFprojection(args.input, args.fixedW, args.outputprefix, normalized=args.normalized)
