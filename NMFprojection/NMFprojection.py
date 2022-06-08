import argparse
import warnings

import pandas as pd
import numpy as np

from sklearn.decomposition import non_negative_factorization


seed = 0


def NMFprojection(X, fixed_W, normalized=False, return_truncated=False):
    if X.index.duplicated().sum() > 0:
        raise ValueError("Gene names are duplicated!")

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
    X_trunc = X.loc[genes]
    fixed_W = fixed_W.loc[genes]
    fixed_H = fixed_W.T
    
    W, _, n_iter = non_negative_factorization(np.array(X_trunc.T, dtype = 'float64'), n_components=fixed_H.shape[0], 
                                          init='custom', random_state=seed, update_H=False, 
                                          H=np.array(fixed_H, dtype = 'float64'))
    H = W.T
    df_H = pd.DataFrame(H, columns=X.columns, index=fixed_W.columns)
    
    if return_truncated:
        return X, X_trunc, df_H, fixed_W
    else:
        return X, df_H, fixed_W
    

def calc_RMSE(X, fixed_W, df_H):
    df_RMSE = pd.DataFrame(np.sqrt(((fixed_W.dot(df_H) - X) ** 2).sum(axis=0) / X.shape[0]), 
                            columns=['Error'], index=X.columns)
    return df_RMSE


def calc_EV(X, fixed_W, df_H):
    l_ev = []
    for c in df_H.index:
        l_ev.append(1 - ((fixed_W[[c]].dot(df_H.loc[[c]]) - X) ** 2).sum().sum() / np.sum(np.array(X) ** 2))

    l_ev.append(1 - ((fixed_W.dot(df_H) - X) ** 2).sum().sum() / np.sum(np.array(X) ** 2))

    df_ev = pd.DataFrame(l_ev, columns=['ExplainedVariance'], index=list(df_H.index) + ['ALL'])
    return df_ev


def calc_hvg_overlap(X_norm, fixed_W, min_mean=0.0125, max_mean=3, min_disp=0.1, n_top_genes=500):
    import scanpy as sc
    tup_blacklist = ('TRAV', 'TRAJ', 'TRBV', 'TRBD', 'TRBJ',
                'TRGV', 'TRGJ', 'TRDV', 'TRDD', 'TRDJ',
                'IGKV', 'IGKJ', 'IGLV', 'IGLJ', 'IGHV', 'IGHD', 'IGHJ')
    if X_norm.index[0] == X_norm.index[0].capitalize():
        tup_blacklist = tuple([x.capitalize() for x in tup_blacklist])

    a = sc.AnnData(X_norm).T
    a = a[:,~a.var.index.str.startswith(tup_blacklist)]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sc.pp.highly_variable_genes(a, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp, n_bins=300)
    df_stats = a.var[['dispersions_norm', 'highly_variable', 'means']]
    df_stats['selected'] = df_stats.index.isin(fixed_W.index)

    if n_top_genes > 0:
        d = df_stats[df_stats.highly_variable]
        df_stats['highly_variable'] = df_stats.index.isin(list(d[d['dispersions_norm'].rank(ascending=False) <= n_top_genes].index))
  
    return df_stats

def main():
    parser = argparse.ArgumentParser(description = 'NMFprojection')
    parser.add_argument('input', 
                        help='input csv/tsv of gene expressions. UMI, scaledTPM, TPM can be used. Row: genes, Columns: samples.')
    parser.add_argument('fixedW', 
                        help='input csv/tsv of precomputed W. Row: genes, Columns: components.')
    parser.add_argument('--outputprefix', default='NMF',
                        help='output prefix. default=NMF')
    parser.add_argument('--normalized', action='store_true', 
                        help='if normalized and log transformed, specify this flag.')
    parser.add_argument('--min_mean', type=float, default=0.0125,
                        help='parameter for calculation of HVGs overlap')
    parser.add_argument('--max_mean', type=float, default=3,
                        help='parameter for calculation of HVGs overlap')
    parser.add_argument('--min_disp', type=float, default=0.1,
                        help='parameter for calculation of HVGs overlap.')
    parser.add_argument('--n_top_genes', default=500,
                        help='parameter for calculation of HVGs overlap.')
    parser.add_argument('--off_calc_hvg_overlap', action='store_true', 
                        help='turn off calc_hvg_overlap')
    parser.add_argument('--save_fullhvgstats', action='store_true', 
                        help='save full stats of hvg (hvg_overlap).')
    args = parser.parse_args()

    X = pd.read_csv(args.input, index_col=0, delim_whitespace=True) # gene x cell or samples
    fixed_W = pd.read_csv(args.fixedW, index_col=0) # gene' x components
    f_outputprefix = args.outputprefix

    X_norm, X_trunc, df_H, fixed_W_trunc = NMFprojection(X, fixed_W, normalized=args.normalized, return_truncated=True)
    df_H.to_csv('{}_projection.csv'.format(f_outputprefix))

    df_RMSE = calc_RMSE(X_trunc, fixed_W_trunc, df_H)
    print("## Stats of RSME")
    print(df_RMSE.describe())
    df_RMSE.to_csv('{}_RMSE.csv'.format(f_outputprefix))

    df_ev = calc_EV(X_trunc, fixed_W_trunc, df_H)
    print("## Stats of Explained Variance")
    print(df_ev)
    df_ev.to_csv('{}_ExplainedVariance.csv'.format(f_outputprefix))

    if not args.off_calc_hvg_overlap:
        df_stats = calc_hvg_overlap(X_norm, fixed_W_trunc, min_mean=args.min_mean, max_mean=args.max_mean, min_disp=args.min_disp,
                                    n_top_genes=args.n_top_genes)
        print('\n## Stats of overlap of HVGs')
        msg = 'Num. genes in fixed W: %s \n' % fixed_W.shape[0]
        msg += 'Num. Retained genes (Prop.): %s (%s)\n' % (fixed_W_trunc.shape[0], fixed_W_trunc.shape[0]/fixed_W.shape[0])
        msg += 'Prop. overlap of HVGs (POH) : {} in {} query HVGs'.format(
            df_stats.loc[df_stats['highly_variable'], 'selected'].sum() / df_stats['highly_variable'].sum(), 
            df_stats.highly_variable.sum())
        print(msg)
        if args.save_fullhvgstats:
            df_stats.to_csv('{}_hvgstats.csv'.format(f_outputprefix))
        with open('{}_hvgstats.txt'.format(f_outputprefix), 'w') as f:
            f.write(msg)

    df_ev.to_csv('{}_ExplainedVariance.csv'.format(f_outputprefix))

if __name__ == "__main__":
    main()
