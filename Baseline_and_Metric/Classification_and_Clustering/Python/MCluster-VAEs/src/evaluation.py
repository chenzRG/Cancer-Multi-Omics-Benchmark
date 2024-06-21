"""
推荐使用R来进行，用python实在太慢了。
因为python中的log rank test和scipy上的检验，都是利用numpy实现的，
底层并不是C++或fortran

这里保留一个简单版本，便于实时评测
"""

import multiprocessing as mp
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import kruskal, chi2_contingency, binomtest
from lifelines.statistics import multivariate_logrank_test
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment
from scipy.special import comb
from sklearn import metrics as M


CLINS = [
    "survival", "pathologic_T", "pathologic_N", "pathologic_M",
    "pathologic_stage", "age_at_initial_pathologic_diagnosis", "gender"
]


def cluster_acc(Y_pred, Y):
    assert Y_pred.size == Y.size
    tab = pd.crosstab(Y_pred, Y).values
    if tab.shape[0] != tab.shape[1]:
        d = max(tab.shape)
        w = np.zeros((d, d), dtype=np.int64)
        w[:tab.shape[0], :tab.shape[1]] = tab
    else:
        w = tab
    rind, cind = linear_sum_assignment(w.max() - w)
    score = w[rind, cind].sum() / Y_pred.size
    return score


def _get_pairs(arr):
    uni = np.unique(arr)
    res = set()
    for a in uni:
        inds = np.nonzero(arr == a)[0]
        if len(inds) > 1:
            s = set(combinations(inds, 2))
            res = res.union(s)
    return res


def cluster_f1(trues, preds):
    trues_pairs = _get_pairs(trues)
    preds_pairs = _get_pairs(preds)
    a = len(trues_pairs.intersection(preds_pairs))
    b = len(trues_pairs.union(preds_pairs))
    return (2 * a) / (a + b)


def cluster_f1_(trues, preds):
    """
    其得到的结果和前一个函数一致，并且拥有更快的速度。
    两者的速度分析在test中有:test_f1.png
    """
    tab = pd.crosstab(trues, preds).values
    a = tab.sum(axis=0)
    b = tab.sum(axis=1)
    nij = comb(tab, 2)
    ai = comb(a, 2)
    bj = comb(b, 2)
    return 2 * nij.sum() / (ai.sum() + bj.sum())


def cluster_metrics(trues, preds, f1=True):
    ind = pd.isnull(trues)
    if ind.sum() > 0:
        trues, preds = trues[~ind], preds[~ind]
    acc = cluster_acc(preds, trues)
    ari = M.adjusted_rand_score(trues, preds)
    nmi = M.normalized_mutual_info_score(trues, preds)
    if not f1:
        return {"acc": acc, "ari": ari, "nmi": nmi}
    f1 = cluster_f1_(trues, preds)
    return {"acc": acc, "ari": ari, "nmi": nmi, "f1": f1}


def clinical_metrics(preds, meta_dict):
    res = {}
    if ("days" in meta_dict) and ("status" in meta_dict):
        # survival log rank test
        # values = df[["days", "status"]].fillna(0.).values  # 来自原论文里的处理
        # 感觉上面的处理还是有些不妥，这里改正一下
        surv = pd.DataFrame({"days": meta_dict["days"],
                             "status": meta_dict["status"]})
        ind = surv[["days", "status"]].notna().all(axis=1).values
        values = surv[["days", "status"]].loc[ind, :].values
        surv_p = logrank_pvalue(values, preds[ind])
        res["survival"] = surv_p

    # others
    for col in [
        "gender",
        "age_at_initial_pathologic_diagnosis",
        "pathologic_M",
        "pathologic_N",
        "pathologic_T",
        "pathologic_stage"
    ]:
        if col not in meta_dict:
            continue
        cli_var = meta_dict[col]
        # 临床特征中可能存在缺失值
        ind = pd.isnull(cli_var)
        if ind.all():
            # 如果一个临床特征全是NaN，则跳过
            continue
        clu_var_nona, cli_var_nona = preds[~ind], cli_var[~ind]
        if col.startswith("age"):
            # 只有age去做kwh
            res[col] = kwh_pvalue(cli_var_nona, clu_var_nona)
        else:
            # 如果临床特征达不到2类及以上，则跳过
            if np.unique(cli_var_nona).shape[0] >= 2:
                res[col] = chi2_pvalue(cli_var_nona, clu_var_nona)

    res = {k: -np.log10(v + 1e-30) for k, v in res.items()}
    return res


def conditional_entropy(probs):
    logp = np.log(probs + 1e-10)
    return -np.mean(np.sum(probs * logp, axis=1))


def chi2_pvalue(values, group):
    tab = pd.crosstab(values, group).values
    return chi2_contingency(tab, correction=False)[1]


def chi2_perm_pvalue(values, group):
    group = np.random.permutation(group)
    return chi2_pvalue(values, group)


def chi2_pvalues_perm_mp(values, group, niters=1000, ncores=20):
    with mp.Pool(ncores) as pool:
        res = pool.starmap(chi2_perm_pvalue,
                           zip([values] * niters, [group] * niters))
    return res


def kwh_pvalue(values, group):
    k = np.unique(group)
    if len(k) < 2:
        return np.NaN
    args = [values[group == i] for i in k]
    return kruskal(*args).pvalue


def kwh_perm_pvalue(values, group):
    group = np.random.permutation(group)
    return kwh_pvalue(values, group)


def kwh_pvalues_perm_mp(values, group, niters=1000, ncores=20):
    with mp.Pool(ncores) as pool:
        res = pool.starmap(kwh_perm_pvalue,
                           zip([values] * niters, [group] * niters))
    return res


def logrank_pvalue(values, group):
    days, status = values[:, 0], values[:, 1]
    return multivariate_logrank_test(days, group, status).p_value


def logrank_perm_pvalue(values, group):
    group = np.random.permutation(group)
    return logrank_pvalue(values, group)


def logrank_pvalues_perm_mp(values, group, niters=1000, ncores=20):
    with mp.Pool(ncores) as pool:
        res = pool.starmap(logrank_perm_pvalue,
                           zip([values] * niters,
                               [group] * niters))
    return res


def permutation_test(values, group, ncores=20, nrepeats=1000, method="chi2"):
    orip = {"chi2": chi2_pvalue,
            "kwh": kwh_pvalue,
            "logrank": logrank_pvalue}[method](values, group)
    # orip = chi2_pvalue(values, group)
    perm_ps = {
        "chi2": chi2_pvalues_perm_mp,
        "kwh": kwh_pvalues_perm_mp,
        "logrank": logrank_pvalues_perm_mp
    }[method](values, group, nrepeats, ncores)
    # perm_ps = chi2_pvalues_perm_mp(values, group, nrepeats, ncores)
    nextreme = (perm_ps <= orip).sum()
    binom_res = binomtest(nextreme, nrepeats)
    return binom_res.pvalue


def plot_for_clin_test(res, orders=None, save_fn=None):
    res = res.groupby(["dat", "method"]).agg({
        "-logp": np.nanmean,
        "pvalue": lambda x: (x < 0.05).sum()
    }).rename(columns={"pvalue": "num.sig"})
    res = res.reset_index()
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    sns.stripplot(data=res, x="-logp", y="method", hue="dat",
                  order=orders, ax=axs[0, 0])
    sns.stripplot(data=res, x="num.sig", y="method", hue="dat",
                  order=orders, ax=axs[0, 1])
    sns.boxplot(data=res, x="-logp", y="method",
                order=orders, ax=axs[1, 0])
    sns.boxplot(data=res, x="num.sig", y="method",
                order=orders, ax=axs[1, 1])
    fig.tight_layout()
    if save_fn is not None:
        fig.savefig(save_fn)
    return fig, axs


def plot_for_surv_box(res, orders=None, save_fn=None):
    res = res[res["clin_vars"] == "survival"]
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
    sns.stripplot(data=res, x="-logp", y="method", hue="dat",
                  order=orders, ax=axs[0])
    sns.boxplot(data=res, x="-logp", y="method",
                order=orders, ax=axs[1])
    fig.tight_layout()
    if save_fn is not None:
        fig.savefig(save_fn)
    return fig, axs
