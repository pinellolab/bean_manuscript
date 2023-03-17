from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.metrics

palette = {
    "mageck sort": "orange",
    "mageck sort, fitdisp": "chocolate",
    "mageck top/bot": "silver",
    "mageck top/bot, fitdisp": "grey",
    "model A": "orchid",
    "model B0": "yellowgreen",
    "model B": "seagreen",
    "model B2": "c",
}


def find_nearest(array, value):
    # https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def plot_auc_(
    ax,
    z,
    pos_idx,
    neg_idx,
    critical_value=1.96,
    plot_cv_legend=False,
    show_auroc=True,
    **kwargs
):
    y_eval = np.concatenate([z[pos_idx], z[neg_idx]])
    label = np.concatenate([np.ones(len(pos_idx)), np.zeros(len(neg_idx))])
    auroc = sklearn.metrics.roc_auc_score(label, y_eval)
    if show_auroc:
        kwargs["label"] = "{} ({:.3f})".format(kwargs["label"], auroc)

    fpr, tpr, thres = sklearn.metrics.roc_curve(label, y_eval)
    critical_thres_idx = [find_nearest(thres, critical_value)]

    xpt = np.linspace(0, 1, 100)
    ax.plot(fpr, tpr, **kwargs)
    if plot_cv_legend:
        kwargs["label"] = "z={}".format(critical_value)
    else:
        kwargs["label"] = ""
    ax.scatter(fpr[critical_thres_idx], tpr[critical_thres_idx], **kwargs)

    if plot_cv_legend:
        kwargs["label"] = "z={}".format(0)
    zero_idx = [find_nearest(thres, 0)]
    ax.scatter(fpr[zero_idx], tpr[zero_idx], marker="x", **kwargs)
    ax.plot(xpt, xpt, linestyle="--", color="black")
    return auroc


def plot_prc_(
    ax,
    z,
    pos_idx,
    neg_idx,
    critical_value=1.96,
    plot_cv_legend=False,
    show_auc=True,
    **kwargs
):
    y_eval = np.concatenate([z[pos_idx], z[neg_idx]])
    label = np.concatenate([np.ones(len(pos_idx)), np.zeros(len(neg_idx))])

    prec, recall, thres = sklearn.metrics.precision_recall_curve(label, y_eval)
    critical_thres_idx = [find_nearest(thres, critical_value)]
    auprc = sklearn.metrics.auc(recall, prec)
    if show_auc:
        kwargs["label"] = "{} ({:.3f})".format(kwargs["label"], auprc)

    ax.plot(recall, prec, **kwargs)
    if plot_cv_legend:
        kwargs["label"] = "z={}".format(critical_value)
    else:
        kwargs["label"] = ""
    ax.scatter(recall[critical_thres_idx], prec[critical_thres_idx], **kwargs)
    if plot_cv_legend:
        kwargs["label"] = "z={}".format(0)
    zero_idx = [find_nearest(thres, 0)]
    ax.scatter(recall[zero_idx], prec[zero_idx], marker="x", **kwargs)
    ax.axhline(
        y=len(pos_idx) / (len(pos_idx) + len(neg_idx)), linestyle="--", color="black"
    )
    return auprc


def get_f_score_inc(z, pos_inc_idx, neg_idx, critical_value=1.96):
    y_eval = np.concatenate([z[pos_inc_idx], z[neg_idx]]) >= critical_value
    label = np.concatenate([np.ones(len(pos_inc_idx)), np.zeros(len(neg_idx))])
    return sklearn.metrics.f1_score(label, y_eval)


def get_f_score_dec(z, pos_dec_idx, neg_idx, critical_value=1.96):
    y_eval = np.concatenate([z[pos_dec_idx], z[neg_idx]]) < -critical_value
    label = np.concatenate([np.ones(len(pos_dec_idx)), np.zeros(len(neg_idx))])
    return sklearn.metrics.f1_score(label, y_eval)


def _get_f_OvO(z, pos_inc_idx, pos_dec_idx, neg_ctrl_idx, critical_value=1.65):
    f_inc = get_f_score_inc(z, pos_inc_idx, neg_ctrl_idx, critical_value=critical_value)
    f_dec = get_f_score_dec(z, pos_dec_idx, neg_ctrl_idx, critical_value=critical_value)
    return (f_inc, f_dec)


def get_f(
    list_zs,
    z_labels,
    pos_inc_idx,
    pos_dec_idx,
    neg_ctrl_idx,
    prefix="",
    critical_value=1.65,
):
    f_incs = list(
        map(
            lambda z: get_f_score_inc(
                z, pos_inc_idx, neg_ctrl_idx, critical_value=critical_value
            ),
            list_zs,
        )
    )
    f_decs = list(
        map(
            lambda z: get_f_score_dec(
                z, pos_dec_idx, neg_ctrl_idx, critical_value=critical_value
            ),
            list_zs,
        )
    )
    df = pd.DataFrame(
        {"{}_inc".format(prefix): f_incs, "{}_dec".format(prefix): f_decs},
        index=z_labels,
    )
    return df


def _get_precision_inc(
    z, fdr, pos_inc_idx, neg_idx, critical_value=0.1, filter_sign=False
):
    y_eval = np.concatenate([fdr[pos_inc_idx], fdr[neg_idx]]) < critical_value
    if filter_sign:
        y_eval = np.logical_and(
            y_eval,
            np.concatenate([z[pos_inc_idx] > 0, np.ones(len(neg_idx), dtype=bool)]),
        )
    if y_eval.sum() == 0:
        return np.nan
    label = np.concatenate([np.ones(len(pos_inc_idx)), np.zeros(len(neg_idx))])
    return sklearn.metrics.precision_score(label, y_eval)


def _get_precision_dec(
    z, fdr, pos_dec_idx, neg_idx, critical_value=0.1, filter_sign=False
):
    y_eval = np.concatenate([fdr[pos_dec_idx], fdr[neg_idx]]) < critical_value
    if filter_sign:
        y_eval = np.logical_and(
            y_eval,
            np.concatenate([z[pos_dec_idx] < 0, np.ones(len(neg_idx), dtype=bool)]),
        )
    if y_eval.sum() == 0:
        return np.nan
    label = np.concatenate([np.ones(len(pos_dec_idx)), np.zeros(len(neg_idx))])
    return sklearn.metrics.precision_score(label, y_eval)


def get_metrics(
    metric_function,
    list_zs: Iterable,
    list_fdrs: Iterable,
    labels: Iterable[str],
    pos_ctrl_idx: Iterable,
    neg_ctrl_idx: Iterable,
    critical_value: float = 0.1,
    filter_sign: bool = False,
    **kwargs
) -> pd.Series:
    """For the list of z-scores and FDRS (paired) from different methods, calculate the metric with `metric_function`."""
    # nan_mask = np.isnan(list_zs[0])
    # for z in list_zs[1:]:
    #     nan_mask = np.logical_or(nan_mask, np.isnan(z))
    # for metrics in [list_zs, list_fdrs, pos_ctrl_idx, neg_ctrl_idx]:

    metrics = list(
        map(
            lambda z, fdr: metric_function(
                z,
                fdr,
                pos_ctrl_idx,
                neg_ctrl_idx,
                critical_value=critical_value,
                filter_sign=filter_sign,
                **kwargs
            ),
            list_zs,
            list_fdrs,
        )
    )
    df = pd.Series(
        metrics,
        index=labels,
    )
    return df


def get_precision(
    list_zs,
    list_fdrs_inc,
    list_fdrs_dec,
    z_labels,
    pos_inc_idx,
    pos_dec_idx,
    neg_ctrl_idx,
    prefix="",
    critical_value=0.1,
    filter_sign=False,
):
    prec_incs = get_metrics(
        _get_precision_inc,
        list_zs,
        list_fdrs_inc,
        z_labels,
        pos_inc_idx,
        neg_ctrl_idx,
        # prefix=prefix,
        critical_value=critical_value,
        filter_sign=filter_sign,
    )
    prec_decs = get_metrics(
        _get_precision_dec,
        list_zs,
        list_fdrs_dec,
        z_labels,
        pos_dec_idx,
        neg_ctrl_idx,
        # prefix=prefix,
        critical_value=critical_value,
        filter_sign=filter_sign,
    )
    df = pd.DataFrame(
        {"{}_inc".format(prefix): prec_incs, "{}_dec".format(prefix): prec_decs},
        index=z_labels,
    )
    return df


def _get_recall_inc(
    z, fdr, pos_inc_idx, neg_idx, critical_value=0.1, filter_sign=False
):
    y_eval = np.concatenate([fdr[pos_inc_idx], fdr[neg_idx]]) < critical_value
    if filter_sign:
        y_eval = np.logical_and(
            y_eval,
            np.concatenate([z[pos_inc_idx] > 0, np.ones(len(neg_idx), dtype=bool)]),
        )
    label = np.concatenate([np.ones(len(pos_inc_idx)), np.zeros(len(neg_idx))])
    return sklearn.metrics.recall_score(label, y_eval)


def _get_recall_dec(
    z, fdr, pos_dec_idx, neg_idx, critical_value=0.1, filter_sign=False
):
    y_eval = np.concatenate([fdr[pos_dec_idx], fdr[neg_idx]]) < critical_value
    if filter_sign:
        y_eval = np.logical_and(
            y_eval,
            np.concatenate([z[pos_dec_idx] < 0, np.ones(len(neg_idx), dtype=bool)]),
        )
    label = np.concatenate([np.ones(len(pos_dec_idx)), np.zeros(len(neg_idx))])
    return sklearn.metrics.recall_score(label, y_eval)


def get_recall(
    list_zs,
    list_fdrs_inc,
    list_fdrs_dec,
    z_labels,
    pos_inc_idx,
    pos_dec_idx,
    neg_ctrl_idx,
    prefix="",
    critical_value=0.1,
    filter_sign=False,
):
    recall_incs = list(
        map(
            lambda z, fdr: _get_recall_inc(
                z,
                fdr,
                pos_inc_idx,
                neg_ctrl_idx,
                critical_value=critical_value,
                filter_sign=filter_sign,
            ),
            list_zs,
            list_fdrs_inc,
        )
    )
    recall_decs = list(
        map(
            lambda z, fdr: _get_recall_dec(
                z,
                fdr,
                pos_dec_idx,
                neg_ctrl_idx,
                critical_value=critical_value,
                filter_sign=filter_sign,
            ),
            list_zs,
            list_fdrs_dec,
        )
    )
    df = pd.DataFrame(
        {"{}_inc".format(prefix): recall_incs, "{}_dec".format(prefix): recall_decs},
        index=z_labels,
    )
    return df


def _get_fdr_f_score_inc(
    z, fdr, pos_inc_idx, neg_idx, critical_value=0.1, filter_sign=False, beta=1
):
    y_eval = np.concatenate([fdr[pos_inc_idx], fdr[neg_idx]]) < critical_value
    if filter_sign:
        y_eval = np.logical_and(
            y_eval,
            np.concatenate([z[pos_inc_idx] > 0, np.ones(len(neg_idx), dtype=bool)]),
        )
    label = np.concatenate([np.ones(len(pos_inc_idx)), np.zeros(len(neg_idx))])
    return sklearn.metrics.fbeta_score(label, y_eval, beta=beta)


def _get_fdr_f_score_dec(
    z, fdr, pos_dec_idx, neg_idx, critical_value=0.1, filter_sign=False, beta=1
):
    y_eval = np.concatenate([fdr[pos_dec_idx], fdr[neg_idx]]) < critical_value
    if filter_sign:
        y_eval = np.logical_and(
            y_eval,
            np.concatenate([z[pos_dec_idx] < 0, np.ones(len(neg_idx), dtype=bool)]),
        )
    label = np.concatenate([np.ones(len(pos_dec_idx)), np.zeros(len(neg_idx))])
    return sklearn.metrics.fbeta_score(label, y_eval, beta=beta)


def get_fdr_f(
    list_zs,
    list_fdrs_inc,
    list_fdrs_dec,
    z_labels,
    pos_inc_idx,
    pos_dec_idx,
    neg_ctrl_idx,
    prefix="",
    critical_value=0.1,
    filter_sign=False,
    beta=1,
):
    f_incs = list(
        map(
            lambda z, fdr: _get_fdr_f_score_inc(
                z,
                fdr,
                pos_inc_idx,
                neg_ctrl_idx,
                critical_value=critical_value,
                filter_sign=filter_sign,
                beta=beta,
            ),
            list_zs,
            list_fdrs_inc,
        )
    )
    f_decs = list(
        map(
            lambda z, fdr: _get_fdr_f_score_dec(
                z,
                fdr,
                pos_dec_idx,
                neg_ctrl_idx,
                critical_value=critical_value,
                filter_sign=filter_sign,
                beta=beta,
            ),
            list_zs,
            list_fdrs_dec,
        )
    )
    df = pd.DataFrame(
        {"{}_inc".format(prefix): f_incs, "{}_dec".format(prefix): f_decs},
        index=z_labels,
    )
    return df


def plot_aucs_OvR(
    axes,
    z,
    pos_inc_idx,
    pos_dec_idx,
    neg_ctrl_idx,
    critical_value=1.96,
    plot_cv_legend=False,
    show_auroc=True,
    **kwargs
):
    auroc_inc = plot_auc_(
        axes[0],
        z,
        pos_inc_idx,
        neg_ctrl_idx,
        critical_value=critical_value,
        plot_cv_legend=plot_cv_legend,
        show_auroc=show_auroc,
        **kwargs
    )
    auroc_dec = plot_auc_(
        axes[1],
        -z,
        pos_dec_idx,
        neg_ctrl_idx,
        critical_value=critical_value,
        plot_cv_legend=plot_cv_legend,
        show_auroc=show_auroc,
        **kwargs
    )

    axes[0].legend(
        bbox_to_anchor=(0, -0.7, 1, 0.5),
        loc="upper center",
        mode="expand",
        title="model (AUROC)",
    )
    axes[1].legend(
        bbox_to_anchor=(0, -0.7, 1, 0.5),
        loc="upper center",
        mode="expand",
        title="model (AUROC)",
    )

    axes[0].set_title("PosCtrl_inc vs NegCtrl")
    axes[1].set_title("PosCtrl_dec vs NegCtrl")
    axes[0].set_xlabel("FPR")
    axes[0].set_ylabel("TPR")
    axes[1].set_xlabel("FPR")
    axes[1].set_ylabel("TPR")
    return (auroc_inc, auroc_dec)


def _get_auroc_inc(z, _, pos_inc_idx, neg_ctrl_idx, *args, **kwargs):
    y_eval = np.concatenate([z[pos_inc_idx], z[neg_ctrl_idx]])
    label = np.concatenate([np.ones(len(pos_inc_idx)), np.zeros(len(neg_ctrl_idx))])
    try:
        nan_mask = ~np.isnan(y_eval)
        auroc = sklearn.metrics.roc_auc_score(label[nan_mask], y_eval[nan_mask])
        auroc = sklearn.metrics.roc_auc_score(label, y_eval)
        return auroc
    except Exception as e:
        print(e)
        return np.nan


def _get_auroc_dec(z, _, pos_dec_idx, neg_ctrl_idx, *args, **kwargs):
    try:
        y_eval = np.concatenate([-z[pos_dec_idx], -z[neg_ctrl_idx]])
    except Exception as e:
        print(e)
    label = np.concatenate([np.ones(len(pos_dec_idx)), np.zeros(len(neg_ctrl_idx))])
    try:
        nan_mask = ~np.isnan(y_eval)
        auroc = sklearn.metrics.roc_auc_score(label[nan_mask], y_eval[nan_mask])
    except Exception as e:
        print(e)
        return np.nan
    return auroc


def get_auroc(
    list_zs, z_labels, pos_inc_idx, pos_dec_idx, neg_ctrl_idx, prefix="", **kwargs
):
    auroc_incs = list(
        map(lambda z: _get_auroc_inc(z, "_", pos_inc_idx, neg_ctrl_idx), list_zs)
    )
    auroc_decs = list(
        map(lambda z: _get_auroc_dec(z, "_", pos_dec_idx, neg_ctrl_idx), list_zs)
    )
    df = pd.DataFrame(
        {"{}_inc".format(prefix): auroc_incs, "{}_dec".format(prefix): auroc_decs},
        index=z_labels,
    )
    return df


def _get_auprc_inc(z, _, pos_inc_idx, neg_ctrl_idx, **kwargs):
    y_eval = np.concatenate([z[pos_inc_idx], z[neg_ctrl_idx]])
    label = np.concatenate([np.ones(len(pos_inc_idx)), np.zeros(len(neg_ctrl_idx))])
    nan_mask = ~np.isnan(y_eval)
    try:
        prec, recall, thres = sklearn.metrics.precision_recall_curve(
            label[nan_mask], y_eval[nan_mask]
        )
        auprc = sklearn.metrics.auc(recall, prec)
        return auprc
    except:
        return np.nan


def _get_auprc_dec(z, _, pos_dec_idx, neg_ctrl_idx, **kwargs):
    y_eval = np.concatenate([-z[pos_dec_idx], -z[neg_ctrl_idx]])
    label = np.concatenate([np.ones(len(pos_dec_idx)), np.zeros(len(neg_ctrl_idx))])
    nan_mask = ~np.isnan(y_eval)
    try:
        prec, recall, thres = sklearn.metrics.precision_recall_curve(
            label[nan_mask], y_eval[nan_mask]
        )
        auprc = sklearn.metrics.auc(recall, prec)
        return auprc
    except:
        return np.nan


def get_auprc(
    list_zs, z_labels, pos_inc_idx, pos_dec_idx, neg_ctrl_idx, prefix="", **kwargs
):
    auprc_incs = list(
        map(lambda z: _get_auprc_inc(z, "_", pos_inc_idx, neg_ctrl_idx), list_zs)
    )
    auprc_decs = list(
        map(lambda z: _get_auprc_dec(-z, "_", pos_dec_idx, neg_ctrl_idx), list_zs)
    )
    df = pd.DataFrame(
        {"{}_inc".format(prefix): auprc_incs, "{}_dec".format(prefix): auprc_decs},
        index=z_labels,
    )
    return df


def plot_prcs_OvR(
    axes,
    z,
    pos_inc_idx,
    pos_dec_idx,
    neg_ctrl_idx,
    critical_value=1.96,
    plot_cv_legend=False,
    show_auc=True,
    **kwargs
):
    auroc_inc = plot_prc_(
        axes[0],
        z,
        pos_inc_idx,
        neg_ctrl_idx,
        critical_value=critical_value,
        plot_cv_legend=plot_cv_legend,
        show_auc=show_auc,
        **kwargs
    )
    auroc_dec = plot_prc_(
        axes[1],
        -z,
        pos_dec_idx,
        neg_ctrl_idx,
        critical_value=-critical_value,
        plot_cv_legend=plot_cv_legend,
        show_auc=show_auc,
        **kwargs
    )

    axes[0].legend(
        bbox_to_anchor=(0, -0.7, 1, 0.5),
        loc="upper center",
        mode="expand",
        title="model (AUPRC)",
    )
    axes[1].legend(
        bbox_to_anchor=(0, -0.7, 1, 0.5),
        loc="upper center",
        mode="expand",
        title="model (AUPRC)",
    )

    axes[0].set_title("PosCtrl_inc vs NegCtrl")
    axes[1].set_title("PosCtrl_dec vs NegCtrl")
    axes[0].set_xlabel("Recall")
    axes[0].set_ylabel("Precision")
    axes[1].set_xlabel("Recall")
    axes[1].set_ylabel("Precision")
    return (auroc_inc, auroc_dec)
