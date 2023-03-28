"""Compare performances of multiple models for a screen object."""
from typing import Sequence, Tuple, List
import os
import argparse
from functools import partial
import numpy as np
import pandas as pd
import evaluate as ae

model_ids = ["Normal", "MixtureNormal", "MixtureNormal+Acc"]
mageck_mle_disp_modes = ["sort", "sort_var"]
mageck_mle_pi_modes = ["", "EM", "EMf"]
mageck_rra_modes = ["bot"]


def parse_args():
    parser = argparse.ArgumentParser(description="Get performance metric and plot.")
    parser.add_argument(
        "screen_name",
        type=str,
    )
    parser.add_argument(
        "--guide-info-path",
        type=str,
        default="resources/gRNA_info/LDLvar_gRNA_bean_accessibility.csv",
    )
    parser.add_argument("--control", type=str, default="splicing")
    return parser.parse_args()


def get_average_metric(df, pos_idx, neg_idx):
    print(df)
    return (df.iloc[:, 0] * len(pos_idx) + df.iloc[:, 1] * len(neg_idx)) / (
        len(pos_idx) + len(neg_idx)
    )


def get_mean_metric(df):
    splicing = df.iloc[:, :2].mean(axis=1)
    strong = df.iloc[:, 2:].mean(axis=1)
    return pd.DataFrame({"splicing": splicing, "mylip_ldlr": strong})


def get_bean_results(result_path_format):
    gene_order = None
    res_tbls = []
    for model_id in model_ids:
        tbl = pd.read_csv(result_path_format.format(model_id))
        if gene_order is None:
            gene_order = tbl["target"]
        tbl = tbl.add_suffix(f"_{model_id}")
        res_tbls.append(tbl.reset_index())
    return pd.concat(res_tbls, axis=1), gene_order


def get_mageck_results(mageck_prefix: str, gene_order: Sequence[str]):
    def read_and_add_suffix(disp_mode, pi_mode, gene_order):
        res_path = f"{mageck_prefix}/{pi_mode}/{disp_mode}.gene_summary.txt"
        tbl = pd.read_table(res_path, sep="\t", index_col=0).reindex(gene_order)
        tbl = tbl.add_suffix(f"_{disp_mode}_{pi_mode}")
        return tbl.reset_index()

    def read_and_add_suffix_rra(mode, gene_order):
        res_path = f"{mageck_prefix}/rra_{mode}.gene_summary.txt"
        tbl = pd.read_table(res_path, sep="\t", index_col=0).reindex(gene_order)
        tbl = tbl.add_suffix(f"_rra_{mode}")
        return tbl.reset_index()

    tbls = [
        read_and_add_suffix(disp_mode, pi_mode, gene_order)
        for disp_mode in mageck_mle_disp_modes
        for pi_mode in mageck_mle_pi_modes
    ]
    labels = [
        f"MAGeCK-MLE_{disp_mode}_{pi_mode}"
        for disp_mode in mageck_mle_disp_modes
        for pi_mode in mageck_mle_pi_modes
    ]

    # RRA
    tbls += [read_and_add_suffix_rra(mode, gene_order) for mode in mageck_rra_modes]
    labels += [f"MAGeCK-RRA_{mode}" for mode in mageck_rra_modes]
    return pd.concat(tbls, axis=1), labels


def get_other_results(
    other_prefix: str, gene_order: Sequence[str]
) -> Tuple[pd.DataFrame, List[str]]:
    NotImplemented


def get_metrics_bidirectional(
    metric_fn,
    res_tbl,
    z_cols,
    fdr_inc_cols,
    fdr_dec_cols,
    labels,
    pos_idx,
    neg_idx,
    ctrl_idx,
    control_label="splicing",
):
    return metric_fn(
        list_zs=[res_tbl[c] for c in z_cols],
        list_fdrs_inc=[res_tbl[c] for c in fdr_inc_cols],
        list_fdrs_dec=[res_tbl[c] for c in fdr_dec_cols],
        z_labels=labels,
        pos_inc_idx=pos_idx,
        pos_dec_idx=neg_idx,
        neg_ctrl_idx=ctrl_idx,
        prefix=control_label,
        # filter_sign=True
    ).round(3)


def get_metrics(
    metric_fn,
    res_tbl,
    z_cols,
    fdr_cols,
    labels,
    pos_idx,
    ctrl_idx,
    control_label="splicing",
):
    return metric_fn(
        list_zs=[res_tbl[c] for c in z_cols],
        list_fdrs=[res_tbl[c] for c in fdr_cols],
        z_labels=labels,
        pos_idx=pos_idx,
        ctrl_idx=ctrl_idx,
        prefix=control_label,
        # filter_sign=True
    ).round(3)


def get_auroc_auprc(
    metric_fn,
    labels,
    res_tbl,
    z_cols,
    pos_idx,
    neg_idx,
    ctrl_idx,
    control_label="splicing",
    **kwargs,
):
    return metric_fn(
        [res_tbl[c] for c in z_cols], labels, pos_idx, neg_idx, ctrl_idx, control_label
    ).round(3)


def main():
    args = parse_args()
    output_path = f"results/model_runs/{args.screen_name}"
    os.makedirs(output_path, exist_ok=True)
    bean_result_path_format = f"results/model_runs/bean/bean_run_result.{args.screen_name}/bean_element_result.{{}}.csv"
    guide_info_path = args.guide_info_path
    mageck_prefix = f"results/model_runs/mageck/{args.screen_name}/"
    other_prefix = NotImplemented

    guide_info = pd.read_csv(guide_info_path)
    target_info = (
        guide_info[["target", "Target gene/variant", "Group2"]]
        .drop_duplicates()
        .set_index("target", drop=True)
    )

    if args.control == "splicing":
        pos_idx = np.where(target_info.Group2 == "PosCtrl_inc")[0]
        neg_idx = np.where(target_info.Group2 == "PosCtrl_dec")[0]
    elif args.control == "strongest_splicing":
        pos_idx = np.where(target_info["Target gene/variant"] == "MYLIP")[0]
        neg_idx = np.where(target_info["Target gene/variant"] == "LDLR")[0]
    else:
        raise ValueError(
            f"Invalid --control {args.control} provided. Use one of 'splicing' or 'strongest_splicing'."
        )
    ctrl_idx = np.where(target_info.Group2 == "NegCtrl")[0]

    def get_cols(bean_pattern, mageck_pattern, mageck_rra_pattern):
        return (
            [bean_pattern.format(m) for m in model_ids]
            + [
                mageck_pattern.format(disp_mode, pi_mode)
                for disp_mode in mageck_mle_disp_modes
                for pi_mode in mageck_mle_pi_modes
            ]
            + [mageck_rra_pattern.format(mode) for mode in mageck_rra_modes]
        )

    bean_results, gene_order = get_bean_results(bean_result_path_format)
    mageck_results, mageck_labels = get_mageck_results(mageck_prefix, gene_order)
    other_results, other_labels = get_other_results(other_prefix, gene_order)
    """This is where other method output should be read in. If the output does not have the same order as they appear in the input rows, use `gene_order` to sort them."""

    all_results = pd.concat([bean_results, mageck_results], axis=1)
    all_results.to_csv(f"results/model_runs/{args.screen_name}/all_scores.csv")

    kwargs = {
        "res_tbl": pd.concat([bean_results, mageck_results], axis=1),
        "z_cols": get_cols("mu_z_{}", "sort_num|z_{}_{}", "pos|lfc_rra_{}"),
        "fdr_inc_cols": get_cols("fdr_inc_{}", "sort_num|fdr_{}_{}", "neg|fdr_rra_{}"),
        "fdr_dec_cols": get_cols("fdr_dec_{}", "sort_num|fdr_{}_{}", "pos|fdr_rra_{}"),
        "labels": model_ids + mageck_labels,
        "pos_idx": pos_idx,
        "neg_idx": neg_idx,
        "ctrl_idx": ctrl_idx,
    }

    precisions = get_metrics_bidirectional(ae.get_precision, **kwargs)
    recalls = get_metrics_bidirectional(ae.get_recall, **kwargs)
    fdr_f1s = get_metrics_bidirectional(partial(ae.get_fdr_f, beta=1), **kwargs)
    fdr_f01s = get_metrics_bidirectional(partial(ae.get_fdr_f, beta=0.1), **kwargs)
    aurocs = get_auroc_auprc(ae.get_auroc, **kwargs)
    auprcs = get_auroc_auprc(ae.get_auprc, **kwargs)

    perf_dict = {
        "Precision": precisions,
        "Recall": recalls,
        "F1": fdr_f1s,
        "F0.1": fdr_f01s,
        "AUROC": aurocs,
        "AUPRC": auprcs,
    }

    avg_perf_dict = {
        k: get_average_metric(v, pos_idx, neg_idx) for k, v in perf_dict.items()
    }
    perf_df = pd.concat(perf_dict.values(), axis=1, keys=perf_dict.keys())
    avg_perf_df = pd.concat(avg_perf_dict.values(), axis=1, keys=avg_perf_dict.keys())

    writer = pd.ExcelWriter(
        f"results/model_runs/{args.screen_name}/{args.control}.metrics.xlsx",
        engine="xlsxwriter",
    )
    avg_perf_df.style.set_properties(**{"width": "px"}).background_gradient().to_excel(
        writer, sheet_name="average_metrics"
    )
    perf_df.style.set_properties(**{"width": "px"}).background_gradient().to_excel(
        writer, sheet_name="metrics"
    )
    writer.save()


if __name__ == "__main__":
    main()
