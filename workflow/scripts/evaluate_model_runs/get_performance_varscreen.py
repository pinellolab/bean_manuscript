"""Compare performances of multiple models for a screen object."""
from typing import Sequence, Tuple, List
import os
import argparse
from functools import partial
import numpy as np
import pandas as pd
import evaluate as ae

model_ids = ["Normal", "MixtureNormal"]
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
    parser.add_argument("--result-suffix", type=str, default="")
    return parser.parse_args()


def get_average_metric(df, pos_idx, neg_idx):
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
        tbl["neg_mu_z"] = -tbl["mu_z"]
        tbl = (
            tbl.set_index("target")
            .reindex(gene_order)
            .reset_index(drop=False)
            .add_suffix(f"_{model_id}")
        )

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


def get_cb2_results(
    cb2_prefix: str, gene_order: Sequence[str]
) -> Tuple[pd.DataFrame, List[str]]:
    res_path = f"{cb2_prefix}/CB2_gene.csv"
    tbl = pd.read_csv(res_path).set_index("gene").reindex(gene_order)
    tbl = tbl.add_suffix("_CB2")
    return tbl.reset_index(), ["CB2"]


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
    bean_result_path_format = f"results/model_runs/bean{args.result_suffix}/bean_run_result.{args.screen_name}/bean_element_result.{{}}.csv"
    guide_info_path = args.guide_info_path
    mageck_prefix = f"results/model_runs/mageck{args.result_suffix}/{args.screen_name}/"
    cb2_prefix = (
        f"results/model_runs/CB2{args.result_suffix}/CB2_run_result.{args.screen_name}/"
    )

    guide_info = pd.read_csv(guide_info_path).rename(
        columns={"Target gene/variant": "target_variant", "Group2": "target_group2"}
    )

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
    cb2_results, cb2_labels = get_cb2_results(cb2_prefix, gene_order)

    all_results = pd.concat([bean_results, mageck_results, cb2_results], axis=1)
    all_results.to_csv(
        f"results/model_runs/{args.screen_name}/all_scores{args.result_suffix}.csv"
    )

    target_info = (
        guide_info[["target", "target_variant", "target_group2"]]
        .drop_duplicates()
        .set_index("target", drop=True)
    ).reindex(gene_order)
    if args.control == "splicing":
        pos_idx = np.where(target_info["target_group2"] == "PosCtrl_inc")[0]
        neg_idx = np.where(target_info["target_group2"] == "PosCtrl_dec")[0]
    elif args.control == "strongest_splicing":
        pos_idx = np.where(target_info["target_variant"] == "MYLIP")[0]
        neg_idx = np.where(target_info["target_variant"] == "LDLR")[0]
    else:
        raise ValueError(
            f"Invalid --control {args.control} provided. Use one of 'splicing' or 'strongest_splicing'."
        )
    ctrl_idx = np.where(target_info["target_group2"] == "NegCtrl")[0]

    kwargs = {
        "res_tbl": all_results,
        "z_cols": get_cols("mu_z_{}", "sort_num|z_{}_{}", "pos|lfc_rra_{}")
        + ["logFC_CB2"],
        "fdr_inc_cols": get_cols(
            "neg_mu_z_{}",
            "sort_num|fdr_{}_{}",
            "neg|fdr_rra_{}",
        )
        + ["fdr_pa_CB2"],
        "fdr_dec_cols": get_cols(
            "mu_z_{}",
            "sort_num|fdr_{}_{}",
            "pos|fdr_rra_{}",
        )
        + ["fdr_pb_CB2"],
        "labels": model_ids + mageck_labels + ["CB2"],
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
        f"results/model_runs/{args.screen_name}/{args.control}{args.result_suffix}.metrics.xlsx",
        engine="xlsxwriter",
    )
    avg_perf_df.style.set_properties(**{"width": "px"}).background_gradient().to_excel(
        writer, sheet_name="average_metrics"
    )
    perf_df.style.set_properties(**{"width": "px"}).background_gradient().to_excel(
        writer, sheet_name="metrics"
    )
    writer.close()


if __name__ == "__main__":
    main()
