"""Compare performances of multiple models for a screen object."""
import os
import argparse
from functools import partial
import numpy as np
import pandas as pd
import evaluate as ae

model_ids = ["MultiMixtureNormal", "MultiMixtureNormal+Acc"]
mageck_mle_disp_modes = ["sort", "sort_var"]
mageck_rra_modes = ["bot"]

def parse_args():
    parser = argparse.ArgumentParser(description="Get performance metric and plot.")
    parser.add_argument(
        "screen_name",
        type=str,
    )
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
    res_tbls = []
    for model_id in model_ids:
        tbl = pd.read_csv(result_path_format.format(model_id))
        tbl = tbl.add_suffix(f"_{model_id}")
        res_tbls.append(tbl.reset_index())
    res_tbl = pd.concat(res_tbls, axis=1)
    res_tbl = res_tbl.rename(
        columns={f"edit_{model_ids[0]}": "variant", f"group_{model_ids[0]}": "group"}
    )
    return res_tbl


def get_mageck_results(mageck_prefix):
    def read_and_add_suffix(disp_mode):
        res_path = f"{mageck_prefix}/{disp_mode}.gene_summary.txt"
        tbl = pd.read_table(res_path, sep="\t")
        tbl = tbl.add_suffix(f"_{disp_mode}")
        return tbl

    tbls = [read_and_add_suffix(disp_mode) for disp_mode in mageck_mle_disp_modes]
    labels = [f"MAGeCK-MLE_{disp_mode}" for disp_mode in mageck_mle_disp_modes]
    res_tbl = pd.concat(tbls, axis=1)
    res_tbl = res_tbl.rename(columns={f"Gene_{mageck_mle_disp_modes[0]}": "variant"})
    return res_tbl, labels


def get_metrics(
    metric_fn,
    res_tbl,
    z_cols,
    fdr_cols,
    labels,
    pos_idx,
    ctrl_idx,
):
    return ae.get_metrics(
        metric_fn,
        list_zs=[res_tbl[c] for c in z_cols],
        list_fdrs=[res_tbl[c] for c in fdr_cols],
        labels=labels,
        pos_ctrl_idx=pos_idx,
        neg_ctrl_idx=ctrl_idx,
        # prefix=control_label,
        # filter_sign=True
    ).round(3)


def get_auroc_auprc(
    metric_fn,
    labels,
    res_tbl,
    z_cols,
    pos_idx,
    ctrl_idx,
    **kwargs,
):
    return ae.get_metrics(
        metric_fn,
        list_zs=[res_tbl[c] for c in z_cols],
        list_fdrs=[None for c in z_cols],
        labels=labels,
        pos_ctrl_idx=pos_idx,
        neg_ctrl_idx=ctrl_idx,
        # prefix=control_label,
        # filter_sign=True
    ).round(3)
    # return metric_fn(
    #     [res_tbl[c] for c in z_cols],
    #     labels,
    #     pos_idx,
    #     ctrl_idx,
    # ).round(3)


def get_performances(
    cmp_results,
    group_col="group",
    pos_group="splicing",
    ctrl_group="negctrl",
    mageck_labels=None,
    mageck_suffixes=None,
):
    pos_idx = np.where(cmp_results[group_col] == pos_group)[0]
    ctrl_idx = np.where(cmp_results[group_col] == ctrl_group)[0]

    def get_cols(bean_pattern, mageck_pattern, suffixes=None):
        bean_cols = [bean_pattern.format(m) for m in model_ids]
        if mageck_suffixes is None:
            mageck_cols = [
                mageck_pattern.format(disp_mode) for disp_mode in mageck_mle_disp_modes
            ]
        else:
            mageck_cols = [
                mageck_pattern.format(disp_mode) + sfx
                for disp_mode in mageck_mle_disp_modes
                for sfx in mageck_suffixes
            ]
            print(mageck_cols)
        return bean_cols + mageck_cols

    kwargs = {
        "res_tbl": cmp_results,
        "z_cols": get_cols("mu_z_{}", "sort_num|z_{}", suffixes=mageck_suffixes),
        "fdr_cols": get_cols(
            "fdr_dec_{}", "sort_num|wald-fdr_{}", suffixes=mageck_suffixes
        ),
        "labels": model_ids + mageck_labels,
        "pos_idx": pos_idx,
        "ctrl_idx": ctrl_idx,
    }

    precisions = get_metrics(ae._get_precision_dec, **kwargs)
    recalls = get_metrics(ae._get_recall_dec, **kwargs)
    fdr_f1s = get_metrics(partial(ae._get_fdr_f_score_dec, beta=1), **kwargs)
    fdr_f01s = get_metrics(partial(ae._get_fdr_f_score_dec, beta=0.1), **kwargs)
    aurocs = get_auroc_auprc(ae._get_auroc_dec, **kwargs)
    auprcs = get_auroc_auprc(ae._get_auprc_dec, **kwargs)

    perf_dict = {
        "Precision": precisions,
        "Recall": recalls,
        "F1": fdr_f1s,
        "F0.1": fdr_f01s,
        "AUROC": aurocs,
        "AUPRC": auprcs,
    }

    perf_df = pd.concat(perf_dict.values(), axis=1, keys=perf_dict.keys())
    return perf_df


def _convert_variant_to_mageck_format(s, guide_len_series: pd.Series):
    """Convert bean outcome for control to mageck format"""
    if "!" not in s:
        return s
    uid, edit = s.split("!")
    guide_len = guide_len_series[uid]
    edit_pos = int(edit.split(":")[0]) - (32 - 6 - guide_len)
    return f"{uid}_{edit_pos}:{edit.split(':')[-1]}"


def main():
    args = parse_args()
    output_path = f"results/model_runs/{args.screen_name}"
    os.makedirs(output_path, exist_ok=True)
    guide_info = pd.read_csv("resources/gRNA_info/LDLRCDS_gRNA_bean.csv", index_col=0)
    bean_result_path_format = f"results/model_runs/bean/bean_run_result.{args.screen_name}/bean_element_result.{{}}.csv"
    bean_results = get_bean_results(bean_result_path_format)
    all_results = None
    all_mageck_labels = []
    # evaluate performance with splicing variants as positive controls, for common variant of pair of methods.
    for annot_type in ["allEdited", "behive"]:
        mageck_prefix = (
            f"results/model_runs/mageck/{args.screen_name}.target_{annot_type}/"
        )
        mageck_results, mageck_labels = get_mageck_results(mageck_prefix)
        mageck_labels = [f"{m}_{annot_type}" for m in mageck_labels]
        all_mageck_labels += mageck_labels
        cmp_results = bean_results.merge(
            mageck_results,
            on="variant",
            how="inner",
        )
        bean_results_mageck_converted = bean_results
        bean_results_mageck_converted.variant = (
            bean_results_mageck_converted.variant.map(
                lambda v: _convert_variant_to_mageck_format(v, guide_info.guide_len)
            )
        )
        all_results = (
            bean_results_mageck_converted.merge(
                mageck_results.set_index("variant")
                .add_suffix(f"_{annot_type}")
                .reset_index(),
                on="variant",
                how="outer",
            )
            if all_results is None
            else all_results.merge(
                mageck_results.set_index("variant")
                .add_suffix(f"_{annot_type}")
                .reset_index(),
                on="variant",
                how="outer",
            )
        )
        for control in ["negctrl", "syn"]:
            perf_df = get_performances(
                cmp_results,
                pos_group="splicing",
                ctrl_group=control,
                mageck_labels=mageck_labels,
            )
            perf_df.to_csv(
                f"results/model_runs/{args.screen_name}/{annot_type}_{control}.metrics.csv"
            )

    # evaluate performance by comparing to Clivar annotations
    all_results = ae.get_dms_df(
        all_results.reset_index(), join_how="left", edit_str_column="variant"
    )
    all_results.to_csv(f"results/model_runs/{args.screen_name}/all_scores.csv")
    perf_df_p = get_performances(
        all_results,
        group_col="clinvar_annot_3",
        pos_group="Pathogenic",
        ctrl_group="Benign/Likely_Benign",
        mageck_labels=all_mageck_labels,
        mageck_suffixes=["_allEdited", "_behive"],
    )
    perf_df_plp = get_performances(
        all_results,
        group_col="clinvar_annot_2",
        pos_group="Pathogenic/Likely_Pathogenic",
        ctrl_group="Benign/Likely_Benign",
        mageck_labels=all_mageck_labels,
        mageck_suffixes=["_allEdited", "_behive"],
    )
    perf_df_p.to_csv(f"results/model_runs/{args.screen_name}/clinvar_pb.metrics.csv")
    perf_df_plp.to_csv(
        f"results/model_runs/{args.screen_name}/clinvar_plpb.metrics.csv"
    )


if __name__ == "__main__":
    main()
