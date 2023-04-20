"""Compare performances of multiple models for a screen object."""
from typing import Sequence, Tuple, List
import os
import argparse
from functools import partial
import numpy as np
import pandas as pd
import evaluate as ae

model_ids = [
    "MultiMixtureNormal",
    "MultiMixtureNormal+Acc",
    "Normal_allEdited",
    "Normal_behive",
]
mageck_mle_disp_modes = ["sort", "sort_var"]
mageck_rra_modes = ["bot"]


def parse_args():
    parser = argparse.ArgumentParser(description="Get performance metric and plot.")
    parser.add_argument(
        "screen_name",
        type=str,
    )
    parser.add_argument("--result-suffix", type=str, default="")
    return parser.parse_args()


def get_average_metric(df, pos_idx, neg_idx):
    return (df.iloc[:, 0] * len(pos_idx) + df.iloc[:, 1] * len(neg_idx)) / (
        len(pos_idx) + len(neg_idx)
    )


def get_bean_results(result_path_format):
    res_tbl = None
    for model_id in model_ids:
        tbl = pd.read_csv(result_path_format.format(model_id))
        if "edit" not in tbl.columns:
            cols = tbl.columns.tolist()
            if cols[1].startswith("target_"):
                cols[1] = "edit"
            elif cols[0].startswith("target_"):
                cols[0] = "edit"
            else:
                raise ValueError(
                    f"Don't know the target element column in {result_path_format.format(model_id)}."
                )
            tbl.columns = cols
        tbl = tbl.add_suffix(f"_{model_id}")
        if res_tbl is None:
            res_tbl = tbl
            res_tbl["edit"] = tbl[f"edit_{model_id}"]
        else:
            res_tbl = res_tbl.merge(
                tbl, left_on="edit", right_on=f"edit_{model_id}", how="outer"
            )
            print(res_tbl)

    res_tbl = res_tbl.rename(
        columns={
            "edit": "variant",
            f"group_{model_ids[0]}": "group",
        }
    )
    return res_tbl


def get_mageck_results(mageck_prefix):
    gene_order = None

    def read_and_add_suffix(disp_mode):
        nonlocal gene_order
        res_path = f"{mageck_prefix}/{disp_mode}.gene_summary.txt"
        tbl = pd.read_table(res_path, sep="\t")
        if gene_order is None:
            gene_order = tbl.Gene.tolist()
        tbl = tbl.set_index("Gene").reindex(gene_order).reset_index()
        tbl = tbl.add_suffix(f"_{disp_mode}")
        return tbl

    def read_and_add_suffix_rra(mode, gene_order):
        res_path = f"{mageck_prefix}/rra_{mode}.gene_summary.txt"
        tbl = pd.read_table(res_path, sep="\t", index_col=0).reindex(gene_order)
        tbl["fdr"] = tbl[["neg|fdr", "pos|fdr"]].min(axis=1)
        tbl = tbl.add_suffix(f"_rra_{mode}")
        return tbl.reset_index()

    tbls = [read_and_add_suffix(disp_mode) for disp_mode in mageck_mle_disp_modes]
    labels = [f"MAGeCK-MLE_{disp_mode}" for disp_mode in mageck_mle_disp_modes]
    tbls += [read_and_add_suffix_rra(mode, gene_order) for mode in mageck_rra_modes]
    labels += [f"MAGeCK-RRA_{mode}" for mode in mageck_rra_modes]

    res_tbl = pd.concat(tbls, axis=1)
    res_tbl = res_tbl.rename(columns={f"Gene_{mageck_mle_disp_modes[0]}": "variant"})
    return res_tbl, labels


def get_cb2_results(cb2_prefix: str) -> Tuple[pd.DataFrame, List[str]]:
    res_path = f"{cb2_prefix}/CB2_with_bcmatch_gene.csv"
    tbl = pd.read_csv(res_path).rename(columns={"gene": "variant"})
    tbl = tbl.set_index("variant").add_suffix("_CB2").reset_index(drop=False)
    return tbl.reset_index(), ["CB2"]


def get_crisphiermix_results(
    crisphiermix_prefix: str,
) -> Tuple[pd.DataFrame, List[str]]:
    res_path = f"{crisphiermix_prefix}/CRISPhieRmix_with_bcmatch.csv"
    tbl = pd.read_csv(res_path).rename(columns={"gene": "variant"})
    tbl = tbl.set_index("variant").add_suffix("_CRISPhieRmix").reset_index(drop=False)
    return tbl.reset_index(), ["CRISPhieRmix"]


def get_metrics(
    metric_fn, res_tbl, z_cols, fdr_cols, labels, pos_idx, ctrl_idx, **kwargs
):
    try:
        [res_tbl[c] for c in z_cols]
    except KeyError as e:
        print(e)
        print(res_tbl.columns.tolist())

    return ae.get_metrics(
        metric_fn,
        list_zs=[res_tbl[c] for c in z_cols],
        list_fdrs=[res_tbl[c] for c in fdr_cols],
        labels=labels,
        pos_ctrl_idx=pos_idx,
        neg_ctrl_idx=ctrl_idx,
        # prefix=control_label,
        # filter_sign=True
        **kwargs,
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
    mageck_suffixes=["_allEdited", "_behive"],
):
    pos_idx = np.where(cmp_results[group_col] == pos_group)[0]
    ctrl_idx = np.where(cmp_results[group_col] == ctrl_group)[0]

    def get_cols(
        bean_pattern,
        mageck_pattern,
        mageck_rra_pattern,
        mageck_suffixes=mageck_suffixes,
    ):
        bean_cols = [bean_pattern.format(m) for m in model_ids]
        mageck_cols = [
            mageck_pattern.format(disp_mode) for disp_mode in mageck_mle_disp_modes
        ] + [mageck_rra_pattern.format(mode) for mode in mageck_rra_modes]
        cols = bean_cols
        for mageck_suffix in mageck_suffixes:
            cols += [m + mageck_suffix for m in mageck_cols]
        return cols

    kwargs = {
        "res_tbl": cmp_results,
        "z_cols": get_cols(
            "mu_z_{}",
            "sort_num|z_{}",
            "pos|lfc_rra_{}",
        )
        + [f"logFC_CB2{sfx}" for sfx in mageck_suffixes]
        + [
            f"score_CRISPhieRmix{sfx}" for sfx in mageck_suffixes
        ],  # CRISPhieRmix doesn't have z-score output
        "fdr_cols": get_cols("fdr_dec_{}", "sort_num|fdr_{}", "pos|fdr_rra_{}")
        + [f"fdr_ts_CB2{sfx}" for sfx in mageck_suffixes]
        + [f"logfdr_CRISPhieRmix{sfx}" for sfx in mageck_suffixes],
        "labels": model_ids
        + mageck_labels
        + [f"CB2{sfx}" for sfx in mageck_suffixes]
        + [f"CRISPhieRmix{sfx}" for sfx in mageck_suffixes],
        "pos_idx": pos_idx,
        "ctrl_idx": ctrl_idx,
        "filter_sign": False,
    }
    try:
        precisions = get_metrics(ae._get_precision_dec, **kwargs)
    except:
        print(cmp_results.columns)
    recalls = get_metrics(ae._get_recall_dec, **kwargs)
    fdr_f1s = get_metrics(partial(ae._get_fdr_f_score_inc, beta=1), **kwargs)
    fdr_f01s = get_metrics(partial(ae._get_fdr_f_score_inc, beta=0.1), **kwargs)
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
    if not isinstance(s, str):
        return s
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
    bean_result_path_format = f"results/model_runs/bean{args.result_suffix}/bean_run_result.{args.screen_name}/bean_element_result.{{}}.csv"
    bean_results = get_bean_results(bean_result_path_format)
    all_results = None
    all_mageck_labels = []
    # evaluate performance with splicing variants as positive controls, for common variant of pair of methods.
    for annot_type in ["allEdited", "behive"]:
        mageck_prefix = f"results/model_runs/mageck{args.result_suffix}/{args.screen_name}.target_{annot_type}/"
        mageck_results, mageck_labels = get_mageck_results(mageck_prefix)
        mageck_labels = [f"{m}_{annot_type}" for m in mageck_labels]
        all_mageck_labels += mageck_labels
        cb2_prefix = f"results/model_runs/CB2/CB2_run_result.{args.screen_name}.target_{annot_type}/"
        cb2_prefix = f"results/model_runs/CRISPhieRmix/CRISPhieRmix_run_result.{args.screen_name}.target_{annot_type}/"
        cb2_results, cb2_labels = get_cb2_results(cb2_prefix)
        crisphiermix_results, crisphiermix_labels = get_crisphiermix_results(
            crisphiermix_prefix
        )
        try:
            cmp_results = (
                bean_results.merge(
                    mageck_results,
                    on="variant",
                    how="inner",
                )
                .merge(cb2_results, on="variant", how="inner")
                .merge(crisphiermix_results, on="variant", how="inner")
            )
        except:
            print(
                bean_results.merge(
                    mageck_results,
                    on="variant",
                    how="inner",
                ).columns
            )
            print(cb2_results.columns)
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
            .merge(
                cb2_results.set_index("variant")
                .add_suffix(f"_{annot_type}")
                .reset_index(),
                on="variant",
                how="outer",
            )
            .merge(
                crisphiermix_results.set_index("variant")
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
            .merge(
                cb2_results.set_index("variant")
                .add_suffix(f"_{annot_type}")
                .reset_index(),
                on="variant",
                how="outer",
            )
            .merge(
                crisphiermix_results.set_index("variant")
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
                mageck_suffixes=[""],
            )
            perf_df.to_csv(
                f"results/model_runs/{args.screen_name}/{annot_type}_{control}{args.result_suffix}.metrics.csv"
            )

    # evaluate performance by comparing to Clivar annotations
    all_results = ae.get_dms_df(
        all_results.reset_index(),
        join_how="outer",
        edit_str_column="variant",
        filter_possible_base_transitions=(("C", "T"), ("G", "A"))
        if "CBE" in args.screen_name
        else (
            ("A", "G"),
            ("T", "C"),
        ),
    )
    all_results.to_csv(
        f"results/model_runs/{args.screen_name}/all_scores{args.result_suffix}.csv"
    )
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
    perf_df_p.to_csv(
        f"results/model_runs/{args.screen_name}/clinvar_pb{args.result_suffix}.metrics.csv"
    )
    perf_df_plp.to_csv(
        f"results/model_runs/{args.screen_name}/clinvar_plpb{args.result_suffix}.metrics.csv"
    )


if __name__ == "__main__":
    main()
