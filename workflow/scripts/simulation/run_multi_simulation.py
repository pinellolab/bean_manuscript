import os
import sys
import subprocess
import numpy as np
import pandas as pd
import argparse
from screen_simulation import SimulatedScreen
import bean as be
import bean.model.post_training.compare as bc

sys.path.append("../evaluate_model_runs/")
import evaluate as ae


def get_fs(plot_df, eval_key, posctrl_es=1.0, negctrl_es=0, critical_value=0.1):
    posctrl_idx = np.where(plot_df.effect_size == posctrl_es)[0]
    assert len(posctrl_idx) > 0, posctrl_es
    negctrl_idx = np.where(plot_df.effect_size == negctrl_es)[0]
    f1 = ae.get_fdr_f_score_inc(
        "",
        plot_df[eval_key],
        posctrl_idx,
        negctrl_idx,
        critical_value=critical_value,
        beta=1,
    )
    f01 = ae.get_fdr_f_score_inc(
        "",
        plot_df[eval_key],
        posctrl_idx,
        negctrl_idx,
        critical_value=critical_value,
        beta=0.1,
    )
    return f1, f01


def get_auprcs(
    plot_df,
    eval_key,
    posctrl_es=1.0,
    negctrl_es=0,
):
    posctrl_idx = np.where(plot_df.effect_size == posctrl_es)[0]
    assert len(posctrl_idx) > 0, posctrl_es
    negctrl_idx = np.where(plot_df.effect_size == negctrl_es)[0]
    return ae._get_auprc_inc(plot_df[eval_key], "", posctrl_idx, negctrl_idx)


def get_fs_ess(plot_df, eval_key, critical_value=0.1):
    ess = plot_df.effect_size.unique()[1:]
    f1s = []
    f01s = []
    for es in ess:
        f1, f01 = get_fs(eval_key, es, critical_value=critical_value)
        f1s.append(f1)
        f01s.append(f01)
    fs = pd.DataFrame({"posCtrls": ess, "F1": f1s, "F0.1": f01s})
    return fs


def get_auprcs_ess(plot_df, eval_key, critical_value=0.1):
    ess = plot_df.effect_size.unique()[1:]
    auprcs = []
    for es in ess:
        auprcs.append(get_auprcs(eval_key, es, critical_value=critical_value))
    auprcs = pd.DataFrame({"posCtrls": ess, "AUPRC": auprcs})
    return auprcs


def get_f_df(eval_keys, critical_value=0.1):
    dfs = []
    for i, k in enumerate(eval_keys):
        df = get_fs_ess(k, critical_value=critical_value)
        df.columns = [df.columns[0]] + [c + "({})".format(k) for c in df.columns[1:]]
        if i == 0:
            dfs.append(df)
        else:
            dfs.append(df.iloc[:, 1:])
    return pd.concat(dfs, axis=1)


def format_sim_screen(sim):
    sim.screen_res[0].guides = sim.screen_res[0].guides.rename(
        columns={"target_id": "target", "guide_id": "name"}
    )

    sim.screen_res[0].condit = sim.screen_res[0].condit.rename(
        columns={"replicate": "rep"}
    )
    if sim.sorting_mode == "topbot":
        sim.screen_res[0].condit["sorting_scheme"] = "topbot"
    else:
        sim.screen_res[0].condit["sorting_scheme"] = "sort"
    sim_obs = sim.screen_res[0][:, sim.screen_res[0].condit.sort != "mid"].copy()
    sim_obs.layers["X_bcmatch"] = sim_obs.layers["X_bcmatch"].astype(float)
    sim_obs.get_guide_edit_rate()
    return sim_obs


def make_simulation_file(kwargs):
    sim = SimulatedScreen.SimulatedScreen(**kwargs)
    sim.simulate_reps()
    sim_obs = format_sim_screen(sim)
    return sim_obs


def run_multi_sim_file(
    screen_prefix, basename, run_result_prefix=None, n_sim=1, **kwargs
):
    """Write multipe simulation objects and run BEAN + MAGeCK on them
    Run BEAN & MAGeCK to save analysis result. run_pyro_mageck_models.py is executed.
    """
    for i in range(n_sim):
        # make data
        if not os.path.exists(screen_prefix):
            os.makedirs(screen_prefix)
        screen_file_path = f"{screen_prefix}/{basename}.{i}.h5ad"
        if not os.path.exists(screen_file_path):
            sim_obs = make_simulation_file(kwargs)
            sim_obs.write(screen_file_path)

        # run model
        if run_result_prefix is None:
            run_result_prefix = screen_prefix
        bean_command_list = [
            "sh",
            "../run_models/run_bean_var.sh",
            screen_file_path,
        ]
        print(f"Running...\n{' '.join(bean_command_list)}")
        subprocess.run(bean_command_list, check=True)


def evaluate(screen_prefix, basename, run_result_prefix=None, n_sim=1):
    """ """
    if run_result_prefix is None:
        run_result_prefix = screen_prefix
    model_ids = ["A", "B16_0", "B16", "B16_1", "B14", "B1", "B18", "B19"]
    plot_dfs = []
    guide_dfs = []
    for i in range(n_sim):
        res_df = bc.get_multiple_result_table(
            adata_path=f"{run_result_prefix}/{basename}.{i}.h5ad",
            prefix=f"{run_result_prefix}/{basename}.{i}",
            model_ids=model_ids,
            path_pattern_string="{}.model{}.result.pkl",
        )
        # read mageck result
        subfolders = ["", "EM", "EM_fit"]
        suffixes = ["", "_EM", "_EMf"]
        mageck_prefixes = [
            f"{run_result_prefix}/{basename}.{i}/{subfolder}"
            for subfolder in subfolders
        ]
        mageck_results = [
            bc.read_mageck_result(mp, suffix=sf)
            for mp, sf in zip(mageck_prefixes, suffixes)
        ]

        plot_df = pd.concat([res_df] + mageck_results, axis=1)
        plot_df["sim_iter"] = i
        plot_dfs.append(plot_df)

        guide_df = be.read_h5ad(f"{run_result_prefix}/{basename}.{i}.h5ad").guides
        guide_df["sim_iter"] = i
        guide_dfs.append(guide_df)
    plot_df_collected = pd.concat(plot_dfs)
    out_path = f"{run_result_prefix}/eval_result.{basename}.collected.csv"
    plot_df_collected.to_csv(out_path)
    print(f"Result saved at {out_path}")

    guide_df_collected = pd.concat(guide_dfs)
    guides_out_path = f"{run_result_prefix}/guide_info.{basename}.collected.csv"
    guide_df_collected.to_csv(guides_out_path)
    print(f"Result saved at {guides_out_path}")


def validate_args(args):
    kwargs = {
        k: v
        for k, v in vars(args).items()
        if (not k in ["screen_prefix", "basename", "n_sim", "run_result_prefix"])
    }
    return kwargs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run model by creating multiple simulation data."
    )

    parser.add_argument(
        "screen_prefix", type=str, help="Path of an ReporterScreen object to write"
    )
    parser.add_argument(
        "basename", type=str, help="Basename of the ReporterScreen object"
    )
    parser.add_argument(
        "--n_sim", "-n", type=int, default=1, help="Number of screens to produce"
    )
    parser.add_argument(
        "--run_result_prefix",
        "-rp",
        type=str,
        default=None,
        help="Relative path to save run results",
    )
    parser.add_argument(
        "--n_reads_per_sample", "-d", help="# reads per sample", default=800000
    )
    parser.add_argument(
        "--sorting_mode",
        default="bins",
        help="Sorting mode. Either ['topbot' or 'bins'].",
    )
    parser.add_argument(
        "--scale_by_accessibility",
        "-sa",
        default=False,
        action="store_true",
        help="scale by accessibility",
    )
    parser.add_argument(
        "--scale_with_variability",
        "-sv",
        default=False,
        action="store_true",
        help="scale by accessibility with variability in the scaled outcome by observed residual.",
    )
    parser.add_argument(
        "--has_reporter",
        "-r",
        default=True,
        help="Screen has reporter readout",
    )

    args = parser.parse_args()
    kwargs = validate_args(args)
    run_multi_sim_file(
        args.screen_prefix, args.basename, args.run_result_prefix, args.n_sim, **kwargs
    )
    evaluate(args.screen_prefix, args.basename, args.run_result_prefix, args.n_sim)
