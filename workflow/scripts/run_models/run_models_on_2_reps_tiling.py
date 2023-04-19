import os
from itertools import combinations
import bean as be
import subprocess
import argparse


def get_whole_reps(bdata):
    good_samples = bdata.samples.groupby("rep")["mask"].sum()
    whole_reps = good_samples[good_samples == 5].index.tolist()
    print(f"Using reps {whole_reps} for 2-replicate sampling...")
    return whole_reps


def get_string_combination(rep_list):
    combs = combinations(rep_list, 2)
    return [f"{c[0]},{c[1]}" for c in combs]


def run_models(bdata_paths, provide_negctrl=False):
    procs = []
    for bdata_path in bdata_paths:
        p = subprocess.Popen(
            [
                "sh",
                f"scripts/run_models/run_bean_tiling{'_negctrl' if provide_negctrl else ''}.sh",
                bdata_path,
            ]
        )
        procs.append(p)

        mageck_command = [
            "sh",
            "scripts/run_models/run_mageck_tiling.sh",
            bdata_path,
            f"results/model_runs/mageck{'_negctrl' if provide_negctrl else ''}",
        ]
        if provide_negctrl:
            mageck_command += [
                "--control-sgrna",
                f"resources/gRNA_info/LDLRCDS_{'CBE' if 'CBE' in bdata_path else 'ABE'}_negctrl_gRNA.txt",
            ]
        p = subprocess.Popen(mageck_command)
        procs.append(p)

        if provide_negctrl:
            p = subprocess.Popen(
                [
                    "sh",
                    "scripts/run_models/run_CRISPhieRmix_tiling_negctrl.sh",
                    f"{'CBE' if 'CBE' in bdata_path else 'ABE'} control",
                ]
            )
        p = subprocess.Popen(
            [
                "sh",
                "scripts/run_models/run_CB2_tiling.sh",
                bdata_path,
            ]
        )
        procs.append(p)
    for p in procs:
        p.wait()


def evaluate_model_runs(bdata_paths, used_negctrl=False):
    procs = []
    for bdata_path in bdata_paths:
        #        pass
        p = subprocess.Popen(
            [
                "python",
                "scripts/evaluate_model_runs/get_performance_tiling.py",
                os.path.basename(bdata_path),
            ]
            + (["--result-suffix", "_negctrl"] if used_negctrl else [])
        )
        procs.append(p)
    for p in procs:
        p.wait()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("bdata_path", type=str, help="Path of an ReporterScreen object")
    parser.add_argument("--use-negctrl", action="store_true", default=False)
    return parser.parse_args()


def main():
    args = parse_args()
    bdata_path = args.bdata_path
    bdata = be.read_h5ad(bdata_path)
    whole_reps = get_whole_reps(bdata)
    rep_combs = get_string_combination(whole_reps)
    print(f"Running for following combinations: {rep_combs}")
    sub_bdata_paths = []
    sub_bdata_prefixes = []
    procs = []
    for rep_comb in rep_combs:
        sub_bdata_prefix = (
            f"{os.path.dirname(bdata_path)}/"
            + os.path.basename(bdata_path).rsplit(".h5ad", 1)[0]
            + f"_{rep_comb.replace(',', '_')}"
        )
        sub_bdata_prefixes.append(sub_bdata_prefix)
        sub_bdata_paths.append(f"{sub_bdata_prefix}.h5ad")
        if not os.path.exists(f"{sub_bdata_prefix}.h5ad"):
            p = subprocess.Popen(
                [
                    "python",
                    "scripts/run_models/subset_screen.py",
                    bdata_path,
                    rep_comb,
                    f"{sub_bdata_prefix}.h5ad",
                ]
            )
            procs.append(p)
    for p in procs:
        p.wait()

    run_models(sub_bdata_paths, args.use_negctrl)
    evaluate_model_runs(sub_bdata_prefixes, args.use_negctrl)


if __name__ == "__main__":
    main()
