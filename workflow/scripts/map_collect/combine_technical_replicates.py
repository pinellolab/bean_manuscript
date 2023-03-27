import sys
import os
import numpy as np
import scipy
import matplotlib.pyplot as plt
import bean as be


def main():
    lib = sys.argv[1]
    if lib not in ["LDLvar", "LDLRCDS"]:
        os.system(
            f"ln -s bean_count_{lib}.h5ad results/mapped/{lib}/bean_count_{lib}_combined.h5ad"
        )
        exit(0)
    bdata_path = f"results/mapped/{lib}/bean_count_{lib}.h5ad"
    bdata = be.read_h5ad(bdata_path)

    bdata.samples[["tech_rep", "bin"]] = bdata.samples.index.to_series().str.split(
        "_", expand=True
    )
    bdata.samples[["rep", "tech_idx"]] = bdata.samples["tech_rep"].str.split(
        ".", expand=True
    )

    rep_with_techreps = bdata.samples.loc[
        ~bdata.samples.tech_idx.isnull(), "rep"
    ].unique()
    bdata.samples.loc[
        bdata.samples.rep.isin(rep_with_techreps) & bdata.samples.tech_idx.isnull(),
        "tech_idx",
    ] = 0
    bdata.samples.loc[~bdata.samples.rep.isin(rep_with_techreps), "tech_idx"] = -1

    bdata.samples.loc[
        bdata.samples.rep.isin(rep_with_techreps), "tech_idx"
    ] = bdata.samples.loc[bdata.samples.rep.isin(rep_with_techreps), "tech_idx"].astype(
        int
    )

    bdata_no_tech_reps = bdata[:, ~bdata.samples.rep.isin(rep_with_techreps)]
    bdata_tech_rep_0 = bdata[:, bdata.samples.tech_idx == 0]
    bdata_tech_rep_1 = bdata[:, bdata.samples.tech_idx == 1]

    assert (
        bdata_tech_rep_0.samples.bin.values == bdata_tech_rep_1.samples.bin.values
    ).all(), (
        bdata_tech_rep_0.samples,
        bdata_tech_rep_1.samples,
    )
    assert (
        bdata_tech_rep_0.samples.rep.values == bdata_tech_rep_1.samples.rep.values
    ).all(), (
        bdata_tech_rep_0.samples,
        bdata_tech_rep_1.samples,
    )
    bdata_tech_rep_1.rename(bdata_tech_rep_0.samples.index, axis=1)

    fig_dir = f"results/mapped/{lib}/figs/"
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    plot_lfc_corr(bdata_tech_rep_0, bdata_tech_rep_1, fig_dir)

    bdata_added = bdata_tech_rep_0 + bdata_tech_rep_1
    bdata_techrep_combined = be.concat([bdata_added, bdata_no_tech_reps], axis=1)
    bdata_techrep_combined.samples = bdata_techrep_combined.samples.drop(
        ["tech_rep", "tech_idx"], axis=1
    )
    bdata_techrep_combined.write(f"results/mapped/{lib}/bean_count_{lib}_combined.h5ad")


def plot_lfc_corr(bdata_techrep_0, bdata_techrep_1, fig_dir):
    bdata_techrep_0.log_norm()
    bdata_techrep_1.log_norm()

    lfc0 = bdata_techrep_0.log_fold_change_reps(
        "bot", "top", rep_col="rep", compare_col="bin"
    )
    lfc1 = bdata_techrep_1.log_fold_change_reps(
        "bot", "top", rep_col="rep", compare_col="bin"
    )
    for rep in bdata_techrep_0.samples.rep.unique():
        lfc0_rep = lfc0[f"{rep}.bot_top.lfc"]
        lfc1_rep = lfc1[f"{rep}.bot_top.lfc"]
        x_valid = ~np.isnan(lfc0_rep)
        y_valid = ~np.isnan(lfc1_rep)
        valid = x_valid & y_valid
        plt.style.use("default")
        fig, ax = plt.subplots(figsize=(3, 3))
        ax.scatter(lfc0_rep[valid], lfc1_rep[valid])
        ax.set_xlabel("Technical replicate 0")
        ax.set_ylabel("Technical replicate 1")
        ax.set_title(
            f"bot/top LFC (R={scipy.stats.pearsonr(lfc0_rep, lfc1_rep)[0]:.3f})"
        )
        ax.set_aspect(1)
        fig.savefig(f"{fig_dir}/lfc_corr_{rep}.pdf")


if __name__ == "__main__":
    main()
