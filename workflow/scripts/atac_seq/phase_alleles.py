import os
import numpy as np
import pandas as pd
import argparse
from bean import Allele, Edit
from tqdm.auto import tqdm
from patsy import dmatrices
import statsmodels.api as sm

edit_dict = {"A": "G", "T": "C"}


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("variant_id")
    return parser.parse_args()


def get_variant_info(var_id, amplicon_start, amplicon_seq):
    try:
        phased_snps = pd.read_table(
            f"results/atac_seq/{var_id}_phasedSNPs.vcf", comment="#", header=None
        )
    except pd.errors.EmptyDataError:
        raise ValueError(f"ERROR: Cannot phase variant {var_id} (No Snp).")
    phased_snps.columns = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "HepG2",
    ]
    phased_snps["GT"] = phased_snps.HepG2.map(lambda s: s.split(":")[0])
    phased_snps = phased_snps.loc[
        phased_snps.GT.isin(["1|0", "0|1"]) & (phased_snps.FILTER == "PASS")
    ]

    if len(phased_snps) == 1:
        raise ValueError(f"ERROR: Cannot phase variant {var_id}.")

    def get_allele(row):
        ids = row.GT.split("|")
        alleles = [row[{"0": "REF", "1": "ALT"}[aid]] for aid in ids]
        return alleles

    phased_snps["alleles"] = phased_snps.apply(get_allele, axis=1)
    phased_snps["amplicon_pos"] = phased_snps.POS - var_pos + amplicon_start
    phased_snps["amplicon_base"] = phased_snps.amplicon_pos.map(
        lambda i: amplicon_seq[i]
    )
    return phased_snps


def _get_allele_from_alignment(
    ref_aligned: str,
    query_aligned: str,
    offset: int,
    strand: int,
    start_pos: int,
    end_pos: int,
    positionwise_quality: np.ndarray = None,
    quality_thres: float = -1,
):
    # Include N, no quality filter
    assert len(ref_aligned) == len(query_aligned)
    allele = Allele()
    ref_gaps = 0
    alt_gaps = 0
    for i in range(len(ref_aligned)):
        if ref_aligned[i] == query_aligned[i]:
            continue
        ref_base = ref_aligned[i]
        alt_base = query_aligned[i]
        if ref_base == "-":
            ref_gaps += 1
        elif alt_base == "-":
            alt_gaps += 1
        ref_pos = i - ref_gaps
        allele.add(
            Edit(
                rel_pos=ref_pos,
                ref_base=ref_base,
                alt_base=alt_base,
                offset=offset,
                strand=strand,
            )
        )
    return allele


def get_edit(row):
    return _get_allele_from_alignment(
        row.Reference_Sequence,
        row.Aligned_Sequence,
        offset=0,
        strand=1,
        start_pos=phased_snps.amplicon_pos.min() - 1,
        end_pos=phased_snps.amplicon_pos.max() + 1,
    )


def filter_allele_by_pos(list_pos, allele, target_pos):
    edits = []
    for edit in allele.split(","):
        if edit == "":
            continue
        pos = int(edit.split(":")[0])
        if pos in list_pos or abs(pos - target_pos) < 5:
            edits.append(edit)
    return ",".join(sorted(edits))


def get_phase(allele, phased_snps, var_pos):
    phase = []
    allele_edits = {}
    if allele != "":
        edits = allele.split(",")
        for edit in edits:
            try:
                allele_edits[int(edit.split(":")[0])] = edit.split(":")[-1]
            except ValueError:
                print(edit)
    for pos in phased_snps.amplicon_pos.tolist():
        if pos == var_pos:
            continue
        if pos in allele_edits.keys():
            phase.append(allele_edits[pos].split(">")[-1])
        else:
            phase.append(
                phased_snps[["amplicon_pos", "amplicon_base"]]
                .set_index("amplicon_pos")
                .loc[pos, "amplicon_base"]
            )
    return phase


def filter_out_indels(allele):
    edits = []
    for edit in allele.edits:
        if edit.ref_base != "-" and edit.alt_base != "-":
            edits.append(str(edit))
    return ",".join(sorted(edits))


def format_allele_freq_table(path, phases, target_pos):
    aftbl = pd.read_csv(path, sep="\t")
    aftbl["allele"] = aftbl.apply(get_edit, axis=1)
    aftbl["be_allele"] = aftbl.allele.map(filter_out_indels)
    aftbl["phased_allele"] = aftbl.be_allele.map(
        lambda a: filter_allele_by_pos(phased_snps.amplicon_pos.tolist(), a, target_pos)
    )
    aftbl["clean"] = aftbl.be_allele.map(
        lambda a: len(a.split(","))
    ) == aftbl.phased_allele.map(lambda s: len(s.split(",")))
    aftbl = aftbl.loc[aftbl.clean, :]
    aftbl = (
        aftbl[["phased_allele", "#Reads", "n_mutated"]]
        .groupby("phased_allele")
        .sum()
        .reset_index()
    )
    aftbl = aftbl.groupby("phased_allele")[["#Reads"]].sum()
    aftbl["phase"] = aftbl.index.map(
        lambda a: get_phase(a, phased_snps, amplicon_start)
    )
    aftbl["assigned_phase"] = aftbl.phase.map(lambda ap: assign_phase(phases, ap))
    return aftbl


def mask_seq(seq, mask_pos):
    seq_array = np.array(list(seq))
    seq_array[mask_pos] = "N"
    return "".join(seq_array.tolist())


def assign_phase(phases, allele_phase):
    allele_phase_str = "".join(allele_phase)
    allele_phase = np.array(allele_phase)
    mask_pos = np.where(allele_phase == "N")[0]
    phases = phases.map(lambda p: np.array(list(p)))
    masked_phases = phases.map(lambda p: mask_seq(p, mask_pos))
    if sum(masked_phases == allele_phase_str) != 1:
        return -1
    else:
        return np.where(masked_phases == allele_phase_str)[0][0]


def get_phased_edit_rate(aftbl, editable_phase, amplicon_start):
    var_amplicon_base = phased_snps.loc[
        phased_snps.amplicon_pos == amplicon_start, "amplicon_base"
    ].item()
    aftbl.loc[~aftbl.assigned_phase.isin(editable_phase), "edited"] = np.nan
    if var_amplicon_base == var_alleles.ref.item():
        aftbl.loc[aftbl.assigned_phase.isin(editable_phase), "edited"] = aftbl.loc[
            aftbl.assigned_phase.isin(editable_phase), :
        ].index.map(
            lambda s: f"{amplicon_start}:{amplicon_start}:+:{var_alleles.ref}>{var_alleles.edited.item()}"
            in s
        )
    elif var_amplicon_base == var_alleles.edited.item():
        aftbl.loc[aftbl.assigned_phase.isin(editable_phase), "edited"] = aftbl.loc[
            aftbl.assigned_phase.isin(editable_phase), :
        ].index.map(
            lambda s: f"{amplicon_start}:{amplicon_start}:+:{var_amplicon_base}>{var_alleles.ref.item()}"
            not in s
        )
    return aftbl.groupby(["assigned_phase", "edited"])["#Reads"].sum()


def get_phased_reads(var_id, phases, target_pos):
    results = {}
    for rep in tqdm([1, 2, 3]):
        for cond in tqdm(["Ser", "SS"]):
            exp_id = f"{rep}-{cond}_20-locus"
            try:
                aftbl_atac = format_allele_freq_table(
                    f"results/atac_seq/CRISPRessoPooled_on_{rep}-{cond}_20-locus_ATAC_Seq_R1/CRISPResso_on_{var_id}_fw/Alleles_frequency_table.txt",
                    phases,
                    target_pos,
                )
            except FileNotFoundError:
                os.system(
                    f"unzip results/atac_seq/CRISPRessoPooled_on_{rep}-{cond}_20-locus_ATAC_Seq_R1/CRISPResso_on_{var_id}_fw/results_atac_seq_.Alleles_frequency_table.zip -d results/atac_seq/CRISPRessoPooled_on_{rep}-{cond}_20-locus_ATAC_Seq_R1/CRISPResso_on_{var_id}_fw/"
                )
                aftbl_atac = format_allele_freq_table(
                    f"results/atac_seq/CRISPRessoPooled_on_{rep}-{cond}_20-locus_ATAC_Seq_R1/CRISPResso_on_{var_id}_fw/Alleles_frequency_table.txt",
                    phases,
                    target_pos,
                )
            try:
                aftbl_gdna = format_allele_freq_table(
                    f"results/atac_seq/CRISPRessoPooled_on_{rep}-{cond}_20-locus_GDNA_R1/CRISPResso_on_{var_id}_fw/Alleles_frequency_table.txt",
                    phases,
                    target_pos,
                )
            except FileNotFoundError:
                os.system(
                    f"unzip results/atac_seq/CRISPRessoPooled_on_{rep}-{cond}_20-locus_GDNA_R1/CRISPResso_on_{var_id}_fw/results_atac_seq_.Alleles_frequency_table.zip -d results/atac_seq/CRISPRessoPooled_on_{rep}-{cond}_20-locus_GDNA_R1/CRISPResso_on_{var_id}_fw/"
                )
                aftbl_gdna = format_allele_freq_table(
                    f"results/atac_seq/CRISPRessoPooled_on_{rep}-{cond}_20-locus_GDNA_R1/CRISPResso_on_{var_id}_fw/Alleles_frequency_table.txt",
                    phases,
                    target_pos,
                )
            results[exp_id] = [
                get_phased_edit_rate(aftbl_atac, phases, amplicon_start),
                get_phased_edit_rate(aftbl_gdna, phases, amplicon_start),
            ]
    res = []
    for k, dfs in results.items():
        print(dfs)
        catdf = pd.concat(dfs, axis=1)
        catdf.columns = ["atac", "gdna"]
        catdf = catdf.T
        catdf.columns = ["unedited", "edited"]
        catdf["exp"] = k
        res.append(catdf)
    results_df = pd.concat(res)
    results_df.columns = ["unedited", "edited", "exp"]
    return results_df


if __name__ == "__main__":
    args = get_parser()

    # load amplicon information
    amplicon_info = pd.read_csv(
        "resources/atac_seq/030123_ATACseq_info.csv",
    ).set_index("gRNA")
    amplicon_start = amplicon_info.loc[args.variant_id, "Varpos_fw"] - 4
    amplicon_seq = amplicon_info.loc[
        args.variant_id, "Amplicon_fw (reference)"
    ].replace("N", "")
    amplicon_len = len(amplicon_seq)
    amplicon_end = amplicon_start + amplicon_len

    # variant into table
    var_info_tbl = pd.read_excel(
        "resources/LDLvar/20221013_LDLvar_simpleZscores_credset.xlsx",
        index_col=0,
        header=1,
    )
    var_alleles = var_info_tbl.loc[args.variant_id, ["A1", "A2"]].reset_index()
    var_alleles["edited"] = var_alleles[args.variant_id].map(
        lambda b: edit_dict[b] if b in edit_dict.keys() else np.nan
    )
    var_path = (
        "/data/pinello/PROJECTS/2021_08_ANBE/data/HepG2_reference/ENCFF713BPG.vcf.gz"
    )
    var_alleles = var_alleles.loc[~var_alleles.edited.isnull()]
    var_alleles.columns = ["allele", "ref", "edited"]
    var_chrom = int(var_info_tbl.loc[args.variant_id, "CHR"])
    var_pos = int(var_info_tbl.loc[args.variant_id, "position_hg19"])
    command = f"bcftools view -r chr{var_chrom}:{var_pos - amplicon_start+20}-{var_pos - amplicon_start + amplicon_len-20} {var_path} > results/atac_seq/{args.variant_id}_phasedSNPs.vcf"
    print(command)
    os.system(command)

    phased_snps = get_variant_info(args.variant_id, amplicon_start, amplicon_seq)
    phased_snps_nonvar = phased_snps.loc[phased_snps.amplicon_pos != amplicon_start]
    phases = pd.DataFrame(
        phased_snps_nonvar.alleles.tolist(),
        columns=[
            f"phase{i}" for i in range(len(phased_snps_nonvar.alleles.tolist()[0]))
        ],
    ).sum(axis=0)
    editable_phase = np.where(
        np.array(
            [
                p in edit_dict.keys()
                for p in phased_snps.loc[
                    phased_snps.amplicon_pos == amplicon_start, "alleles"
                ].item()
            ]
        )
    )[0]
    print(f"Phasing for {args.variant_id}:\n{phases}")
    result_df = get_phased_reads(args.variant_id, phases, amplicon_start)
    result_df.to_csv(f"results/atac_seq/{args.variant_id}_phased_counts.csv")

    y, X = dmatrices(
        "edited + unedited ~ exp + is_atac", data=result_df, return_type="dataframe"
    )
    exog = sm.add_constant(result_df.is_atac.values.astype(int), prepend=False)
    glm_binom = sm.GLM(y, X, family=sm.families.Binomial())
    glm_res = glm_binom.fit()
    glm_res.summary().tables[1].to_csv(
        f"results/atac_seq/{args.variant_id}_fit_result.csv"
    )
