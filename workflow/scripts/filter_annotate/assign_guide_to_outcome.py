import argparse
import os
import re
from typing import List

import bean as be
import numpy as np
import pandas as pd
from bean.annotate.translate_allele import CDS, RefBaseMismatchException


def edit_range(seq, start=6 + 3, end=6 + 7, ref_base="A", alt_base="G"):
    return seq[:start] + seq[start:end].replace(ref_base, alt_base) + seq[end:]


def get_allele(seq, offset, edit_start=6 + 3, edit_end=6 + 7, strand=1):
    allele = be.Allele()
    for i in range(edit_start, edit_end):
        if seq[i] == "A":
            allele.add(be.Edit(i, "A", "G", offset=offset, strand=strand))
    return allele


def _get_allele(row):
    if row.strand == "neg":
        strand = -1
        offset = row.start_pos + 32 - 6 - 1
    else:
        strand = 1
        offset = row.start_pos - (32 - 6 - row.guide_len)
    return get_allele(
        row.Reporter,  # TODO: will the column name be changed?
        offset,
        edit_start=32 - 6 - row.guide_len + 3,
        edit_end=32 - 6 - row.guide_len + 7,
        strand=strand,
    )


def edit_alleles(allele_str, include_syn=True):
    ldlr_cds = CDS()
    ldlr_cds.edit_allele(allele_str)
    return ldlr_cds.get_aa_change(include_syn)


def translate_lenient(s, include_syn=True):
    try:
        return edit_alleles(s, include_syn)
    except RefBaseMismatchException as e:
        print(e)
        return "ref_mismatch"


def _format_behive_edit(edit_string: str):
    formatted_edits = []
    severity = []
    for edit_str in edit_string.split(","):
        edit_str = edit_str.strip()
        if len(edit_str.split("-")) == 4:  # Noncoding edit
            chrom, pos, ref, alt = edit_str.split("-")
            assert chrom == "19"  # TODO: should be relaxed later on
            formatted_edits.append(f"{pos}:{ref}>{alt}")
            severity.append(int(ref != alt))
        elif re.fullmatch("[*A-Z]\d+[A-Z]", edit_str, flags=0):  # amino acid edit
            ref = edit_str[0]
            pos = edit_str[1:-1]
            alt = edit_str[-1]
            formatted_edits.append(f"A{pos}:{ref}>{alt}")
            severity.append(0.5 * float(ref != alt))
        else:
            print(
                f"Edit string {edit_str} didn't match known pattern in BE-Hive input."
            )
            formatted_edits.append(edit_str)
    if len(formatted_edits) > 1:
        return formatted_edits[np.argmax(severity)]
    else:
        return formatted_edits[0]


def format_behive(behive_pred_path: str) -> pd.DataFrame:
    behive_pred = pd.read_excel(behive_pred_path)
    behive_pred["edit_str"] = behive_pred.apply(
        lambda row: row.aa_str if isinstance(row.aa_str, str) else row.variant_str,
        axis=1,
    )
    behive_pred.edit_str = behive_pred.edit_str.map(_format_behive_edit)
    return behive_pred


def classify_behive_outcome(
    behive_pred: pd.DataFrame, guide_info: pd.DataFrame
) -> pd.DataFrame:
    behive_pred["name"] = behive_pred["id"].map(lambda s: s.replace("-", "' "))
    behive_pred["n_edits"] = behive_pred.edit_str.map(lambda s: len(s.split(",")))
    guide_info = pd.merge(
        guide_info,
        behive_pred[["name", "edit_str", "consequence_terms"]],
        on="name",
        how="left",
    )
    guide_info = guide_info.rename(
        columns={"edit_str": "target_behive", "consequence_terms": "group_behive"}
    )
    return guide_info


def get_behive_outcome(
    behive_pred_path: str,
    guide_info: pd.DataFrame,
    ctrl_behive_pred_path: str = "resources/gRNA_info/target_prediction/control_gRNA_BEHive_consequence.csv",
) -> pd.DataFrame:
    """Read predicted BE-Hive annotations from precomputed files."""
    behive_pred = format_behive(behive_pred_path)
    guide_info = classify_behive_outcome(behive_pred, guide_info)
    ctrl_behive_pred = pd.read_csv(ctrl_behive_pred_path, index_col=0)
    guide_info.loc[
        guide_info.Region.isin(["CBE control", "ABE control"]), "target_behive"
    ] = ctrl_behive_pred.loc[
        guide_info[guide_info.Region.isin(["CBE control", "ABE control"])].name,
        "annotation_str",
    ]
    guide_info.loc[guide_info.Region == "CBE control", "group_behive"] = "CBE control"
    guide_info.loc[guide_info.Region == "ABE control", "group_behive"] = "ABE control"
    return guide_info


def get_nt_severity(nt_allele, region, is_splice):
    if is_splice:
        return 1.5
    if len(nt_allele.edits) == 0:
        return -1
    if "UTR" in region:
        return 0.2
    return 0.5


group_from_severity = {
    1.5: "splicing",
    -1: "no_edit",
    0.2: "UTR",
    0.5: "intron",
    0: "synonymous",
    1: "missense",
    2: "nonsense",
}


def get_most_severe_edit(row, splice_positions: List):
    if row.severity == 2:
        for aa_edit in row.coding.edits:
            if aa_edit.alt == "*":
                return aa_edit
    if row.severity == 1.5:
        for nt_edit in row.noncoding.edits:
            if nt_edit.pos in splice_positions:
                return nt_edit
    if row.severity == 1:
        for aa_edit in row.coding.edits:
            if aa_edit.alt != aa_edit.ref:
                return aa_edit
    if row.severity == 0.5 or row.severity == 0.2:
        return next(iter(row.noncoding.edits))
    if row.severity == 0:
        return next(iter(row.coding.edits))
    if row.severity == -1:
        return ""


def get_allEdited_outcome(
    bdata: be.ReporterScreen, splice_sites: pd.DataFrame
) -> pd.DataFrame:
    # http://localhost:20081/notebooks/PROJECTS/2021_08_ANBE/mageck_results/LDLRCDS_fullsort/allEdited/all_pos_edited.ipynb
    bdata.guides.loc[
        bdata.guides.Region.isin(["ABE control", "CBE control"]), "start_pos"
    ] = 0
    bdata.guides["allele"] = bdata.guides.apply(_get_allele, axis=1)
    aa_allele = pd.DataFrame({"aa_allele": bdata.guides.allele.map(translate_lenient)})
    aa_allele["coding"] = aa_allele.aa_allele.map(
        lambda a: be.AminoAcidAllele() if isinstance(a, str) else a.aa_allele
    )
    aa_allele["noncoding"] = aa_allele.aa_allele.map(
        lambda a: be.Allele() if isinstance(a, str) else a.nt_allele
    )
    aa_allele["is_splice"] = aa_allele.noncoding.map(
        lambda a: any([e.pos in splice_sites.pos.tolist() for e in a.edits])
    )
    aa_allele["aa_severity"] = aa_allele.coding.map(lambda a: a.get_most_severe())
    aa_allele["Region"] = bdata.guides["Region"]
    aa_allele["noncoding_severity"] = aa_allele.apply(
        lambda row: get_nt_severity(row.noncoding, row.Region, row.is_splice), axis=1
    )
    aa_allele["severity"] = aa_allele.apply(
        lambda row: max(row.aa_severity, row.noncoding_severity), axis=1
    )
    aa_allele["most_severe_edit"] = aa_allele.apply(
        lambda row: get_most_severe_edit(row, splice_sites.pos.tolist()), axis=1
    )
    aa_allele["most_severe_edit_abs"] = aa_allele.most_severe_edit.map(
        lambda a: a.get_abs_edit() if isinstance(a, be.Edit) else a
    )

    bdata.guides["most_severe_edit"] = aa_allele.most_severe_edit_abs

    # Assign unique ID for edits introduced by non-targeting controls
    ctrl_with_edit = bdata.guides.Region.isin(["ABE control", "CBE control"]) & (
        bdata.guides.most_severe_edit != ""
    )
    bdata.guides.loc[ctrl_with_edit, "most_severe_edit"] = (
        bdata.guides.loc[ctrl_with_edit, :].index
        + "_"
        + bdata.guides.loc[ctrl_with_edit, "most_severe_edit"]
    )
    bdata.guides["most_severe_edit"] = bdata.guides.most_severe_edit.map(str)
    bdata.guides["severity"] = aa_allele.severity
    bdata.guides["group"] = bdata.guides.severity.map(lambda s: group_from_severity[s])
    bdata.guides = bdata.guides.drop(columns="allele")
    bdata.guides.loc[bdata.guides.Region == "ABE control", "group"] = "ABE control"
    bdata.guides.loc[bdata.guides.Region == "CBE control", "group"] = "CBE control"
    bdata.guides = bdata.guides.rename(
        columns={"most_severe_edit": "target_allEdited", "group": "group_allEdited"}
    )

    return bdata.guides


def main():
    def validate_args(args):
        if args.mode == "all" or args.mode == "both":
            assert os.path.exists(args.splice_site_csv)
        elif args.mode == "pred" or args.mode == "both":
            assert os.path.exists(args.pred_path)

    parser = argparse.ArgumentParser(description="Assign guide to editing outcome.")
    parser.add_argument(
        "mode",
        type=str,
        help="Mode of editing. Needs to be either [all, pred]. `all` assumes all nucleotides in editing window are edited. `behive` needs external prediction (ex. from BE-Hive) per guide, where the prediction outcome should be fed in as argument.",
    )
    parser.add_argument("bdata_path", type=str, help="Path of an ReporterScreen object")
    parser.add_argument(
        "output_path", type=str, help="Output path where csv file will be written."
    )
    parser.add_argument(
        "--write-bdata",
        default=False,
        action="store_true",
        help="Output is bdata with outcome assignment. If false, guide information dataframe is written.",
    )
    parser.add_argument(
        "--splice-site-csv",
        "-s",
        type=str,
        help="Path of text file that has targetable position at splicing site. Needs to have `pos` as position.",
    )
    parser.add_argument(
        "--pred_path", "-p", type=str, help="Path of predicted editing outcome"
    )  # TODO: needs to be simplified to plain table instead of excel
    args = parser.parse_args()
    validate_args(args)

    bdata = be.read_h5ad(args.bdata_path)
    if args.mode == "all" or args.mode == "both":
        splice_sites = pd.read_csv(args.splice_site_csv)
        guide_info = get_allEdited_outcome(bdata, splice_sites)
        if args.write_bdata:
            bdata.guides = guide_info
    if args.mode == "pred" or args.mode == "both":
        guide_info = get_behive_outcome(args.pred_path, bdata.guides)
        if args.write_bdata:
            if "name" in guide_info.columns:
                guide_info = guide_info.set_index("name")
            bdata.guides = guide_info
    if args.write_bdata:
        try:
            bdata.write(args.output_path)
        except TypeError as e:
            print(e)
            print(bdata.obs.columns)
            bdata.guides.to_csv("tmp")
    else:
        guide_info.to_csv(args.output_path)


if __name__ == "__main__":
    main()
