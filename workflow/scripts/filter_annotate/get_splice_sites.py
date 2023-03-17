import argparse
import re
import numpy as np
import pandas as pd


def parse_args():
    argparser = argparse.ArgumentParser(
        "Get splice site position",
        usage="Get splice site position from exon fasta and target editing base.",
    )

    argparser.add_argument(
        "exon_fa_path", help="File path to fasta file with exon position information."
    )
    argparser.add_argument("edited_base", help="Edited base, either A or C.")
    argparser.add_argument(
        "output_path", help="output path of the splice site csv file."
    )

    return argparser.parse_args()


def get_splice_positions(exon_fa_path):
    splice_donor_pos = []
    splice_acceptor_pos = []
    p = re.compile("range=chr(\d+):(\d+)-(\d+)")
    with open(exon_fa_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                result = p.search(line)
                exon_start = int(result[2])
                exon_end = int(result[3])
                splice_donor_pos.append(exon_end)
                splice_acceptor_pos.append(exon_start)
    splice_donor_pos = np.array(splice_donor_pos)
    splice_acceptor_pos = np.array(splice_acceptor_pos)
    splice_donor_pos = splice_donor_pos[:-1]
    splice_acceptor_pos = splice_acceptor_pos[1:]
    return splice_donor_pos, splice_acceptor_pos


def get_targetable_splice_positions(
    splice_donor_pos, splice_acceptor_pos, edited_base="A"
):
    """
    Splice donor: GT
    Splice acceptor: AG
    """
    splice_site_dfs = []
    revcomp_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
    splice_donor_consensus = {"G": 1, "T": 2}
    splice_acceptor_consensus = {"A": -2, "G": -1}
    for splice_site_pos, rel_basepos_map, label in zip(
        [splice_donor_pos, splice_acceptor_pos],
        [splice_donor_consensus, splice_acceptor_consensus],
        ["SD", "SA"],
    ):
        splice_site_dfs.extend(
            pd.DataFrame(
                {
                    "pos": splice_site_pos + rel_basepos_map[base],
                    "type": label,
                    "target_base": base,
                }
            )
            for base in [edited_base, revcomp_map[edited_base]]
            if base in rel_basepos_map
        )
    return pd.concat(splice_site_dfs)


if __name__ == "__main__":
    args = parse_args()
    sd_pos, sa_pos = get_splice_positions(args.exon_fa_path)
    splice_target_df = get_targetable_splice_positions(sd_pos, sa_pos, args.edited_base)
    splice_target_df.to_csv(args.output_path)
