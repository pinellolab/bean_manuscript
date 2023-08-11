import os
import shlex
import subprocess
import pandas as pd


reps = [1, 2, 3, "C"]
subs = ["Ser", "SS"]
libraries = ["ATAC_Seq", "GDNA"]

amplicons_tbl = pd.read_csv(
    "resources/atac_seq/030123_amplicon_info_CRISPRessoInput_fw.tsv",
    sep="\t",
    header=None,
)
samples = ",".join(amplicons_tbl.values[:, 0].tolist())
amplicon_seqs = ",".join(amplicons_tbl.values[:, 1].tolist())
processes = []
for rep in reps:
    for sub in subs:
        for lib in libraries:
            file_path = f"results/raw/atac_seq/{rep}-{sub}_20-locus_{lib}_R1.fastq.gz"
            # for i, sample in enumerate(samples):
            # amplicon_seq = amplicon_seqs[i]
            output_path = f"results/atac_seq/crispresso_runs_indiv//CRISPRessoPooled_on_{rep}-{sub}_20-locus_{lib}_R1.html"
            if not os.path.exists(output_path):
                crispresso_command = f"CRISPResso --fastq_r1 ../../raw/atac_seq/{rep}-{sub}_20-locus_{lib}_R1.fastq.gz --amplicon_seq {amplicon_seqs} --amplicon_name {samples} --base_editor_output --conversion_nuc_from A --conversion_nuc_to G --exclude_bp_from_left 20 --exclude_bp_from_right 20 --needleman_wunsch_aln_matrix_loc /data/pinello/PROJECTS/2021_08_ANBE/data/endo_site/aln_mat_ABE.txt -q 30 --min_bp_quality_or_N 20 "

                sample_command_list = shlex.split(
                    crispresso_command.replace("'", "\\'")
                )
                p = subprocess.Popen(
                    sample_command_list,
                    cwd="results/atac_seq/crispresso_runs_indiv/",
                )
                processes.append(p)
            else:
                print(f"Output file {output_path} exists.")
print(f"submitted {len(processes)} jobs.")
while processes:
    for p in processes:
        if p.poll() is not None:
            processes.remove(p)
            print("{} done, status {}".format(p.args[2], p.returncode))
            break
