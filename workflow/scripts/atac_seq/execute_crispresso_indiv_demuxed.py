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
var_ids = [
    amplicon_name.split("_fw")[0]
    for amplicon_name in amplicons_tbl.values[:, 0].tolist()
]
amplicon_seqs = amplicons_tbl.values[:, 1].tolist()
processes = []
for i, var_id in enumerate(var_ids):
    for rep in reps:
        for sub in subs:
            for lib in libraries:
                file_path = f"../../../results/raw/atac_seq/demuxed/{rep}-{sub}_20-locus_{lib}.{var_id}.fq"
                if var_id == "rs4719853_Maj_ABE_214_trimmed2":
                    # file_path = f"../../../results/raw/atac_seq/demuxed/{rep}-{sub}_20-locus_{lib}.rs4719853_Maj_ABE_214.fq"
                    os.system(
                        f"rm results/raw/atac_seq/demuxed/{rep}-{sub}_20-locus_{lib}.{var_id}.fq; ln -s {rep}-{sub}_20-locus_{lib}.rs4719853_Maj_ABE_214.fq results/raw/atac_seq/demuxed/{rep}-{sub}_20-locus_{lib}.{var_id}.fq"
                    )
                amplicon_seq = amplicon_seqs[i]
                output_path = f"results/atac_seq/crispresso_runs_indiv_demuxed/CRISPResso_on_{rep}-{sub}_20-locus_{lib}.{var_id}.html"
                if not os.path.exists(output_path):
                    crispresso_command = f"CRISPResso --fastq_r1 {file_path} --amplicon_seq {amplicon_seq} --amplicon_name {var_id} --base_editor_output --conversion_nuc_from A --conversion_nuc_to G --needleman_wunsch_aln_matrix_loc /data/pinello/PROJECTS/2021_08_ANBE/data/endo_site/aln_mat_ABE.txt -q 30 --min_bp_quality_or_N 20 "
                    if var_id != "rs4719853_Maj_ABE_214_trimmed":
                        crispresso_command += (
                            "--exclude_bp_from_left 20 --exclude_bp_from_right 20 "
                        )
                    else:
                        crispresso_command += (
                            "--exclude_bp_from_left 0 --exclude_bp_from_right 0 "
                        )

                    sample_command_list = shlex.split(
                        crispresso_command.replace("'", "\\'")
                    )
                    p = subprocess.Popen(
                        sample_command_list,
                        cwd="results/atac_seq/crispresso_runs_indiv_demuxed/",
                    )
                    processes.append(p)
                    # processes.append(sample_command_list)
                else:
                    pass  # print(f"Output file {output_path} exists.")


print(f"submitted {len(processes)} jobs.")

while processes:
    for p in processes:
        if p.poll() is not None:
            processes.remove(p)
            print("{} done, status {}".format(p.args[2], p.returncode))
            break
