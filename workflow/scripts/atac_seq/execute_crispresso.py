import os
import shlex
import subprocess
import pandas as pd


reps = [1, 2, 3, "C"]
subs = ["Ser", "SS"]
libraries = ["ATAC_Seq", "GDNA"]


processes = []
for rep in reps:
    for sub in subs:
        for lib in libraries:
            file_path = f"results/raw/atac_seq/{rep}-{sub}_20-locus_{lib}_R1.fastq.gz"
            output_path = f"results/atac_seq/crispresso_runs/CRISPRessoPooled_on_{rep}-{sub}_20-locus_{lib}_R1.html"
            if not os.path.exists(output_path):
                crispresso_command = f"CRISPRessoPooled -r1 ../../../{file_path} -f ../../../resources/atac_seq/030123_amplicon_info_CRISPRessoInput_fw.tsv --base_editor_output --conversion_nuc_from A --conversion_nuc_to G --exclude_bp_from_left 20 --exclude_bp_from_right 20 --needleman_wunsch_aln_matrix_loc /data/pinello/PROJECTS/2021_08_ANBE/data/endo_site/aln_mat_ABE.txt -q 28 --min_bp_quality_or_N 20 --file_prefix=results/atac_seq/ --min_reads_to_use_region=1 --skip_failed"
                sample_command_list = shlex.split(
                    crispresso_command.replace("'", "\\'")
                )
                p = subprocess.Popen(
                    sample_command_list,
                    cwd="results/atac_seq/202305054_crispresso_runs_allow_indels/",
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
