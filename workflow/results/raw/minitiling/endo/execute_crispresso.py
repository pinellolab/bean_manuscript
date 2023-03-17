import os
import shlex
import subprocess
import pandas as pd

output_dir = "q28_n20_guides_full"
os.makedirs(output_dir, exist_ok = True)
amplicon_seq = {"A":"CTTGTCCAATGTGAGGCGGTggaggcggaggcgggcgtcgggaggacggggcttgtgtacgagcggggcggggctggcgcggaagtctgagcctcaccttgtccggggcgaggcggatgcaggggaggcctggcgttcctccgcggttcctgTCACAAAGGCGACGACAAG",
        "B":"CGAAGGACTGGAGTGGGAATCagagcttcacgggttaaaaagccgatgtcacatcggccgttcgaaactcctcctcttgcagtgaggtgaagacatttgaaaatcaccccactgcaaactcctccccctgctagaaacctcacattgAAATGCTGTAAATGACGTGGGC",
        "C":"CTGGAGCAAGCCTTACCTGCagtccccgccgcggcgaggagcaaggcgacggtccagcgcaatttccagccccagggccccatgctcgcagcctctgccaggcagtgtcccgacccggatcacgacctgctgtgtcctagctggaaaccctggcttcccgcgattgcactcggGGCCCACGTCATTTACAGCATT",
        "D":"CACCTCAAAGACGGCCAAGGagaaggggtgggccagcctcttttcatcctccaagatggtcttccggttgcccccgttgacatcgatgcttgagatggagtgaagtttggagtcaacccagtagaggcggccactgaggagatctagacACACAAGTGGATAAGGAGAGGTCA"}

tbl = pd.read_excel("../0322_LDLRendo_ANBE_libraries.xlsx")
def arg_join(str_list):
    str_list = [x.replace(" ", "-") for x in str_list]
    return ','.join(str_list)

sublib_info = tbl.groupby("Sublibrary").agg({"name":arg_join,"sequence":arg_join}) #.to_csv("0322_LDLRendo_ANBE_libraries_crispresso_args.tsv",sep="\t")

sublibrary_experiments = {"A":["Q5"],
        "B":["Q5"],
        "C":["Q5", "Q5U"],
        "D":["Q5", "Q5U"]}

replicates = ["rep1", "rep2", "rep3", "WT_Ctrl"]
fastq_base_name="HepG2_LDLREndo_{}_{}_Endogenous_gDNA_R1"
#crispresso_command = "CRISPResso --fastq_r1 {fastq} --amplicon_seq {amplicon} -gn {guide_names} -g {guide_sequence} --base_editor_output --conversion_nuc_from A --conversion_nuc_to G --needleman_wunsch_gap_open -20 --needleman_wunsch_gap_extend -10 --needleman_wunsch_gap_incentive -10 -wc 10 -w 10 --exclude_bp_from_left 3 --exclude_bp_from_right 0 --needleman_wunsch_aln_matrix_loc /data/pinello/PROJECTS/2021_08_ANBE/data/endo_site/aln_mat_ABE.txt -q 28 --min_bp_quality_or_N 30 -o {outdir} & "
crispresso_command = "CRISPResso --fastq_r1 {fastq} --amplicon_seq {amplicon} -gn {guide_names} -g {guide_sequence} --base_editor_output --conversion_nuc_from A --conversion_nuc_to G --needleman_wunsch_gap_open -20 --needleman_wunsch_gap_extend -10 --needleman_wunsch_gap_incentive -10 --plot_window_size=10 --quantification_window_center -10 --exclude_bp_from_left 0 --exclude_bp_from_right 0 --needleman_wunsch_aln_matrix_loc /data/pinello/PROJECTS/2021_08_ANBE/data/endo_site/aln_mat_ABE.txt -q 28 --min_bp_quality_or_N 20 -o {outdir} -v 1"
crispresso_command = "CRISPResso --fastq_r1 {fastq} --amplicon_seq {amplicon} -gn {guide_names} -g {guide_sequence} --base_editor_output --conversion_nuc_from A --conversion_nuc_to G --needleman_wunsch_gap_open -20 --needleman_wunsch_gap_extend -10 --needleman_wunsch_gap_incentive -10 --plot_window_size=10 --quantification_window_center -10 --exclude_bp_from_left 0 --exclude_bp_from_right 0 --needleman_wunsch_aln_matrix_loc /data/pinello/PROJECTS/2021_08_ANBE/data/endo_site/aln_mat_ABE.txt -q 28 --min_bp_quality_or_N 20 -o {outdir} -v 1"

jobs_count = 0
processes = []
for sublibrary in ["A", "B", "C", "D"]:
    for exp in sublibrary_experiments[sublibrary]:
        sample = "Sub{}_{}".format(sublibrary, exp)
        for rep in replicates:
            fastq_base = fastq_base_name.format(sample, rep)
            fastq_str = "fastq/" + fastq_base + ".fastq.gz"

            sample_command = crispresso_command.format(
                    fastq = fastq_str, 
                    amplicon = amplicon_seq[sublibrary], 
                    guide_names = sublib_info.loc[sublibrary, "name"], 
                    guide_sequence = sublib_info.loc[sublibrary, "sequence"],
                    outdir=output_dir,
                    base_name=fastq_base)
            sample_command_list = shlex.split(sample_command.replace("'", "\\'"))
            print(sample)
            jobs_count += 1
            p = subprocess.Popen(sample_command_list)
            processes.append(p)
    
while processes: 
    for p in processes: 
        if p.poll() is not None:
           processes.remove(p) 
           print('{} done, status {}'.format(p.args[2], p.returncode))
           break
            
print("submitted {} jobs.".format(jobs_count))
