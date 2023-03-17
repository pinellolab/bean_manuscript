A_seq=CTTGTCCAATGTGAGGCGGTggaggcggaggcgggcgtcgggaggacggggcttgtgtacgagcggggcggggctggcgcggaagtctgagcctcaccttgtccggggcgaggcggatgcaggggaggcctggcgttcctccgcggttcctgTCACAAAGGCGACGACAAG
B_seq=CGAAGGACTGGAGTGGGAATCagagcttcacgggttaaaaagccgatgtcacatcggccgttcgaaactcctcctcttgcagtgaggtgaagacatttgaaaatcaccccactgcaaactcctccccctgctagaaacctcacattgAAATGCTGTAAATGACGTGGGC
C_seq=CTGGAGCAAGCCTTACCTGCagtccccgccgcggcgaggagcaaggcgacggtccagcgcaatttccagccccagggccccatgctcgcagcctctgccaggcagtgtcccgacccggatcacgacctgctgtgtcctagctggaaaccctggcttcccgcgattgcactcggGGCCCACGTCATTTACAGCATT
D_seq=CACCTCAAAGACGGCCAAGGagaaggggtgggccagcctcttttcatcctccaagatggtcttccggttgcccccgttgacatcgatgcttgagatggagtgaagtttggagtcaacccagtagaggcggccactgaggagatctagacACACAAGTGGATAAGGAGAGGTCA

pids=()
for cond in SubA_Q5; do for rep in rep1 rep2 rep3 WT_Ctrl; do CRISPResso --fastq_r1 results/raw/endo/fastq/HepG2_LDLREndo_${cond}_${rep}_Endogenous_gDNA_R1.fastq.gz --amplicon_seq $A_seq --base_editor_output --conversion_nuc_from A --conversion_nuc_to G --needleman_wunsch_gap_open -20 --needleman_wunsch_gap_extend -10 --needleman_wunsch_gap_incentive -10 -qwc 4-150 --exclude_bp_from_left 3 --exclude_bp_from_right 0 --needleman_wunsch_aln_matrix_loc /data/pinello/PROJECTS/2021_08_ANBE/data/endo_site/aln_mat_ABE.txt -q 28 --min_bp_quality_or_N 20 -o results/mapped/minitiling/endo/ & 
pids+=($!)
done; done

for cond in SubB_Q5; do for rep in rep1 rep2 rep3 WT_Ctrl; do CRISPResso --fastq_r1 results/raw/endo/fastq/HepG2_LDLREndo_${cond}_${rep}_Endogenous_gDNA_R1.fastq.gz --amplicon_seq $B_seq --base_editor_output --conversion_nuc_from A --conversion_nuc_to G --needleman_wunsch_gap_open -20 --needleman_wunsch_gap_extend -10 --needleman_wunsch_gap_incentive -10 -qwc 4-150 --exclude_bp_from_left 3 --exclude_bp_from_right 0 --needleman_wunsch_aln_matrix_loc /data/pinello/PROJECTS/2021_08_ANBE/data/endo_site/aln_mat_ABE.txt -q 28 --min_bp_quality_or_N 20 -o results/mapped/minitiling/endo/ & 
pids+=($!)
done; done

for cond in SubC_Q5 SubC_Q5U ;do for rep in rep1 rep2 rep3 WT_Ctrl; do CRISPResso --fastq_r1 results/raw/endo/fastq/HepG2_LDLREndo_${cond}_${rep}_Endogenous_gDNA_R1.fastq.gz --amplicon_seq $C_seq --base_editor_output --conversion_nuc_from A --conversion_nuc_to G --needleman_wunsch_gap_open -20 --needleman_wunsch_gap_extend -10 --needleman_wunsch_gap_incentive -10 -qwc 4-150 --exclude_bp_from_left 3 --exclude_bp_from_right 0 --needleman_wunsch_aln_matrix_loc /data/pinello/PROJECTS/2021_08_ANBE/data/endo_site/aln_mat_ABE.txt -q 28 --min_bp_quality_or_N 20 -o results/mapped/minitiling/endo/ & 
pids+=($!)
done; done

for cond in SubD_Q5 SubD_Q5U ;do for rep in rep1 rep2 rep3 WT_Ctrl; do CRISPResso --fastq_r1 results/raw/endo/fastq/HepG2_LDLREndo_${cond}_${rep}_Endogenous_gDNA_R1.fastq.gz --amplicon_seq $D_seq --base_editor_output --conversion_nuc_from A --conversion_nuc_to G --needleman_wunsch_gap_open -20 --needleman_wunsch_gap_extend -10 --needleman_wunsch_gap_incentive -10 -qwc 4-150 --exclude_bp_from_left 3 --exclude_bp_from_right 0 --needleman_wunsch_aln_matrix_loc /data/pinello/PROJECTS/2021_08_ANBE/data/endo_site/aln_mat_ABE.txt -q 28 --min_bp_quality_or_N 20 -o results/mapped/minitiling/endo/ & 
pids+=($!)
done; done

for pid in ${pids[*]}; do
    wait $pid
done