rule run_crispresso:
    run:
        shell("python scripts/atac_seq/execute_crispresso_indiv_demuxed.py")

rule analyze_output:
    output:
        
    run:
        shell("python scripts/atac_seq/analyze_outputs.py")