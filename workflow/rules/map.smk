import itertools
ruleorder: map_plasmid > map_samples

LIBS = ["LDLvar", "LDLRCDS"]
REPS = {
    "LDLvar":list(range(1, 16)), #1-15
    "LDLRCDS":list(range(1, 10)) #1-9
}

def reps_to_bins(rep):
    if rep in [1, 2, 3, 4]:
        return ["top", "bot", "bulk"]
    else:
        return ["top", "high", "bulk", "low", "bot"]

def get_reps_bins(lib):
    reps = REPS[lib]
    return [f"rep{rep}_{sort_bin}" for sort_bin in reps_to_bins(rep) for rep in reps]

rule get_fastq_data:
    run:
        pass
    output:
        R1s = expand(
            ["results/raw/{wildcards.lib}/{wildcards.lib}_{rep_bin}_R1.fastq.gz"],
            rep_bin = get_reps_bins(wildcards.lib)
        )
        R2s = expand(
            ["results/raw/{wildcards.lib}/{wildcards.lib}_{rep_bin}_R2.fastq.gz"],
            rep_bin = get_reps_bins(wildcards.lib)
        )
        plasimd_R1 = "results/raw/{wildcards.lib}/{wildcards.lib}_{wildcards.rep_bin}_R1.fastq.gz",
        plasmid_R2 = "results/raw/{wildcards.lib}/{wildcards.lib}_{wildcards.rep_bin}_R2.fastq.gz",
        sample_list_file = 'results/raw/{wildcards.lib}/sample_list.csv',
        guide_start_seqs = "results/raw/{wildcards.lib}/guide_start_seqs.txt",
        guide_end_seqs = "results/raw/{wildcards.lib}/guide_end_seqs.txt"

rule write_sample_list_file:
    input:
        R1s = expand(
            ["results/raw/{wildcards.lib}/{wildcards.lib}_{rep_bin}_R1.fastq.gz"],
            rep_bin = get_reps_bins(wildcards.lib)
        )
        R2s = expand(
            ["results/raw/{wildcards.lib}/{wildcards.lib}_{rep_bin}_R2.fastq.gz"],
            rep_bin = get_reps_bins(wildcards.lib)
        )
    output:
        sample_list='results/raw/{wildcards.lib}/sample_list.csv'
    run:
        f = open(output.sample_list)
        for rep_bin in get_reps_bins(wildcards.lib):
            f.write(f"results/raw/{wildcards.lib}/{wildcards.lib}_{rep_bin}_R1.fastq.gz,results/raw/{wildcards.lib}/{wildcards.lib}_{rep_bin}_R2.fastq.gz,{rep_bin}\n")
        f.close()

rule map_plasmid:
    input:
        guide_info = 'resources/gRNA_info/{lib}_gRNA_beret.csv',
        plasimd_R1 = "results/raw/{lib}/{lib}_{rep_bin}_R1.fastq.gz",
        plasmid_R2 = "results/raw/{lib}/{lib}_{rep_bin}_R2.fastq.gz",
    params:
        output_dir = 'results/mapped/{lib}'
    output:
        out_h5ad = 'results/mapped/{lib}/beret_count_plasmid_{lib}.h5ad'
    run:
        shell('mkdir -p {params.output_dir}')
        shell(
            "beret-count --R1 results/raw/{wildcards.lib}/{wildcards.lib}_plasmid_R1.fastq --R2 results/raw/{wildcards.lib}_plasmid_R2.fastq -b A -f {input.guide_info} -o {params.output_dir} -a")

rule map_samples:
    input:
        guide_info = 'resources/gRNA_info/{lib}_gRNA_beret.csv',
        sample_list = 'results/raw/{lib}/sample_list.csv',# ?
        guide_start_seqs = "results/raw/{lib}/guide_start_seqs.txt",
        guide_end_seqs = "results/raw/{lib}/guide_end_seqs.txt",

    params:
        output_dir = 'results/mapped/{lib}/'
    output:
        out_h5ad = 'results/mapped/{lib}/beret_count_{lib}.h5ad'
    run:
        shell('mkdir -p {params.output_dir}')
        shell(
            "cd /data/pinello/PROJECTS/2021_08_ANBE/data/{wildcards.seq_data}/crisprep_counts/{wildcards.lib}/; "+
            "beret-count-samples --input sample_list.csv -b A -f {input.guide_info} -o {params.output_dir} -a -t 12 --name {wildcards.seq_data}_{wildcards.lib} --guide_start_seqs_file={input.guide_start_seqs} --guide_end_seqs_file={input.guide_end_seqs}{[wildcards.]}GGAAAGGACGAAACACCG")

rule map_all:
    input:
        expand(['results/mapped/{lib}/beret_count_{lib}.h5ad'], lib=LIBS),
        expand(['results/mapped/{lib}/beret_count_plasmid_{lib}.h5ad'], lib=LIBS)
