import sys

lib = sys.argv[1]
path = sys.argv[2]

ALL_REPS = {  # including technical replicates
    "LDLvar": list(range(1, 16)) + [1.1, 2.1, 3.1, 4.1, 5.1],  # 1-15
    "LDLRCDS": list(range(1, 10)) + [1.1, 2.1, 3.1, 4.1, 5.1, 6.1],  # 1-9
    "LDLRCDS_CBE_SpRY": list(range(1, 5)),
    "LDLRCDS_CBE_CasNG": list(range(1, 5)),
}


def reps_to_bins(lib, rep):
    if lib in ["LDLvar", "LDLRCDS"] and rep in [1, 2, 3, 4, 1.1, 2.1, 3.1, 4.1]:
        return ["top", "bot", "bulk"]
    else:
        return ["top", "high", "bulk", "low", "bot"]


def get_reps_bins(lib):
    reps = ALL_REPS[lib]
    res = []
    for rep in reps:
        res += [f"rep{rep}_{sort_bin}" for sort_bin in reps_to_bins(lib, rep)]
    return res


f = open(path, "w")
for rep_bin in get_reps_bins(lib):
    f.write(
        f"results/raw/{lib}/{lib}_{rep_bin}_R1.fastq.gz,results/raw/{lib}/{lib}_{rep_bin}_R2.fastq.gz,{rep_bin}\n"
    )
f.close()
