import sys


ALL_REPS = {  # including technical replicates
    "LDLvar": list(range(1, 16)) + [1.1, 2.1, 3.1, 4.1, 5.1],  # 1-15
    "LDLRCDS": list(range(1, 10)) + [1.1, 2.1, 3.1, 4.1, 5.1, 6.1],  # 1-9
}


def reps_to_bins(rep):
    if rep in [1, 2, 3, 4, 1.1, 2.1, 3.1, 4.1]:
        return ["top", "bot", "bulk"]
    else:
        return ["top", "high", "bulk", "low", "bot"]


def get_reps_bins(lib):
    reps = ALL_REPS[lib]
    res = []
    for rep in reps:
        res += [f"rep{rep}_{sort_bin}" for sort_bin in reps_to_bins(rep)]
    return res


def rep_to_guide_start_seq(rep):
    return "GGAAAGGACGAAACACCG" if rep in [1, 2, 3, 4] else ""


def write_guide_start_seqs_file(path, lib):
    with open(path, "w") as f:
        for r in ALL_REPS[lib]:
            for b in reps_to_bins(r):
                f.write(f"rep{r}_{b},{rep_to_guide_start_seq(r)}\n")


if __name__ == "__main__":
    lib = sys.argv[1]
    path = sys.argv[2]
    write_guide_start_seqs_file(path, lib)
