import sys
import bean as be

bdata_path = sys.argv[1]
reps = sys.argv[2].split(",")
outfile_path = sys.argv[3]

bdata = be.read_h5ad(bdata_path)
bdata_sub = bdata[:, bdata.samples.rep.isin(reps)]
bdata_sub.write(outfile_path)
