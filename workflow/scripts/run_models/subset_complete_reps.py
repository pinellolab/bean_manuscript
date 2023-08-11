import sys
import bean as be

bdata_path = sys.argv[1]
out_bdata_path = sys.argv[2]

bdata = be.read_h5ad(bdata_path)

good_samples = bdata.samples.groupby("rep")["mask"].sum()
whole_reps = good_samples[good_samples == 5].index.tolist()
bdata_complete = bdata[:, bdata.samples.rep.isin(whole_reps)]
bdata_complete.write(out_bdata_path)
