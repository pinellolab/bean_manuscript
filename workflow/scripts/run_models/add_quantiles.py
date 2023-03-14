import sys
import bean as be
import pandas as pd

bdata_path = sys.argv[1]
out_bdata_path = sys.argv[2]
bdata = be.read_h5ad(bdata_path)
guide_info_with_acc = pd.read_csv(
    "resources/gRNA_info/LDLvar_gRNA_beret_accessibility.csv", index_col=0
)
bdata.guides = pd.concat(
    [bdata.guides, guide_info_with_acc.loc[bdata.guides.index, ["genomic_pos", "chr"]]],
    axis=1,
)
bdata.get_guide_edit_rate()
lq_map = {"top": 0.8, "high": 0.6, "bulk": 0.0, "low": 0.2, "bot": 0.0}
uq_map = {"top": 1.0, "high": 0.8, "bulk": 1.0, "low": 0.4, "bot": 0.2}
bdata.samples["lower_quantile"] = bdata.samples["bin"].map(lq_map)
bdata.samples["upper_quantile"] = bdata.samples["bin"].map(uq_map)
bdata.write(out_bdata_path)
