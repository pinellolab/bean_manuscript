import sys
import bean as be
import pandas as pd

bdata_path = sys.argv[1]
out_bdata_path = sys.argv[2]
bdata = be.read_h5ad(bdata_path)
if "LDLvar" in bdata_path:
    guide_info_with_acc = pd.read_csv(
        "resources/gRNA_info/LDLvar_gRNA_bean_accessibility.csv", index_col=0
    )
    bdata.guides = pd.concat(
        [
            bdata.guides,
            guide_info_with_acc.loc[bdata.guides.index, ["genomic_pos", "chr"]],
        ],
        axis=1,
    )
    bdata.guides = bdata.guides.rename(
        columns={
            "target_pos": "spacer_target_pos",
            "Target gene/variant": "target_variant",
            "Group": "target_group",
            "Group2": "target_group2",
        }
    )
else:
    if "targetPos" in bdata.guides.columns:
        bdata.guides = bdata.guides.rename(columns={"targetPos": "genomic_pos"})
        assert "genomic_pos" in bdata.guides.columns
        bdata.guides["chr"] = "chr19"
bdata.get_guide_edit_rate()
lq_map = {"top": 0.8, "high": 0.6, "bulk": 0.0, "low": 0.2, "bot": 0.0}
uq_map = {"top": 1.0, "high": 0.8, "bulk": 1.0, "low": 0.4, "bot": 0.2}
bdata.samples["lower_quantile"] = bdata.samples["bin"].map(lq_map)
bdata.samples["upper_quantile"] = bdata.samples["bin"].map(uq_map)
bdata.write(out_bdata_path)
