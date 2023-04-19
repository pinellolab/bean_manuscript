set -e
bdata_path=$1

mask_col=mask
pids=()

mageck_result_prefix=$2
trailing_args=${@:3}


## MAGeCK input file /PHShome/jr1025/projects/ANBE/mageck_results/
mageck_path=${mageck_result_prefix}/$(basename "${bdata_path}" .h5ad)/
mkdir -p $mageck_path
mageck_outfile=$mageck_path/$(basename "${bdata_path}" .h5ad).mageck_input.txt
sgrna_eff=$mageck_path/$(basename "${bdata_path}" .h5ad).mageck_sgrna_eff.txt
dm_topbot=${mageck_path}/$(basename "${bdata_path}" .h5ad).mageck_dm_topbot.txt
dm_topbot_complete=${mageck_path}/$(basename "${bdata_path}" .h5ad).mageck_dm_topbot_complete.txt
dm_sort=${mageck_path}/$(basename "${bdata_path}" .h5ad).mageck_dm_sort.txt

if [ ! -f $mageck_outfile ] || [ ! -f $dm_topbot ] || [ ! -f $dm_sort ] || [ ! -f $sgrna_eff ]; then
    echo "Writing MAGeCK input files to " $mageck_path " ..."
    mageck_outfile=$mageck_path/$(basename "${bdata_path}" .h5ad).bcmatch.mageck_input.txt
    python scripts/run_models/make_mageck_input.py $bdata_path -p $mageck_path -m $mask_col --use_bcmatch 
fi

## MAGeCK MLE
echo "Submitting MAGeCK runs..."
if [ ! -f $mageck_path/topbot.gene_summary.txt ]; then
    mageck mle -k $mageck_outfile -d $dm_topbot -n $mageck_path/topbot --threads=1 $trailing_args &
    pids+=($!) 
fi
if [ ! -f $mageck_path/sort.gene_summary.txt ]; then
    mageck mle -k $mageck_outfile -d $dm_sort -n $mageck_path/sort --threads=1 $trailing_args &
    pids+=($!) 
fi
if [ ! -f $mageck_path/topbot_var.gene_summary.txt ]; then
    mageck mle -k $mageck_outfile -d $dm_topbot --genes-varmodeling 1000 -n $mageck_path/topbot_var --threads=1 $trailing_args &
    pids+=($!) 
fi
if [ ! -f $mageck_path/sort_var.gene_summary.txt ]; then
    mageck mle -k $mageck_outfile -d $dm_sort --genes-varmodeling 1000 -n $mageck_path/sort_var --threads=1 $trailing_args &
    pids+=($!)
fi

## MAGeCK MLE with EM, fixed eff
mkdir -p $mageck_path/EM/
if [ ! -f $mageck_path/EM/topbot.gene_summary.txt ]; then
    mageck mle -k $mageck_outfile -d $dm_topbot -n $mageck_path/EM/topbot --threads=1 --sgrna-efficiency $sgrna_eff $trailing_args &
    pids+=($!)
fi
if [ ! -f $mageck_path/EM/sort.gene_summary.txt ]; then
    mageck mle -k $mageck_outfile -d $dm_sort -n $mageck_path/EM/sort --threads=1 --sgrna-efficiency $sgrna_eff $trailing_args &
    pids+=($!)
fi
if [ ! -f $mageck_path/EM/topbot_var.gene_summary.txt ]; then
    mageck mle -k $mageck_outfile -d $dm_topbot --genes-varmodeling 1000 -n $mageck_path/EM/topbot_var --threads=1 --sgrna-efficiency $sgrna_eff $trailing_args &
    pids+=($!)
fi
if [ ! -f $mageck_path/EM/sort_var.gene_summary.txt ]; then
    mageck mle -k $mageck_outfile -d $dm_sort --genes-varmodeling 1000 -n $mageck_path/EM/sort_var --threads=1 --sgrna-efficiency $sgrna_eff $trailing_args &
    pids+=($!)
fi

## MAGeCK MLE with EM
mkdir -p $mageck_path/EMf/
if [ ! -f $mageck_path/EMf/topbot.gene_summary.txt ]; then
    mageck mle -k $mageck_outfile -d $dm_topbot -n $mageck_path/EMf/topbot --threads=1 --sgrna-efficiency $sgrna_eff --update-efficiency $trailing_args &
    pids+=($!)
fi
if [ ! -f $mageck_path/EMf/sort.gene_summary.txt ]; then
    mageck mle -k $mageck_outfile -d $dm_sort -n $mageck_path/EMf/sort --threads=1 --sgrna-efficiency $sgrna_eff --update-efficiency $trailing_args &
    pids+=($!)
fi
if [ ! -f $mageck_path/EMf/topbot_var.gene_summary.txt ]; then
    mageck mle -k $mageck_outfile -d $dm_topbot --genes-varmodeling 1000 -n $mageck_path/EMf/topbot_var --threads=1 --sgrna-efficiency $sgrna_eff --update-efficiency $trailing_args &
    pids+=($!)
fi
if [ ! -f $mageck_path/EMf/sort_var.gene_summary.txt ]; then
    mageck mle -k $mageck_outfile -d $dm_sort --genes-varmodeling 1000 -n $mageck_path/EMf/sort_var --threads=1 --sgrna-efficiency $sgrna_eff --update-efficiency $trailing_args &
    pids+=($!)
fi

## MAGeCK RRA
bot_samples=$(awk -F'\t' '{if ($3=='0') {print $1}}' $dm_topbot_complete | paste -sd "," -)
top_samples=$(awk -F'\t' '{if ($3=='1') {print $1}}' $dm_topbot_complete | paste -sd "," -)
if [ $rerun ] || [ ! -f "$mageck_path/rra_top.gene_summary.txt" ]; then
    echo "Running MAGeCK RRA..."
    mageck test -k $mageck_outfile -t $top_samples -c $bot_samples --paired -n $mageck_path/rra_top $trailing_args
    pids+=($!)
fi
if [ $rerun ] || [ ! -f "$mageck_path/rra_bot.gene_summary.txt" ]; then
    echo "Running MAGeCK RRA..."
    mageck test -k $mageck_outfile -t $bot_samples -c $top_samples --paired -n $mageck_path/rra_bot $trailing_args
    pids+=($!)
fi

## CRISPhieRmix

## ACE

## BAGEL

for pid in ${pids[*]}; do
    wait $pid
done