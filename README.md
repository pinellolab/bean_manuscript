# bean_manuscript
Scripts for plots used in manuscript. Adapts (WorkflowHub standards)[https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#uploading-workflows-to-workflowhub].
```
├── .gitignore
├── README.md
├── images
│   └── rulegraph.svg
└── workflow
    ├── rules
    │   ├── map_collect.smk
    │   ├── filter_annotate.smk
    │   ├── run_models.smk
    │   └── generate_figures.smk
    |── config
    |── resources
    |   ├── gRNA_info/
    |   ├── LDLR/
    |   ├── LDLvar/
    |   ├── accessibility/
    |   └── ...
    ├── results
    │   ├── raw/
    |   ├── mapped/
    |   ├── collected/
    |   ├── filtered_annotated/
    |   └── model_runs/
    ├── scripts
    │   ├── map_collect
    │   ├── filter_annotate
    │   ├── run_models
    │   └── evaluate_model_runs
    ├── envs
    │   ├── map_collect.yaml
    │   └── analysis.yaml
    ├── notebooks
    ├── report
    │   └── qc_report.html
    ├── Snakefile
    └── documentation.md

```
