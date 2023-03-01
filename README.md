# anbe_manuscript
Scripts for plots used in manuscript. Follows (WorkflowHub standards)[https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#uploading-workflows-to-workflowhub].
```
├── .gitignore
├── README.md
├── LICENSE.md
├── CODE_OF_CONDUCT.md
├── CONTRIBUTING.md
├── .tests
│   ├── integration
│   └── unit
├── images
│   └── rulegraph.svg
├── workflow
│   ├── rules
|   │   ├── map_collect.smk
|   │   ├── filter_annotate.smk
|   │   ├── run_models.smk
|   │   └── generate_figures.smk
│   ├── envs
|   │   ├── map_collect.yaml
|   │   └── analysis.yaml
│   ├── scripts
|   │   └── ...
│   ├── notebooks
│   ├── report
|   │   ├── plot1.rst
|   │   ├── plot2.rst
|   │   ├── plot3.rst
|   │   ├── plot4.rst
|   │   └── plot5.rst
│   ├── Snakefile
|   └── documentation.md
├── config
├── results
│   ├── raw/
|   ├── mapped/
|   ├── collected/
|   ├── filtered_annotated/
|   └── ...
└── resources
    ├── library/
    └── ...

```
