# Validation commands

Run from the project root after the public release has been staged and
relative paths and repository metadata have been configured:

```text
python release/communications_biology_code/figure_scripts/make_figure1.py
python release/communications_biology_code/figure_scripts/make_figure2.py
python release/communications_biology_code/figure_scripts/make_figure3.py
python release/communications_biology_code/figure_scripts/make_figure4.py
python release/communications_biology_code/figure_scripts/make_figure5.py
python release/communications_biology_code/figure_scripts/make_supplementary.py
```

Then compare the generated number and source maps with the frozen QC files. The release does not include restricted inputs and cannot reproduce the full statistical analyses without the original provider files.
