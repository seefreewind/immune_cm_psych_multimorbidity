#!/usr/bin/env python3
from make_figures import configure_style, make_supplementary, write_qc_files

if __name__ == "__main__":
    configure_style()
    make_supplementary()
    write_qc_files()
