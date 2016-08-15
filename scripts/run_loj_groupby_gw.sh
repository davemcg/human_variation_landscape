#!/bin/bash

gzcat 100bp_gene_windows_1bp_slide.20160701_ensembl_homo_sapiens_variation.loj.dat.gz | ~/git/human_variation_landscape/scripts/loj_groupby_gw.py
