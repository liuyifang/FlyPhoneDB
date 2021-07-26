#!/bin/bash

Rscript --vanilla FlyPhone_parallel_batch.R \
		--matrix ./test_dataset_60M/2021-03-25_matrix.csv \
		--metadata ./test_dataset_60M/2021-03-25_metadata.csv \
		--lrpair ./annotation/Ligand_receptor_pair_high_confident_2021vs1_clean.txt \
		--corecomponents ./annotation/Pathway_core_components_2021vs1_clean.txt \
		--cores 8 \
		--output ./output_test_dataset_60M/
