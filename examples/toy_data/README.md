# Toy data (synthetic)

Random numbers made up to show the **file formats** and to let you test the correlation scripts in
minutes, without running the pipeline. Nothing here is real data and the figures it produces
mean nothing. It deliberately contains one all-zero phage (`phage_06`) to show how constant
features are dropped.

| File | What it imitates |
|------|------------------|
| `03_metaphlan/merged_abundance_table.txt` | output of step 3 (merged MetaPhlAn table) |
| `04_humann/merged/merged_pathabundance_cpm.tsv` | output of step 4 (merged HUMAnN pathways, CPM) |
| `bacteriophage_abundance_tpm.csv` | imitates the output of step 14 (or your own table): one row per phage, one column per sample |
| `cutibacterium_abundance.csv`, `abundant_pathways.csv` | what `analysis/prepare_correlation_inputs.R` writes from the two files above |

Run everything (from the repository root, with the `atacama-r` environment active):

```bash
METAPHLAN_DIR=examples/toy_data/03_metaphlan HUMANN_DIR=examples/toy_data/04_humann \
PHAGE_TPM=examples/toy_data/bacteriophage_abundance_tpm.csv CORR_DIR=results/15_correlation_toy \
bash scripts/15_correlation_analysis.sh
```

## Downstream documents: `downstream/`

Synthetic tables in the formats read by `Microbial_community_analysis.Rmd` and `Co-occurrence_network.Rmd`
(24 samples: 12 plant samples of 4 plants and 12 soil samples of 4 regions, 40 taxa, random counts). They only
show the file formats and let you check your R setup; the figures mean nothing. See `downstream_input/README.md`.

```bash
DOWNSTREAM_DIR=examples/toy_data/downstream DOWNSTREAM_OUT_DIR=results/downstream_toy \
bash scripts/downstream_analysis.sh
```

## Amplicon workflow: `amplicon/`

Two synthetic genus tables (`genus_tables/studyA.tsv`, `studyB.tsv`: a `Genus` column, one column per sample, several ASVs
per genus, counts of at least 1 as in the processing document) and the matching `metadata.csv` (`Sample`, `Origin` =
`Plant`/`Soil`). They stand in for the output of `amplicon_process_dataset.sh` and let you run the merging, filtering and
differential abundance documents without any sequencing data (random numbers, the results mean nothing):

```bash
mkdir -p results/amplicon_toy && cp -r examples/toy_data/amplicon/genus_tables results/amplicon_toy/
AMPLICON_DIR=results/amplicon_toy AMPLICON_METADATA=examples/toy_data/amplicon/metadata.csv \
bash scripts/amplicon_analysis.sh
```
