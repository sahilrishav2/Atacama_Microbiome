# Input of the downstream community analysis

`Microbial_community_analysis.Rmd` (alpha and beta diversity) and `Co-occurrence_network.Rmd` (co-occurrence network)
read the tables listed here. They work on any taxon abundance table, so the same documents are used for **amplicon**
and for **shotgun** data. Put your files in this folder (they are also found in sub-folders, as long as the file name
is unique) or set another folder with the `data_dir` parameter of the documents / `DOWNSTREAM_DIR` of the scripts.

All tables are CSV files with a header row. Taxa are rows, samples are columns. Metadata files have one row per sample.

| File | Used for | Content |
|------|----------|---------|
| `countdata.csv` | alpha diversity | taxa (first column) x all samples |
| `metadata.csv` | alpha diversity, network; also the amplicon filtering and differential abundance documents | **sample IDs in the first column** (header `Sample`), a column `Origin` (the alpha diversity groups; `Plant` and `Soil` for the differential abundance contrast) |
| `countdata_plant.csv` | beta diversity, plant samples | taxa x the plant (rhizosphere) samples |
| `metadata_plant.csv` | beta diversity, plant samples | columns `Sample` and `Plant` (at most 4 plants: the colour list has 4 entries) |
| `countdata_soil.csv` | beta diversity, soil samples | taxa x the soil samples |
| `metadata_soil.csv` | beta diversity, soil samples | columns `Sample` and `Region` (at most 12 regions: the colour list has 12 entries) |
| `relative_abundance.csv` | co-occurrence network | taxa x all samples, relative abundance |

**Sample order.** The columns of every count table must be in the same order as the samples of its metadata file
(`Sample` column, or the first column of `metadata.csv`). The documents stop with a message if they are not.

## Where the tables come from

- **Amplicon data:** `scripts/amplicon_analysis.sh` (see the README, "16S amplicon workflow") makes the filtered count table and
  its relative abundance table. Run `DOWNSTREAM_SOURCE=amplicon bash scripts/prepare_downstream_input.sh` to subset and order
  them for each metadata file you put here. You can also put your own final count tables here.
- **Shotgun data:** `scripts/prepare_downstream_input.sh` (default `DOWNSTREAM_SOURCE=shotgun`) converts the merged MetaPhlAn table of step 3 into these
  tables, in the order of the metadata files you put here (rank: `DOWNSTREAM_RANK`, default species). Only the
  metadata files are yours to provide; sample IDs must be the run accessions of `config/accessions.txt`.
  MetaPhlAn gives relative abundances (percent), not read counts: that is fine for Shannon and Simpson diversity,
  for the compositional PCA and for the network, but not for methods that need counts.

## Try it

`examples/toy_data/downstream/` contains a small synthetic set of all these files (random numbers, formats only):

```bash
mamba activate atacama-downstream
bash scripts/downstream_analysis.sh        # reads downstream_input/
# or, with the example data:
DOWNSTREAM_DIR=examples/toy_data/downstream DOWNSTREAM_OUT_DIR=results/downstream_toy bash scripts/downstream_analysis.sh
```
