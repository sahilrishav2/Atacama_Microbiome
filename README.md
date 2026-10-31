The aim of this study is to identify the plant-growth promoting microbes of the Atacama desert that can survive harsh environmental conditions. We performed large scale analysis of all the publicly available microbiome datasets (16S amplicon and Whole metagenome) of this region. The codes used for the analysis are provided here.

# Atacama Microbiome: Shotgun Metagenomics Workflow

A step-by-step, resumable shell workflow for shotgun metagenomic data. It goes from raw SRA
downloads to read QC, taxonomic profiling (MetaPhlAn), functional annotation (HUMAnN), assembly,
binning, bin quality assessment, viral protein screening of high-quality bins, strain-level
identification (StrainScan and StrainGE, including custom database construction), and R scripts that
correlate *Cutibacterium acnes* with bacteriophages and bacteriophages with bacterial pathways.

## Workflow

```mermaid
flowchart TD
    A["SRA accessions<br/>config/accessions.txt"] --> B["01 prefetch + fasterq-dump"]
    B --> C["02 FastQC / MultiQC / fastp"]
    C --> D["03 MetaPhlAn<br/>taxonomic profiles"]
    D -->|"--taxonomic-profile"| HU["04 HUMAnN<br/>functional annotation"]
    C --> HU
    C --> E["05 MEGAHIT<br/>assembly"]
    E --> F["06 BWA + MetaBAT2<br/>binning"]
    F --> G["07 CheckM2<br/>bin quality"]
    G -->|"completeness >= 50%<br/>contamination <= 5%"| H["09 DIAMOND blastp<br/>vs viral proteins"]
    V["08 build DIAMOND<br/>viral database"] --> H
    C --> SS["11 StrainScan<br/>strain identification"]
    S["10 build StrainScan<br/>database"] --> SS
    C --> SG["13 StrainGE<br/>strain identification"]
    SD["12 build StrainGE<br/>database"] --> SG
    D --> R["14 R: Spearman correlations"]
    HU --> R
    P["your phage abundance table<br/>(TPM)"] --> R
```

| Step | Script | Main tools | Key output |
|------|--------|-----------|------------|
| 1 | `01_download_sra.sh` | sra-tools, pigz | `data/raw/<run>_{1,2}.fastq.gz` |
| 2 | `02_quality_control.sh` | FastQC, MultiQC, fastp | `data/trimmed/`, `results/02_qc/` |
| 3 | `03_taxonomic_profiling.sh` | MetaPhlAn | `results/03_metaphlan/merged_abundance_table.txt` |
| 4 | `04_functional_annotation_humann.sh` | HUMAnN | `results/04_humann/merged/merged_pathabundance_cpm.tsv` (+ gene families, pathway coverage) |
| 5 | `05_assembly.sh` | MEGAHIT | `results/05_assembly/contigs/<run>.fa` |
| 6 | `06_binning.sh` | BWA, samtools, MetaBAT2 | `results/06_binning/<run>/bins/bin.*.fa` |
| 7 | `07_bin_quality_checkm2.sh` | CheckM2 | `results/07_checkm2/all_samples_quality_report.tsv` |
| 8 | `08_build_viral_diamond_db.sh` | DIAMOND | `databases/ncbi_viral/diamond_viral_db.dmnd` |
| 9 | `09_viral_search_bins.sh` | DIAMOND | `results/09_viral_hits/all_viral_hits.tsv` |
| 10 | `10_strainscan_build_db.sh` | StrainScan | `results/10_strainscan/<db_name>/` |
| 11 | `11_strainscan_identify_strains.sh` | StrainScan | `results/11_strainscan_identify/<sample>/<database>/` |
| 12 | `12_strainge_build_db.sh` | StrainGE | `results/12_strainge_db/<db_name>/<db_name>_pan-genome-db.hdf5` |
| 13 | `13_strainge_identify_strains.sh` | StrainGE | `results/13_strainge_identify/<sample>/<sample>_<database>.tsv` |
| 14 | `14_correlation_analysis.sh` | R (corrplot, ggcorrplot, psych) | `results/14_correlation/*.png`, `*.csv` |

`scripts/setup_humann_databases.sh` is a one-time download used by step 4.

## Requirements

- Linux x86_64 (or WSL2), bash >= 4, and [conda/mamba](https://github.com/conda-forge/miniforge).
- Disk space: the HUMAnN databases are large (tens of GB), plus the MetaPhlAn, CheckM2 and NCBI viral databases and your reads.
- HUMAnN needs a fair amount of memory and time per sample; set `THREADS` to what your machine can give it.

## Installation

Each group of tools has its own environment because their dependencies often conflict.

```bash
git clone https://github.com/sahilrishav2/Atacama_Microbiome.git
cd Atacama_Microbiome

mamba env create -f envs/core.yml        # steps 1, 2, 5, 6, 8, 9
mamba env create -f envs/metaphlan.yml   # step 3
mamba env create -f envs/humann.yml      # step 4
mamba env create -f envs/checkm2.yml     # step 7
mamba env create -f envs/strainge.yml    # steps 12, 13
mamba env create -f envs/r-analysis.yml  # step 14
```

**One-time database setup**

```bash
# MetaPhlAn marker database. The version is pinned in config/config.sh (METAPHLAN_INDEX)
# because HUMAnN 3.x only accepts profiles made with the vJun23 (or MetaPhlAn 3) database.
mamba activate atacama-metaphlan && metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202403

# HUMAnN databases: ChocoPhlAn, UniRef90, utility mapping (large, takes a while)
mamba activate atacama-humann && bash scripts/setup_humann_databases.sh

# CheckM2 database
mamba activate atacama-checkm2 && checkm2 database --download

# NCBI viral proteins for DIAMOND: handled by scripts/08_build_viral_diamond_db.sh
```

To keep databases outside the repository, set `HUMANN_DB_DIR` (HUMAnN) and `METAPHLAN_DB_DIR` (MetaPhlAn,
also add `--bowtie2db <dir>` to `metaphlan --install`) before running.

**StrainScan** (steps 10, 11) is run from a local clone: follow its
[README](https://github.com/liaoherui/StrainScan) to install it and its dependencies, then set
`STRAINSCAN_DIR` in `config/config.sh`.

**StrainGE** (steps 12, 13): `straingst` comes with `envs/strainge.yml`, but `prepare_strainge_db.py` is a
separate download described in the
[StrainGE database tutorial](https://strainge.readthedocs.io/en/latest/create-database.html). Save it in a folder
and set `STRAINGE_DIR` in `config/config.sh` to that folder (default `~/StrainGE`).

## Usage

1. Put your run accessions in `config/accessions.txt` (one per line). Reads must be paired-end.
2. Adjust threads, memory, thresholds, and paths in `config/config.sh`.
3. Run the scripts **from anywhere**, in order, with the right environment active:

```bash
mamba activate atacama-core
bash scripts/01_download_sra.sh
bash scripts/02_quality_control.sh

mamba activate atacama-metaphlan
bash scripts/03_taxonomic_profiling.sh

mamba activate atacama-humann
bash scripts/04_functional_annotation_humann.sh

mamba activate atacama-core
bash scripts/05_assembly.sh
bash scripts/06_binning.sh

mamba activate atacama-checkm2
bash scripts/07_bin_quality_checkm2.sh

mamba activate atacama-core
bash scripts/08_build_viral_diamond_db.sh   # once
bash scripts/09_viral_search_bins.sh

# Strain-level analysis (independent of steps 3-9, see "Strain-level analysis" below)
bash scripts/10_strainscan_build_db.sh      # in the environment where StrainScan runs
bash scripts/11_strainscan_identify_strains.sh
mamba activate atacama-strainge
bash scripts/12_strainge_build_db.sh
bash scripts/13_strainge_identify_strains.sh

mamba activate atacama-r
bash scripts/14_correlation_analysis.sh     # needs your phage table, see "Step 14 inputs"
```

Every setting in `config/config.sh` can be overridden per run without editing it, for example
`THREADS=8 bash scripts/05_assembly.sh` or `RAW_DIR=/scratch/raw bash scripts/02_quality_control.sh`.

All scripts **skip samples that already have output**, so you can re-run after a failure. If a step
crashed halfway, delete that sample's partial output files first. Each script stops with a hint if a
required tool is missing, which usually means the wrong environment is active.

### Step 4: functional annotation

HUMAnN is run on the trimmed reads of each sample with the MetaPhlAn profile from step 3 as its
taxonomic profile, so MetaPhlAn is not run a second time. Per sample you get gene families, pathway
abundance and pathway coverage in `results/04_humann/per_sample/`. They are merged across samples and
converted to copies per million (CPM) in `results/04_humann/merged/`. To also regroup gene families
(for example to KEGG orthologs), set `HUMANN_REGROUP`:

```bash
HUMANN_REGROUP="uniref90_ko uniref90_level4ec" bash scripts/04_functional_annotation_humann.sh
```

### Strain-level analysis (steps 10-13)

These steps are independent of steps 3-9. Two tools are supported; each has one step to build a database and
one to identify strains in your samples.

| | Build a database | Identify strains | Database format |
|--|------------------|------------------|-----------------|
| StrainScan | step 10 | step 11 | a folder |
| StrainGE (StrainGST) | step 12 | step 13 | `<name>_pan-genome-db.hdf5` |

**Building a database.** Put the reference genomes of one species in a folder, then:

```bash
# StrainScan: a folder of genome FASTA files (data/strain_genomes/C.acnes, or STRAIN_GENOMES_DIR)
bash scripts/10_strainscan_build_db.sh

# StrainGE: genomes (*.fa) AND the assembly report (*.txt) in the same folder
# (data/strainge_genomes/C.acnes, or STRAINGE_GENOMES_DIR); name set by STRAINGE_DB_NAME
bash scripts/12_strainge_build_db.sh
```

Step 12 runs the standard StrainGE recipe: `prepare_strainge_db.py`, `straingst kmerize` for every genome,
`straingst kmersim --all-vs-all -S jaccard -S subset`, `straingst cluster -d -C 0.99 -c 0.90` (change with
`STRAINGE_CLUSTER_ARGS`) and `straingst createdb`. The k-mer size (`STRAINGE_KMER`, default 23) is used for the
database and for the reads, so the two always match.

**Identifying strains.** Steps 11 and 13 search *every read pair* found in `STRAIN_FASTQ_DIR` (default
`data/trimmed`) against *every database* in a list file. Read pairs can be named `<name>_1/_2`, `<name>_R1/_R2` or
`<name>_R1_001/_R2_001` (`.fastq` or `.fq`, optionally gzipped), so you can also point `STRAIN_FASTQ_DIR` at your own
FASTQ files. The list file has one database per line (blank lines and `#` comments are ignored):

```
# config/strainscan_databases.txt       StrainScan: the database folder
/full_path/E.meliloti
/full_path/A.humicireducens

# config/strainge_databases.txt         StrainGE: the path WITHOUT "_pan-genome-db.hdf5"
/full_path/E.meliloti                   (= /full_path/E.meliloti_pan-genome-db.hdf5)
```

If the list is empty, the database built by step 10 / step 12 is used. Typical calls, including the log location
and, for StrainGE, the k-mer size and an optional label added to the output names:

```bash
STRAIN_FASTQ_DIR=/path/to/fastq_files STRAINSCAN_DB_LIST=list.txt \
STRAINSCAN_LOG=/path/to/processing.log bash scripts/11_strainscan_identify_strains.sh

STRAIN_FASTQ_DIR=/path/to/fastq_files STRAINGE_DB_LIST=list.txt STRAINGE_KMER=23 \
STRAIN_OUTPUT_SUFFIX=100X_110S STRAINGE_LOG=/path/to/processing.log bash scripts/13_strainge_identify_strains.sh
```

Both steps write `processing.log` (the output of every tool run plus GNU `time -v` statistics) and
`resource_usage.tsv` (user, system and wall-clock time and peak memory of every run) into their results folder,
followed by a summary of totals. Finished runs are skipped when you run the step again, and a failed run does not stop
the others; the step reports how many failed and exits with an error at the end. StrainGE also keeps the
k-merized reads (`<sample>.hdf5`) so that a new database can be searched without repeating that step; delete them
to save space.

### Step 14: correlation analyses

Step 14 runs three R scripts from `analysis/`; each can also be used on its own (`--help` lists the options).

| R script | What it does | Output |
|----------|--------------|--------|
| `prepare_correlation_inputs.R` | Takes one species (default *C. acnes*) from the merged MetaPhlAn table and the most abundant unstratified pathways from the merged HUMAnN table | `cutibacterium_abundance.csv`, `abundant_pathways.csv` |
| `correlate_cacnes_phages.R` | Spearman correlation of *C. acnes* with every phage, one-row circle plot (red = negative, blue = positive) | `cacnes_phage_correlations.csv`, `cacnes_phage_corrplot.png` |
| `correlate_phages_pathways.R` | Spearman correlation of every phage with every pathway, square-tile plot with the coefficients printed | `phage_pathway_{rho,pvalue,long}.csv`, `phage_pathway_corrplot.png` |

**Step 14 inputs.** The two MetaPhlAn/HUMAnN tables come from steps 3 and 4. The bacteriophage abundance
table is **not produced by this workflow**; supply it yourself (set `PHAGE_TPM`, default
`data/phage/bacteriophage_abundance_tpm.csv`). All tables are CSV files with feature names in the first
column and one numeric column per sample, and the same sample IDs (the run accessions) in every table:

```
"","SRR9729086","SRR9729087", ...
"phage_01",12.5,0,...
"phage_02",3.1,8.7,...
```

The pathway selection in `prepare_correlation_inputs.R` is a parameter, not a fixed rule: by default the
20 pathways with the highest mean CPM are kept (`CORR_TOP_PATHWAYS`, or `--top-n` / `--min-prevalence` on the
script). To correlate a different species, set `CORR_SPECIES` (MetaPhlAn name without `s__`) and `CORR_LABEL`.

**Try it without running the pipeline.** `examples/toy_data/` holds small *synthetic* tables in the right formats:

```bash
mamba activate atacama-r
METAPHLAN_DIR=examples/toy_data/03_metaphlan HUMANN_DIR=examples/toy_data/04_humann \
PHAGE_TPM=examples/toy_data/bacteriophage_abundance_tpm.csv CORR_DIR=results/14_correlation_toy \
bash scripts/14_correlation_analysis.sh
```

## Repository layout

```
.
├── config/
│   ├── config.sh            # all paths, thread counts, thresholds, database versions
│   ├── accessions.txt       # SRA run accessions to process
│   ├── strainscan_databases.txt   # databases searched by step 11 (optional)
│   └── strainge_databases.txt     # databases searched by step 13 (optional)
├── envs/                    # conda environment files
├── scripts/                 # numbered pipeline steps (01 ... 14) + setup_humann_databases.sh
├── analysis/                # R scripts used by step 14
├── examples/toy_data/       # synthetic example tables (formats only, not real data)
├── data/                    # created at run time, git-ignored
├── results/                 # created at run time, git-ignored
└── databases/               # created at run time, git-ignored
```

## Methods notes

- **Trimming:** fastp with paired-end adapter detection, front and tail quality cutting, and mean quality 20.
- **Taxonomy:** MetaPhlAn on trimmed pairs with `--stat_q 0.15`, eukaryotes and archaea ignored, marker database
  `mpa_vJun23_CHOCOPhlAnSGB_202403`.
- **Functional annotation:** HUMAnN 3.9 on the trimmed read pairs (concatenated, as HUMAnN takes one file) with the
  MetaPhlAn profile as `--taxonomic-profile`, ChocoPhlAn pangenomes and UniRef90. Because step 3 ignores archaea and
  eukaryotes, their reads are not mapped to pangenomes and go straight to the translated UniRef90 search. Gene
  families and pathway abundances are merged and normalised to CPM.
- **Assembly and binning:** per-sample MEGAHIT assemblies. Trimmed reads are mapped back with BWA-MEM,
  contig depths come from `jgi_summarize_bam_contig_depths`, and bins from MetaBAT2.
- **Bin selection:** bins with completeness >= 50% and contamination <= 5% (CheckM2) are screened for
  viral proteins. Change `MIN_COMPLETENESS` / `MAX_CONTAMINATION` in `config/config.sh`.
- **Viral screening:** protein predictions from CheckM2 (`protein_files/`) are searched against NCBI
  RefSeq viral proteins with `diamond blastp --sensitive`. A hit indicates homology to viral proteins
  (e.g. prophage or viral-like genes), not by itself a complete virus, so inspect the hits before drawing conclusions.
- **Strain identification:** StrainScan (`StrainScan.py -i R1 -j R2 -d <database>`) and StrainGE's StrainGST
  (reads k-merized with `straingst kmerize`, then `straingst run` against a pan-genome database). Reads are searched
  against each listed database separately; compare results only between runs with the same database and settings.
- **Correlations:** Spearman rank correlation on samples present in both tables. Features with zero variance are
  removed. For each analysis the p-values are corrected across all tested pairs with Benjamini-Hochberg (`q_value_BH`
  in the output tables); the plots show the correlation coefficients only. *C. acnes* vs phages uses pairwise-complete
  observations and does not plot the bacterium's correlation with itself. With fewer than 10 shared samples the
  scripts warn that the correlations are unreliable, and below 3 they stop.

## Troubleshooting

| Problem | Likely cause and fix |
|---------|----------------------|
| `command(s) not found` | The wrong conda environment is active; see the `Env:` line at the top of the script. |
| Step 4 stops with "was not made with the MetaPhlAn vJun23 (or v3) database" | The step 3 profiles came from another marker database (MetaPhlAn >= 4.2 defaults to vJan25, which HUMAnN 3.9 rejects). Delete the old profiles in `results/03_metaphlan/` and re-run step 3 with the pinned `METAPHLAN_INDEX`. |
| HUMAnN complains about the ChocoPhlAn version | The ChocoPhlAn folder is from an older release. Delete `databases/humann/chocophlan` and re-run `scripts/setup_humann_databases.sh`. |
| Steps 11/13: "no read pairs found" | `STRAIN_FASTQ_DIR` holds no matching pairs. File names must follow `_1/_2`, `_R1/_R2` or `_R1_001/_R2_001`; unpaired R1 files are reported and skipped. |
| Step 12: "no assembly report" or "prepare_strainge_db.py not found" | `STRAINGE_GENOMES_DIR` needs the genomes (`*.fa`) and the assembly report (`*.txt`); download `prepare_strainge_db.py` and set `STRAINGE_DIR` (see Installation). |
| Step 13: "StrainGE database not found" | List the database path without `_pan-genome-db.hdf5` (or give the `.hdf5` file), or build it with step 12. |
| A sample is skipped with "not found (run step N)" | The input of that step is missing for the sample; run the earlier step first. |
| Step 14: "Need at least 3 shared samples" | Sample IDs differ between the tables. For MetaPhlAn names, adjust `--metaphlan-suffix` of `prepare_correlation_inputs.R`. |
| Plot labels are cut off | Set `--width` / `--height` for the R scripts. |

## Tools

[SRA Toolkit](https://github.com/ncbi/sra-tools) ·
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) ·
[MultiQC](https://github.com/MultiQC/MultiQC) ·
[fastp](https://github.com/OpenGene/fastp) ·
[MetaPhlAn](https://github.com/biobakery/MetaPhlAn) ·
[HUMAnN](https://github.com/biobakery/humann) ·
[MEGAHIT](https://github.com/voutcn/megahit) ·
[BWA](https://github.com/lh3/bwa) ·
[samtools](https://www.htslib.org/) ·
[MetaBAT2](https://bitbucket.org/berkeleylab/metabat) ·
[CheckM2](https://github.com/chklovski/CheckM2) ·
[DIAMOND](https://github.com/bbuchfink/diamond) ·
[StrainScan](https://github.com/liaoherui/StrainScan) ·
[StrainGE](https://github.com/broadinstitute/strainge) ·
R packages [corrplot](https://github.com/taiyun/corrplot), [ggcorrplot](https://github.com/kassambara/ggcorrplot) and [psych](https://cran.r-project.org/package=psych)

Please cite the original publications of these tools when you use this workflow, in particular
MetaPhlAn 4 (Blanco-Míguez et al., *Nat Biotechnol* 2023), HUMAnN 3 / bioBakery 3 (Beghini et al., *eLife* 2021),
StrainScan (Liao et al., *Microbiome* 2023) and StrainGE (van Dijk et al., *Genome Biol* 2022).

## License

MIT, see [LICENSE](LICENSE).
