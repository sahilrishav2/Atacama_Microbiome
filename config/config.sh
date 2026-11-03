#!/usr/bin/env bash
# Central configuration. Edit here; every script in scripts/ sources this file.
# Every variable can also be overridden on the command line, e.g.:
#   THREADS=8 bash scripts/05_assembly.sh
#   RAW_DIR=/scratch/raw bash scripts/02_quality_control.sh
# shellcheck disable=SC2034   # variables are consumed by the scripts that source this file

# The scripts use arrays and mapfile, which need bash >= 4 (macOS ships 3.2: 'brew install bash').
if (( BASH_VERSINFO[0] < 4 )); then
  echo "ERROR: bash >= 4 is required (found ${BASH_VERSION})" >&2
  return 1 2>/dev/null || exit 1
fi

# ---------------- Resources ----------------
# Defaults follow the machine you run on; set THREADS / SORT_MEM yourself for clusters.
detect_cpus() { nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4; }
THREADS="${THREADS:-$(detect_cpus)}"
FASTP_THREADS="${FASTP_THREADS:-$(( THREADS > 16 ? 16 : THREADS ))}"   # fastp uses at most 16 threads
CHECKM2_THREADS="${CHECKM2_THREADS:-$THREADS}"
DIAMOND_THREADS="${DIAMOND_THREADS:-$THREADS}"
HUMANN_THREADS="${HUMANN_THREADS:-$THREADS}"
SORT_MEM="${SORT_MEM:-8G}"               # fasterq-dump -m (memory for sorting)

# ---------------- Input ----------------
ACCESSIONS_FILE="${ACCESSIONS_FILE:-config/accessions.txt}"

# ---------------- Data directories (large, git-ignored) ----------------
SRA_DIR="${SRA_DIR:-data/sra}"
RAW_DIR="${RAW_DIR:-data/raw}"
TRIM_DIR="${TRIM_DIR:-data/trimmed}"
TMP_DIR="${TMP_DIR:-data/tmp}"
HUMANN_INPUT_DIR="${HUMANN_INPUT_DIR:-data/humann_input}"   # read pairs concatenated for HUMAnN (temporary)

# ---------------- Result directories (numbered by script) ----------------
QC_RAW_DIR="${QC_RAW_DIR:-results/02_qc/raw}"
QC_TRIM_DIR="${QC_TRIM_DIR:-results/02_qc/trimmed}"
FASTP_REPORT_DIR="${FASTP_REPORT_DIR:-results/02_qc/fastp}"
METAPHLAN_DIR="${METAPHLAN_DIR:-results/03_metaphlan}"
HUMANN_DIR="${HUMANN_DIR:-results/04_humann}"
MEGAHIT_DIR="${MEGAHIT_DIR:-results/05_assembly/megahit}"      # raw MEGAHIT output folders
ASSEMBLY_DIR="${ASSEMBLY_DIR:-results/05_assembly/contigs}"    # <sample>.fa (final.contigs.fa copied here)
BINNING_DIR="${BINNING_DIR:-results/06_binning}"
CHECKM2_DIR="${CHECKM2_DIR:-results/07_checkm2}"
VIRAL_HITS_DIR="${VIRAL_HITS_DIR:-results/09_viral_hits}"
STRAINSCAN_OUT_DIR="${STRAINSCAN_OUT_DIR:-results/10_strainscan}"           # StrainScan database (step 10)
STRAINSCAN_ID_DIR="${STRAINSCAN_ID_DIR:-results/11_strainscan_identify}"   # StrainScan strain identification (step 11)
STRAINGE_OUT_DIR="${STRAINGE_OUT_DIR:-results/12_strainge_db}"             # StrainGE database (step 12)
STRAINGE_ID_DIR="${STRAINGE_ID_DIR:-results/13_strainge_identify}"         # StrainGE strain identification (step 13)
CORR_DIR="${CORR_DIR:-results/15_correlation}"

# ---------------- MetaPhlAn (script 03) ----------------
# HUMAnN 3.x only accepts MetaPhlAn profiles made with the vJun23 (or MetaPhlAn 3) marker database,
# whereas current MetaPhlAn releases default to a newer one. The index is therefore pinned here.
# Install it once with:  metaphlan --install --index "$METAPHLAN_INDEX"
METAPHLAN_INDEX="${METAPHLAN_INDEX:-mpa_vJun23_CHOCOPhlAnSGB_202403}"
METAPHLAN_DB_DIR="${METAPHLAN_DB_DIR:-}"        # empty = MetaPhlAn's default database location

# ---------------- HUMAnN (script 04) ----------------
HUMANN_DB_DIR="${HUMANN_DB_DIR:-databases/humann}"
HUMANN_NUC_DB="${HUMANN_NUC_DB:-$HUMANN_DB_DIR/chocophlan}"          # ChocoPhlAn pangenomes
HUMANN_PROT_DB="${HUMANN_PROT_DB:-$HUMANN_DB_DIR/uniref}"            # UniRef90 (DIAMOND)
HUMANN_UTIL_DB="${HUMANN_UTIL_DB:-$HUMANN_DB_DIR/utility_mapping}"   # only needed for HUMANN_REGROUP
HUMANN_REGROUP="${HUMANN_REGROUP:-}"            # optional, e.g. "uniref90_ko uniref90_level4ec"
HUMANN_KEEP_TEMP="${HUMANN_KEEP_TEMP:-0}"       # 1 = keep HUMAnN's (large) temporary folders
HUMANN_KEEP_INPUT="${HUMANN_KEEP_INPUT:-0}"     # 1 = keep the concatenated read files
HUMANN_SKIP_PROFILE_CHECK="${HUMANN_SKIP_PROFILE_CHECK:-0}"   # 1 = do not check the MetaPhlAn database version

# ---------------- Bin filtering for the viral search ----------------
MIN_COMPLETENESS="${MIN_COMPLETENESS:-50}"      # percent, bins must be >= this
MAX_CONTAMINATION="${MAX_CONTAMINATION:-5}"     # percent, bins must be <= this

# ---------------- Viral protein database (script 08) ----------------
VIRAL_DB_DIR="${VIRAL_DB_DIR:-databases/ncbi_viral}"
VIRAL_FAA="${VIRAL_FAA:-$VIRAL_DB_DIR/viral.protein.faa}"
VIRAL_DMND="${VIRAL_DMND:-$VIRAL_DB_DIR/diamond_viral_db}"

# ---------------- StrainScan (script 10) ----------------
STRAINSCAN_DIR="${STRAINSCAN_DIR:-$HOME/StrainScan}"                       # where you cloned StrainScan
STRAIN_GENOMES_DIR="${STRAIN_GENOMES_DIR:-data/strain_genomes/C.acnes}"    # folder of genome FASTAs
STRAIN_DB_NAME="${STRAIN_DB_NAME:-C.acnes}"

# ---------------- StrainGE (steps 12 and 13) ----------------
# Install StrainGE with conda (envs/strainge.yml). prepare_strainge_db.py is a separate download from the
# StrainGE database tutorial: save it in STRAINGE_DIR.
STRAINGE_DIR="${STRAINGE_DIR:-$HOME/StrainGE}"
# Folder with the reference genomes (*.fa) AND the assembly report (*.txt), as prepare_strainge_db.py expects
STRAINGE_GENOMES_DIR="${STRAINGE_GENOMES_DIR:-data/strainge_genomes/C.acnes}"
STRAINGE_DB_NAME="${STRAINGE_DB_NAME:-C.acnes}"
STRAINGE_KMER="${STRAINGE_KMER:-23}"                       # k-mer size, used for the database AND for the reads
STRAINGE_CLUSTER_ARGS="${STRAINGE_CLUSTER_ARGS:--d -C 0.99 -c 0.90}"   # options of 'straingst cluster'

# ---------------- Strain identification (steps 11 and 13) ----------------
# Read pairs are looked up in STRAIN_FASTQ_DIR (default: the trimmed reads of step 2). Any of the names
# <name>_1/_2, <name>_R1/_R2, <name>_R1_001/_R2_001 (.fastq/.fq, optionally .gz, optionally .trimmed) work,
# so the identification steps can also be used on your own FASTQ files.
STRAIN_FASTQ_DIR="${STRAIN_FASTQ_DIR:-$TRIM_DIR}"
# Databases to search, one per line (blank lines and # comments are ignored). If the list is empty or the file
# is missing, the database built by step 10 (StrainScan) / step 12 (StrainGE) is used.
#   StrainScan: path of the database folder           e.g. /full_path/E.meliloti
#   StrainGE:   database path WITHOUT the suffix      e.g. /full_path/E.meliloti  (= /full_path/E.meliloti_pan-genome-db.hdf5)
STRAINSCAN_DB_LIST="${STRAINSCAN_DB_LIST:-config/strainscan_databases.txt}"
STRAINGE_DB_LIST="${STRAINGE_DB_LIST:-config/strainge_databases.txt}"
# Optional label added to the output names (e.g. 100X_110S), to keep runs with different settings apart
STRAIN_OUTPUT_SUFFIX="${STRAIN_OUTPUT_SUFFIX:-}"
# Run details (tool output, GNU time statistics) are appended to processing.log in each results folder;
# override with STRAINSCAN_LOG / STRAINGE_LOG.

# ---------------- Bacteriophage abundance (script 14) ----------------
# DIAMOND blastx of the trimmed reads against the viral proteins of step 8; counts per phage -> TPM.
PHAGE_DIR="${PHAGE_DIR:-results/14_phage_abundance}"
PHAGE_EVALUE="${PHAGE_EVALUE:-1e-5}"               # DIAMOND e-value cut-off
PHAGE_MIN_IDENTITY="${PHAGE_MIN_IDENTITY:-30}"     # minimum % identity of a read's best hit
PHAGE_MIN_FRAGMENTS="${PHAGE_MIN_FRAGMENTS:-10}"   # a phage needs at least this many read pairs over all samples
PHAGE_NAME_FILTER="${PHAGE_NAME_FILTER:-phage}"    # regex an organism name must match to count as a phage
# Table used by the correlation step: one phage per row, one sample per column (CSV). Override it to use your own.
PHAGE_TPM="${PHAGE_TPM:-$PHAGE_DIR/bacteriophage_abundance_tpm.csv}"

# ---------------- Correlation analysis (script 15, R) ----------------
# The bacteriophage abundance table comes from step 14 (PHAGE_TPM above); point PHAGE_TPM at your own CSV to use that instead.
CORR_SPECIES="${CORR_SPECIES:-Cutibacterium_acnes}"   # species name as in MetaPhlAn, without the s__ prefix
CORR_LABEL="${CORR_LABEL:-Cutibacterium acnes}"       # label shown in the plot
CORR_TOP_PATHWAYS="${CORR_TOP_PATHWAYS:-20}"          # how many of the most abundant pathways to correlate

# ---------------- Downstream community analysis (R Markdown documents in the repository root) ----------------
# Microbial_community_analysis.Rmd and Co-occurrence_network.Rmd work on amplicon and on shotgun data.
DOWNSTREAM_DIR="${DOWNSTREAM_DIR:-downstream_input}"          # count tables + metadata files, see downstream_input/README.md
DOWNSTREAM_OUT_DIR="${DOWNSTREAM_OUT_DIR:-results/downstream}"  # rendered reports, figures, network.graphml
DOWNSTREAM_RANK="${DOWNSTREAM_RANK:-species}"                 # taxonomic rank taken from MetaPhlAn (prepare_downstream_input.sh)

# ---------------- 16S amplicon workflow (scripts/amplicon_*.sh; documents in Amplicon_data/) ----------------
AMPLICON_DIR="${AMPLICON_DIR:-results/amplicon}"                      # results: <dataset>/, genus_tables/, merged tables, reports/
AMPLICON_METADATA="${AMPLICON_METADATA:-$DOWNSTREAM_DIR/metadata.csv}"  # columns Sample and Origin, covers all datasets
SILVA_DB="${SILVA_DB:-databases/silva/silva_nr99_v138.1_wSpecies_train_set.fa.gz}"   # taxonomy reference for DADA2
SILVA_URL="${SILVA_URL:-https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz?download=1}"
AMPLICON_COLLAPSE="${AMPLICON_COLLAPSE:-mean}"           # merging: how rows of one genus are combined (mean as originally, or sum)
AMPLICON_PREVALENCE="${AMPLICON_PREVALENCE:-0.9}"        # filtering: a genus must be present in this fraction of the samples
AMPLICON_GROUP_COLUMN="${AMPLICON_GROUP_COLUMN:-Origin}" # differential abundance: metadata column with the groups
AMPLICON_GROUP_A="${AMPLICON_GROUP_A:-Plant}"            # contrast: group A - group B
AMPLICON_GROUP_B="${AMPLICON_GROUP_B:-Soil}"
AMPLICON_FDR="${AMPLICON_FDR:-0.05}"                     # adjusted p-value cut-off for the table of significant genera
DOWNSTREAM_SOURCE="${DOWNSTREAM_SOURCE:-shotgun}"        # prepare_downstream_input.sh: shotgun (MetaPhlAn) or amplicon

# ---------------- Helpers ----------------
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }

die() { echo "ERROR: $*" >&2; exit 1; }

# require_cmd CMD...  : stop with a hint if a tool is not on PATH (usually the wrong conda env)
require_cmd() {
  local missing=() c
  for c in "$@"; do
    command -v "$c" > /dev/null 2>&1 || missing+=("$c")
  done
  if (( ${#missing[@]} > 0 )); then
    echo "ERROR: command(s) not found: ${missing[*]}" >&2
    echo "       Activate the conda environment named in the header of this script (see README)." >&2
    exit 1
  fi
}

# have_files FILE...  : true if every argument is an existing, non-empty file
have_files() {
  local f
  for f in "$@"; do
    [[ -s "$f" ]] || return 1
  done
}

# Fills the SAMPLES array from $ACCESSIONS_FILE (ignores blank lines and # comments)
load_samples() {
  if [[ ! -s "$ACCESSIONS_FILE" ]]; then
    echo "ERROR: accession file '$ACCESSIONS_FILE' is missing or empty" >&2
    exit 1
  fi
  mapfile -t SAMPLES < <(sed -e 's/#.*//' -e 's/[[:space:]]//g' -e '/^$/d' "$ACCESSIONS_FILE")
  if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "ERROR: no accessions found in '$ACCESSIONS_FILE'" >&2
    exit 1
  fi
}

# load_db_list LIST_FILE FALLBACK : fills the DBS array with the entries of LIST_FILE (blank lines and
# # comments ignored). If the file is missing or has no entries, FALLBACK is used.
load_db_list() {
  local list="$1" fallback="$2"
  DBS=()
  if [[ -f "$list" ]]; then
    mapfile -t DBS < <(sed -e 's/#.*//' -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e '/^$/d' "$list")
  fi
  if [[ ${#DBS[@]} -eq 0 ]]; then
    log "No databases listed in '$list'; using the database built by the matching build step: $fallback"
    DBS=("$fallback")
  fi
}

# load_read_pairs DIR : finds read pairs in DIR (searched recursively) and fills the parallel arrays
# PAIR_NAMES, PAIR_R1 and PAIR_R2. R1 files are <name>_R1_001 / <name>_R1 / <name>_1 followed by
# [.trimmed].fastq|fq[.gz]; the R2 file is the same name with 1 -> 2. R1 files without an R2 are reported and skipped.
load_read_pairs() {
  local dir="$1" r1 base r2 tag t1 t2 re
  PAIR_NAMES=(); PAIR_R1=(); PAIR_R2=()
  [[ -d "$dir" ]] || die "read directory '$dir' does not exist"
  while IFS= read -r -d '' r1; do
    base="$(basename "$r1")"
    for tag in "_R1_001:_R2_001" "_R1:_R2" "_1:_2"; do
      t1="${tag%%:*}"; t2="${tag##*:}"
      re="^(.+)${t1}(\\.trimmed)?\\.(fastq|fq)(\\.gz)?\$"
      if [[ "$base" =~ $re ]]; then
        r2="$(dirname "$r1")/${BASH_REMATCH[1]}${t2}${BASH_REMATCH[2]}.${BASH_REMATCH[3]}${BASH_REMATCH[4]}"
        if [[ -s "$r2" ]]; then
          PAIR_NAMES+=("${BASH_REMATCH[1]}"); PAIR_R1+=("$r1"); PAIR_R2+=("$r2")
        else
          log "WARNING: no R2 file found for $r1, skipping"
        fi
        break
      fi
    done
  done < <(find "$dir" -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) -print0 | sort -z)
}

# usage_table_init FILE : creates the per-run resource table if it does not exist yet
usage_table_init() {
  [[ -s "$1" ]] || printf 'tool\tsample\tdatabase\tstep\tuser_s\tsystem_s\twall_s\tmax_rss_kb\texit_code\n' > "$1"
}

# run_logged LOG USAGE_TSV TOOL SAMPLE DATABASE STEP -- COMMAND...
# Runs COMMAND with its output appended to LOG and one line of resource use (user/system/wall seconds, peak
# memory in kB; GNU time is used when /usr/bin/time exists) appended to USAGE_TSV. Returns the command's exit code.
run_logged() {
  local logf="$1" usage="$2" tool="$3" sample="$4" db="$5" step="$6"
  shift 7
  local tfile rc=0 user=NA sys=NA wall=NA rss=NA
  tfile="$(mktemp)"
  { echo "### $(date '+%F %T') | $tool | $sample | $db | $step"; echo "### $*"; } >> "$logf"
  if [[ -x /usr/bin/time ]] && /usr/bin/time -v true > /dev/null 2>&1; then
    /usr/bin/time -v -o "$tfile" "$@" >> "$logf" 2>&1 || rc=$?
    cat "$tfile" >> "$logf"
    IFS=$'\t' read -r user sys wall rss < <(awk '
      /User time \(seconds\)/            { sub(/.*: /, ""); u = $0 }
      /System time \(seconds\)/          { sub(/.*: /, ""); s = $0 }
      /Elapsed \(wall clock\) time/      { sub(/.*: /, ""); n = split($0, a, ":"); w = (n == 3) ? a[1]*3600 + a[2]*60 + a[3] : a[1]*60 + a[2] }
      /Maximum resident set size/        { sub(/.*: /, ""); r = $0 }
      END { printf "%s\t%s\t%s\t%s\n", (u == "" ? "NA" : u), (s == "" ? "NA" : s), (w == "" ? "NA" : w), (r == "" ? "NA" : r) }' "$tfile")
  else
    "$@" >> "$logf" 2>&1 || rc=$?
  fi
  rm -f "$tfile"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$tool" "$sample" "$db" "$step" "$user" "$sys" "$wall" "$rss" "$rc" >> "$usage"
  return "$rc"
}

# summarize_usage USAGE_TSV LOG : prints totals over all recorded runs and appends them to LOG
summarize_usage() {
  [[ -s "$1" ]] || return 0
  {
    echo "========================================"
    echo "RESOURCE SUMMARY (all runs recorded in $1)"
    awk -F'\t' 'NR > 1 {
        n++; if ($9 != 0) f++
        if ($5 != "NA") u += $5
        if ($6 != "NA") s += $6
        if ($7 != "NA") w += $7
        if ($8 != "NA") { m += $8; k++; if ($8 + 0 > mx) mx = $8 + 0 }
      }
      END {
        printf "Runs recorded: %d (failed: %d)\n", n, f
        printf "Total user time (seconds): %.2f\n", u
        printf "Total system time (seconds): %.2f\n", s
        printf "Total wall-clock time (seconds): %.2f\n", w
        printf "Peak memory of a single run (kbytes): %d\n", mx
        printf "Average peak memory per run (kbytes): %.0f\n", (k ? m / k : 0)
      }' "$1"
  } | tee -a "$2" >&2
}
