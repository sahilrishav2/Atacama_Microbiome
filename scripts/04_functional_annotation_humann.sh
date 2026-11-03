#!/usr/bin/env bash
# Step 4 - Functional annotation with HUMAnN, using the MetaPhlAn profiles from step 3 as the
# taxonomic profile (--taxonomic-profile). HUMAnN therefore does not run MetaPhlAn again: it builds a
# sample-specific pangenome database from the species MetaPhlAn found, maps the reads to it, and
# sends the unmapped reads to a translated search against UniRef90.
# Per sample: gene families, pathway abundance and pathway coverage. These are merged across samples
# and converted to copies per million (CPM).
# Env: atacama-humann   (run scripts/setup_humann_databases.sh once)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
load_samples
require_cmd humann humann_join_tables humann_renorm_table
[[ -d "$HUMANN_NUC_DB" && -d "$HUMANN_PROT_DB" ]] \
  || die "HUMAnN databases not found ($HUMANN_NUC_DB, $HUMANN_PROT_DB). Run scripts/setup_humann_databases.sh first."

PER_SAMPLE="$HUMANN_DIR/per_sample"
MERGED="$HUMANN_DIR/merged"
mkdir -p "$PER_SAMPLE" "$MERGED" "$HUMANN_INPUT_DIR"

# HUMAnN 3.x only accepts MetaPhlAn profiles made with the v3 or vJun23 marker database.
# The database name is written in the first lines of every MetaPhlAn profile.
profile_ok() {
  local re='mpa_vJun23|mpa_v3[01]'
  [[ "$(head -n 5 "$1")" =~ $re ]]
}

TEMP_ARGS=()
[[ "$HUMANN_KEEP_TEMP" == "1" ]] || TEMP_ARGS=(--remove-temp-output)

for SAMPLE in "${SAMPLES[@]}"; do
  R1="$TRIM_DIR/${SAMPLE}_1.trimmed.fastq.gz"
  R2="$TRIM_DIR/${SAMPLE}_2.trimmed.fastq.gz"
  PROFILE="$METAPHLAN_DIR/${SAMPLE}_metagenome.txt"
  INPUT="$HUMANN_INPUT_DIR/${SAMPLE}.fastq.gz"

  complete=true
  for TABLE in genefamilies pathabundance pathcoverage; do
    [[ -s "$PER_SAMPLE/${SAMPLE}_${TABLE}.tsv" ]] || complete=false
  done
  if $complete; then
    log "$SAMPLE: HUMAnN output exists, skipping"
    continue
  fi

  have_files "$PROFILE" || { log "$SAMPLE: MetaPhlAn profile $PROFILE not found (run step 3), skipping"; continue; }
  have_files "$R1" "$R2" || { log "$SAMPLE: trimmed reads not found in $TRIM_DIR (run step 2), skipping"; continue; }

  if [[ "$HUMANN_SKIP_PROFILE_CHECK" != "1" ]] && ! profile_ok "$PROFILE"; then
    die "$SAMPLE: $PROFILE was not made with the MetaPhlAn vJun23 (or v3) database that HUMAnN 3.x requires.
       Delete the old profiles in $METAPHLAN_DIR and re-run step 3 with METAPHLAN_INDEX=mpa_vJun23_CHOCOPhlAnSGB_202403,
       or set HUMANN_SKIP_PROFILE_CHECK=1 to skip this check."
  fi

  # HUMAnN takes a single input file, so the read pair is concatenated (the documented way to use paired reads)
  if [[ ! -s "$INPUT" ]]; then
    log "$SAMPLE: concatenating read pair"
    cat "$R1" "$R2" > "$INPUT.partial"
    mv "$INPUT.partial" "$INPUT"
  fi

  log "$SAMPLE: HUMAnN"
  humann --input "$INPUT" --input-format fastq.gz \
    --output "$PER_SAMPLE" --output-basename "$SAMPLE" \
    --taxonomic-profile "$PROFILE" \
    --nucleotide-database "$HUMANN_NUC_DB" --protein-database "$HUMANN_PROT_DB" \
    --threads "$HUMANN_THREADS" ${TEMP_ARGS[@]+"${TEMP_ARGS[@]}"}

  [[ "$HUMANN_KEEP_INPUT" == "1" ]] || rm -f "$INPUT"
done

compgen -G "$PER_SAMPLE/*_pathabundance.tsv" > /dev/null || die "no HUMAnN results found in $PER_SAMPLE"

log "Merging per-sample tables"
for TABLE in genefamilies pathabundance pathcoverage; do
  humann_join_tables --input "$PER_SAMPLE" --output "$MERGED/merged_${TABLE}.tsv" --file_name "$TABLE"
done

log "Normalising to copies per million (CPM)"
humann_renorm_table --input "$MERGED/merged_genefamilies.tsv" \
  --output "$MERGED/merged_genefamilies_cpm.tsv" --units cpm
humann_renorm_table --input "$MERGED/merged_pathabundance.tsv" \
  --output "$MERGED/merged_pathabundance_cpm.tsv" --units cpm

# Optional: regroup UniRef90 gene families into other systems (needs the utility_mapping database)
read -ra REGROUP_LIST <<< "$HUMANN_REGROUP"
for GROUP in ${REGROUP_LIST[@]+"${REGROUP_LIST[@]}"}; do
  require_cmd humann_regroup_table
  log "Regrouping gene families to $GROUP"
  humann_regroup_table --input "$MERGED/merged_genefamilies.tsv" --groups "$GROUP" \
    --output "$MERGED/merged_${GROUP}.tsv"
  humann_renorm_table --input "$MERGED/merged_${GROUP}.tsv" \
    --output "$MERGED/merged_${GROUP}_cpm.tsv" --units cpm
done

log "Done. Merged tables are in $MERGED"
