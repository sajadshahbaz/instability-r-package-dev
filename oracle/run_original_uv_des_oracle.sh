#!/usr/bin/env bash
set -euo pipefail

readonly ORACLE_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
readonly REPOSITORY_ROOT="$(cd -- "$ORACLE_SCRIPT_DIR/.." && pwd -P)"
readonly MANIFEST="$ORACLE_SCRIPT_DIR/uv_des_manifest.sha256"
readonly COMPARATOR="$ORACLE_SCRIPT_DIR/compare_original_uv_des.R"
readonly REQUIRED_LOCALE="en_US.UTF-8"

WORKSPACE=""
PRE_HEAD=""
PRE_STATUS_FILE=""
PRE_PROTECTED_HASHES=""
RUN_SUCCEEDED=0

die() {
  printf 'CRITICAL FAILURE: %s\n' "$*" >&2
  exit 1
}

validate_workspace() {
  [[ -n "$WORKSPACE" ]] || die "workspace path is empty"
  [[ "$WORKSPACE" != "/" ]] || die "workspace cannot be /"
  [[ "$WORKSPACE" != "/tmp" ]] || die "workspace cannot be /tmp"
  [[ -d "$WORKSPACE" ]] || die "workspace does not exist"
  local resolved
  resolved="$(cd -- "$WORKSPACE" && pwd -P)"
  [[ "$resolved" == /tmp/nm2-uvdes-oracle.* ]] || die "workspace was not created by the expected mktemp template"
  [[ "$resolved" != "$REPOSITORY_ROOT" ]] || die "workspace cannot be repository root"
  case "$resolved/" in
    "$REPOSITORY_ROOT"/*) die "workspace cannot be inside repository" ;;
  esac
  [[ ! -e "$resolved/04_results" ]] || die "new workspace unexpectedly contains 04_results"
  WORKSPACE="$resolved"
}

protected_paths() {
  awk '{sub(/^[[:xdigit:]]+[[:space:]]+[* ]?/, ""); print}' "$MANIFEST"
}

capture_protected_hashes() {
  local destination="$1"
  (
    cd -- "$REPOSITORY_ROOT"
    protected_paths | while IFS= read -r path; do
      [[ -f "$path" ]] || exit 1
      sha256sum -- "$path"
    done
  ) > "$destination"
}

post_run_protection_check() {
  local check_failed=0
  if [[ -z "$WORKSPACE" || ! -d "$WORKSPACE" ]]; then
    printf 'CRITICAL FAILURE: no valid workspace available for post-run protection check\n' >&2
    return 1
  fi

  local post_head_file="$WORKSPACE/post_git_head.txt"
  local post_status_file="$WORKSPACE/post_git_status_short.txt"
  local post_hashes="$WORKSPACE/post_protected_sha256.txt"

  git -C "$REPOSITORY_ROOT" rev-parse HEAD > "$post_head_file" || check_failed=1
  git -C "$REPOSITORY_ROOT" status --short > "$post_status_file" || check_failed=1
  capture_protected_hashes "$post_hashes" || check_failed=1

  if [[ -n "$PRE_HEAD" ]] && ! cmp -s -- "$WORKSPACE/pre_git_head.txt" "$post_head_file"; then
    printf 'CRITICAL FAILURE: Git HEAD changed during oracle execution\n' >&2
    check_failed=1
  fi
  if [[ -n "$PRE_STATUS_FILE" ]] && ! cmp -s -- "$PRE_STATUS_FILE" "$post_status_file"; then
    printf 'CRITICAL FAILURE: git status --short changed during oracle execution\n' >&2
    check_failed=1
  fi
  if [[ -n "$PRE_PROTECTED_HASHES" ]] && ! cmp -s -- "$PRE_PROTECTED_HASHES" "$post_hashes"; then
    printf 'CRITICAL FAILURE: a scientific source, formal input, or frozen reference changed\n' >&2
    check_failed=1
  fi
  return "$check_failed"
}

on_exit() {
  local original_status=$?
  trap - EXIT
  if [[ -z "$WORKSPACE" || ! -d "$WORKSPACE" ]]; then
    printf 'ORACLE FAILED before temporary workspace creation\n' >&2
    (( original_status != 0 )) && exit "$original_status"
    exit 1
  fi
  local protection_status=0
  post_run_protection_check || protection_status=$?
  if (( original_status != 0 || protection_status != 0 || RUN_SUCCEEDED != 1 )); then
    printf 'ORACLE FAILED; diagnostic workspace preserved: %s\n' "$WORKSPACE" >&2
    exit 1
  fi
  printf 'ORACLE PASSED; workspace preserved: %s\n' "$WORKSPACE"
}
trap on_exit EXIT

command -v git >/dev/null 2>&1 || die "git is unavailable"
command -v sha256sum >/dev/null 2>&1 || die "sha256sum is unavailable"
command -v Rscript >/dev/null 2>&1 || die "Rscript is unavailable"
command -v locale >/dev/null 2>&1 || die "locale is unavailable"
[[ -f "$MANIFEST" ]] || die "manifest is missing"
[[ -f "$COMPARATOR" ]] || die "comparator is missing"

WORKSPACE="$(mktemp -d /tmp/nm2-uvdes-oracle.XXXXXXXXXX)"
validate_workspace

PRE_HEAD="$(git -C "$REPOSITORY_ROOT" rev-parse HEAD)" || die "cannot capture Git HEAD"
PRE_STATUS_FILE="$WORKSPACE/pre_git_status_short.txt"
PRE_PROTECTED_HASHES="$WORKSPACE/pre_protected_sha256.txt"
printf '%s\n' "$PRE_HEAD" > "$WORKSPACE/pre_git_head.txt"
git -C "$REPOSITORY_ROOT" status --short > "$PRE_STATUS_FILE"
capture_protected_hashes "$PRE_PROTECTED_HASHES" || die "cannot checksum protected files"

(
  cd -- "$REPOSITORY_ROOT"
  sha256sum --check --strict "$MANIFEST"
) > "$WORKSPACE/manifest_verification.txt" || die "manifest verification failed"

LC_ALL="$REQUIRED_LOCALE" locale charmap > "$WORKSPACE/locale_charmap.txt" 2>&1 ||
  die "required locale is unavailable: $REQUIRED_LOCALE"

mkdir -p -- "$WORKSPACE/00_config" "$WORKSPACE/01_data/metadata" "$WORKSPACE/01_data/processed"
cp -- "$REPOSITORY_ROOT/00_config/config.yaml" "$WORKSPACE/00_config/config.yaml"
cp -- "$REPOSITORY_ROOT/01_data/metadata/Simplified_Metadata_Table.core.csv" \
  "$WORKSPACE/01_data/metadata/Simplified_Metadata_Table.core.csv"
cp -- "$REPOSITORY_ROOT/01_data/processed/vst_counts_filtered.tsv" \
  "$WORKSPACE/01_data/processed/vst_counts_filtered.tsv"

for relative_path in \
  00_config/config.yaml \
  01_data/metadata/Simplified_Metadata_Table.core.csv \
  01_data/processed/vst_counts_filtered.tsv; do
  source_hash="$(sha256sum -- "$REPOSITORY_ROOT/$relative_path" | awk '{print $1}')"
  copied_hash="$(sha256sum -- "$WORKSPACE/$relative_path" | awk '{print $1}')"
  [[ "$source_hash" == "$copied_hash" ]] || die "copied input hash mismatch: $relative_path"
done

(
  export LC_ALL="$REQUIRED_LOCALE"
  export TZ="Europe/Warsaw"
  export NM2_ROOT="$WORKSPACE"
  cd -- "$WORKSPACE"
  [[ "$(pwd -P)" == "$NM2_ROOT" ]] || exit 1

  Rscript --vanilla -e 'writeLines(capture.output(sessionInfo())); cat("RNGkind:", paste(RNGkind(), collapse=", "), "\n"); cat("Rscript:", Sys.which("Rscript"), "\n"); cat("data.table:", as.character(packageVersion("data.table")), "\n"); cat("ggplot2:", as.character(packageVersion("ggplot2")), "\n"); cat("BLAS:", extSoftVersion()[["BLAS"]], "\n"); cat("LAPACK:", La_version(), "\n")' \
    > "$WORKSPACE/environment.txt"

  Rscript --vanilla "$REPOSITORY_ROOT/02_methods/R/core_pipeline/01_pool_vs_strat_rebuild_uv_des.R" \
    > "$WORKSPACE/stage_01_console.txt" 2>&1
  Rscript --vanilla "$REPOSITORY_ROOT/02_methods/R/core_pipeline/02_plots_pool_vs_strat.R" \
    > "$WORKSPACE/stage_02_console.txt" 2>&1
  Rscript --vanilla "$REPOSITORY_ROOT/02_methods/R/core_pipeline/03_tail_divergence.R" \
    > "$WORKSPACE/stage_03_console.txt" 2>&1
  Rscript --vanilla "$REPOSITORY_ROOT/02_methods/R/core_pipeline/04_make_core_table.R" \
    > "$WORKSPACE/stage_04_console.txt" 2>&1
)

for incidental in \
  04_results/logs/01_pool_vs_strat_log.txt \
  04_results/logs/02_plots_pool_vs_strat_log.txt \
  04_results/logs/03_tail_divergence_log.txt \
  04_results/qc_figures/pool_vs_strat/uv_instability_vs_effect_divergence.pdf \
  04_results/qc_figures/pool_vs_strat/uv_pooled_vs_stratified.pdf \
  04_results/qc_figures/pool_vs_strat/des_instability_vs_effect_divergence.pdf \
  04_results/qc_figures/pool_vs_strat/des_pooled_vs_stratified.pdf \
  04_results/figures/uv_tail_jaccard_curve.png \
  04_results/figures/uv_rank_displacement_hist.png \
  04_results/figures/des_tail_jaccard_curve.png \
  04_results/figures/des_rank_displacement_hist.png; do
  if [[ ! -s "$WORKSPACE/$incidental" ]]; then
    printf 'WARNING: missing or empty incidental output: %s\n' "$incidental" >&2
  fi
done

awk '/04_results\/tables\// {print}' "$MANIFEST" > "$WORKSPACE/reference_table_manifest.sha256"
while read -r expected relative_path; do
  filename="${relative_path##*/}"
  generated="$WORKSPACE/04_results/tables/$filename"
  [[ -s "$generated" ]] || die "missing or empty scientific output: $filename"
  actual="$(sha256sum -- "$generated" | awk '{print $1}')"
  [[ "$actual" == "$expected" ]] || die "SHA-256 mismatch for scientific output: $filename"
  printf 'BYTE IDENTICAL: %s\n' "$filename"
done < "$WORKSPACE/reference_table_manifest.sha256" > "$WORKSPACE/byte_comparison.txt"

Rscript --vanilla "$COMPARATOR" "$REPOSITORY_ROOT" "$WORKSPACE" \
  > "$WORKSPACE/parsed_comparison.txt" 2>&1

for relative_path in \
  00_config/config.yaml \
  01_data/metadata/Simplified_Metadata_Table.core.csv \
  01_data/processed/vst_counts_filtered.tsv; do
  source_hash="$(sha256sum -- "$REPOSITORY_ROOT/$relative_path" | awk '{print $1}')"
  copied_hash="$(sha256sum -- "$WORKSPACE/$relative_path" | awk '{print $1}')"
  [[ "$source_hash" == "$copied_hash" ]] || die "input changed during execution: $relative_path"
done

RUN_SUCCEEDED=1
