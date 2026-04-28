#!/usr/bin/env bash
set -euo pipefail

SHEET=""
DB=""
OUT=""
THREADS=1
DO_JSON=0
DO_SCREEN=0

usage() {
  echo "Usage: $0 --sheet samples.csv --db ref.msh --out outdir [--threads N] [--json] [--screen]"
}

# ---------------- Parse args ----------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sheet)   SHEET="$2"; shift 2;;
    --db)      DB="$2";    shift 2;;
    --out)     OUT="$2";   shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --json)    DO_JSON=1; shift;;
    --screen)  DO_SCREEN=1; shift;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

[[ -z "$SHEET" || -z "$DB" || -z "$OUT" ]] && { usage; exit 2; }

mkdir -p "$OUT"

echo "[mash_wrapper] Starting run..."
echo "[mash_wrapper] DB     : $DB"
echo "[mash_wrapper] Sheet  : $SHEET"
echo "[mash_wrapper] Outdir : $OUT"
echo "[mash_wrapper] Threads: $THREADS"

# ---------------- Run per-sample wrapper ----------------
tail -n +2 "$SHEET" | while IFS=, read -r SAMPLE R1 R2; do
  [[ -z "$SAMPLE" || -z "$R1" ]] && continue
  echo "[mash_wrapper] → Sample: $SAMPLE"

  args=( -rs "$DB" -t "$THREADS" -o "$OUT/${SAMPLE}" )
  (( DO_JSON == 1 )) && args+=( -j )
  (( DO_SCREEN == 1 )) && args+=( -ms )

  if [[ -n "${R2:-}" && "$R2" != " " ]]; then
    mash_wrapper.py "${args[@]}" -r "$R1" "$R2"
  else
    mash_wrapper.py "${args[@]}" -r "$R1"
  fi
done

# ---------------- Post-run summary ----------------
echo "[mash_wrapper] Completed all samples."

# 1️⃣ Count result types
TAB_COUNT=$(find "$OUT" -type f -name '*.tab' | wc -l)
JSON_COUNT=$(find "$OUT" -type f -name '*.json' | wc -l)
PDF_COUNT=$(find "$OUT" -type f -name '*.pdf' | wc -l)

echo "[mash_wrapper] Result summary:"
echo "   - TAB files : $TAB_COUNT"
echo "   - JSON files: $JSON_COUNT"
echo "   - PDF files : $PDF_COUNT"

# 2️⃣ List PDFs if any
if [[ "$PDF_COUNT" -gt 0 ]]; then
  echo "[mash_wrapper] ✅ PDF report(s) detected:"
  find "$OUT" -type f -name '*.pdf' -printf '      %p\n'
else
  echo "[mash_wrapper] ⚠️  No PDF reports detected. If you expected one, check that Python plotting libs are installed (matplotlib, seaborn)."
fi

echo "[mash_wrapper] Done."