#!/usr/bin/env bash
# Usage: ./setup.sh "/path/to/your/matlab/project" "Your Name"
# Creates ~/Desktop/eeg-meditation-connectivity with only .m files under src/,
# adds README, LICENSE, .gitignore and a size-guard hook, makes the first commit,
# and prints exactly what was committed.
set -euo pipefail

SRC="${1:?Give the folder that holds your .m scripts as the first argument}"
HOLDER="${2:?Give the copyright holder name (for LICENSE) as the second argument}"
KIT="$(cd "$(dirname "$0")" && pwd)"
DEST="$HOME/Desktop/eeg-meditation-connectivity"
LIMIT=$((10 * 1024 * 1024))

command -v git >/dev/null || { echo "git is not installed."; exit 1; }
if [ -z "$(git config user.name)" ] || [ -z "$(git config user.email)" ]; then
  echo "Set your git identity first:"
  echo '  git config --global user.name  "Your Name"'
  echo '  git config --global user.email "you@example.com"'
  exit 1
fi
[ -d "$SRC" ] || { echo "Source folder not found: $SRC"; exit 1; }
if [ -e "$DEST" ] && [ -n "$(ls -A "$DEST" 2>/dev/null)" ]; then
  echo "$DEST already exists and is not empty. Move it aside first."; exit 1
fi
mkdir -p "$DEST/src"

echo "== Copying .m files from $SRC =="
copied=0; skipped=0
while IFS= read -r -d '' f; do
  rel="${f#"$SRC"/}"
  size=$(wc -c < "$f" | tr -d ' ')
  if [ "$size" -gt "$LIMIT" ]; then
    echo "  skip (>10 MB): $rel"; skipped=$((skipped+1)); continue
  fi
  mkdir -p "$DEST/src/$(dirname "$rel")"
  cp "$f" "$DEST/src/$rel"
  copied=$((copied+1))
done < <(find "$SRC" \
  \( -type d \( -name results -o -name results_clean -o -name results_state \
     -o -name data -o -name FULL_BACKUP -o -name .git \) -prune \) \
  -o \( -type f -name '*.m' -print0 \))
echo "  copied: $copied   skipped: $skipped"
[ "$copied" -gt 0 ] || { echo "No .m files found. Check the source path."; exit 1; }

cp "$KIT/templates/README.md"     "$DEST/README.md"
cp "$KIT/templates/gitignore.txt" "$DEST/.gitignore"
sed -e "s/__YEAR__/$(date +%Y)/" -e "s/__HOLDER__/$HOLDER/" \
    "$KIT/templates/LICENSE.txt" > "$DEST/LICENSE"
mkdir -p "$DEST/.githooks"
cp "$KIT/templates/pre-commit.txt" "$DEST/.githooks/pre-commit"
chmod +x "$DEST/.githooks/pre-commit"

cd "$DEST"
git init -q -b main 2>/dev/null || { git init -q; git checkout -q -b main; }
git config core.hooksPath .githooks
git add -A

echo
echo "== Files staged for the first commit =="
git ls-files -s | awk '{print $4}' | while read -r p; do
  printf '%10s  %s\n' "$(wc -c < "$p" | tr -d ' ')" "$p"
done
bad=$(git ls-files | grep -Ev '\.m$|^README\.md$|^LICENSE$|^\.gitignore$|^\.githooks/pre-commit$' || true)
if [ -n "$bad" ]; then
  echo; echo "!! Unexpected files staged — NOT committing:"; echo "$bad"; exit 1
fi

git commit -q -m "Initial commit: EEG meditation connectivity MATLAB pipeline"
echo
echo "== Committed $(git ls-files | wc -l | tr -d ' ') files =="
git log --oneline -1
