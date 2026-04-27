#!/usr/bin/env bash

# -------------------------------------------------------------------
# lowercase_directories.sh
#
# Safely rename Git-tracked directories so that directory names are
# lowercase.
#
# This changes directory names only. It does not rename files.
#
# Usage:
#   bash lowercase_directories.sh --dry-run
#   bash lowercase_directories.sh --apply
#
# Recommended workflow:
#   git checkout main
#   git pull origin main
#   git checkout -b lowercase-directories
#   bash lowercase_directories.sh --dry-run
#   bash lowercase_directories.sh --apply
#   git status
#   git diff --summary
#   git commit -m "Standardise directory names to lowercase"
#   git push -u origin lowercase-directories
# -------------------------------------------------------------------

MODE="${1:-}"

if [[ "$MODE" != "--dry-run" && "$MODE" != "--apply" ]]; then
    echo "Usage: bash lowercase_directories.sh --dry-run"
    echo "   or: bash lowercase_directories.sh --apply"
    exit 1
fi

# Make sure we are inside a Git repository.
REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

echo "Repository root:"
echo "$REPO_ROOT"
echo

# Make Git pay attention to case-only renames.
git config core.ignorecase false

# Refuse to run if there are uncommitted changes.
if [[ -n "$(git status --porcelain)" ]]; then
    echo "ERROR: Your working tree is not clean."
    echo
    echo "Commit, stash, or discard your changes before running this script."
    echo
    git status --short
    exit 1
fi

# Create a temporary file containing all Git-tracked directories.
TMP_DIR_LIST="$(mktemp)"

git ls-files \
    | awk '
        {
            path = $0
            while (path ~ /\//) {
                sub("/[^/]*$", "", path)
                print path
            }
        }
    ' \
    | sort -u \
    | awk '
        {
            depth = gsub("/", "/")
            print depth "\t" $0
        }
    ' \
    | sort -rn \
    | cut -f2- \
    > "$TMP_DIR_LIST"

echo "Directory rename preview:"
echo

changes_found=0

while IFS= read -r dir; do
    lower_dir="$(dirname "$dir")/$(basename "$dir" | tr '[:upper:]' '[:lower:]')"

    # Fix dirname output for top-level directories.
    if [[ "$(dirname "$dir")" == "." ]]; then
        lower_dir="$(basename "$dir" | tr '[:upper:]' '[:lower:]')"
    fi

    if [[ "$dir" != "$lower_dir" ]]; then
        changes_found=1
        echo "$dir  ->  $lower_dir"
    fi
done < "$TMP_DIR_LIST"

echo

if [[ "$changes_found" -eq 0 ]]; then
    echo "No uppercase directory names found."
    rm -f "$TMP_DIR_LIST"
    exit 0
fi

if [[ "$MODE" == "--dry-run" ]]; then
    echo "Dry run only. No files were changed."
    echo
    echo "To apply these changes, run:"
    echo "bash lowercase_directories.sh --apply"
    rm -f "$TMP_DIR_LIST"
    exit 0
fi

echo "Applying directory renames..."
echo

while IFS= read -r dir; do
    # The directory may already have moved because a parent directory was renamed.
    if [[ ! -d "$dir" ]]; then
        continue
    fi

    parent="$(dirname "$dir")"
    base="$(basename "$dir")"
    lower_base="$(echo "$base" | tr '[:upper:]' '[:lower:]')"

    if [[ "$base" == "$lower_base" ]]; then
        continue
    fi

    if [[ "$parent" == "." ]]; then
        target="$lower_base"
    else
        target="$parent/$lower_base"
    fi

    # Skip if already correct.
    if [[ "$dir" == "$target" ]]; then
        continue
    fi

    # Use a temporary name so case-only renames work on macOS/Windows too.
    tmp="${dir}.__tmp_lowercase_rename__"

    echo "Renaming: $dir -> $target"

    git mv "$dir" "$tmp"
    git mv "$tmp" "$target"

done < "$TMP_DIR_LIST"

rm -f "$TMP_DIR_LIST"

echo
echo "Done."
echo
echo "Review changes with:"
echo "git status"
echo "git diff --summary"