#!/bin/bash
#SBATCH --job-name=git_update
#SBATCH --output=git_update-%j.out
#SBATCH --error=git_update-%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH -n 8
#SBATCH --mem=16G

REPO="/uoa/home/r02hw22/sharedscratch/repos/PhD---omics-analysis-portfolio"
COMMIT_MSG="${1:-Update repository files}"

cd "$REPO"

echo "Repository: $REPO"
echo "Branch: $(git branch --show-current)"

# Refuse to proceed if another Git operation still owns/left a lock.
if [[ -e .git/index.lock ]]; then
    echo "ERROR: .git/index.lock exists."
    echo "Confirm no Git process is running, remove the stale lock, and resubmit."
    exit 1
fi

# Prevent a batch job from hanging while waiting for interactive HTTPS credentials.
export GIT_TERMINAL_PROMPT=0

echo "Staging changes..."
git add -A

echo "Staged changes:"
git diff --cached --stat

if git diff --cached --quiet; then
    echo "No staged changes to commit."
    exit 0
fi

echo "Creating commit..."
git commit -m "$COMMIT_MSG"

BRANCH=$(git branch --show-current)
echo "Pushing branch: $BRANCH"

# For SSH remotes, avoid hanging on an interactive password/passphrase prompt.
REMOTE=$(git remote get-url origin)
if [[ "$REMOTE" == git@* || "$REMOTE" == ssh://* ]]; then
    export GIT_SSH_COMMAND="ssh -o BatchMode=yes"
fi

git push origin "$BRANCH"

echo "Completed successfully."
