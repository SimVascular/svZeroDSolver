#!/usr/bin/env bash
set -euo pipefail

: "${DEFAULT_BRANCH:?DEFAULT_BRANCH is required}"
: "${SOURCE_SHA:?SOURCE_SHA is required}"
: "${CURRENT_VERSION:?CURRENT_VERSION is required}"
: "${RELEASE_VERSION:?RELEASE_VERSION is required}"
: "${GITHUB_OUTPUT:?GITHUB_OUTPUT is required}"

abort() {
  printf '%s\n' "$*" >&2
  exit 1
}

[[ $(git rev-parse HEAD) == "$SOURCE_SHA" ]] || abort "HEAD must match SOURCE_SHA."
git diff --quiet --ignore-submodules=none && git diff --cached --quiet --ignore-submodules=none \
  || abort "Release recording requires a clean tracked worktree and index."

python .github/scripts/release_version.py set \
  --expected-current "$CURRENT_VERSION" --version "$RELEASE_VERSION"
git add -- pyproject.toml
expected_tree=$(git write-tree)
tag="v$RELEASE_VERSION"
branch_ref="refs/remotes/origin/$DEFAULT_BRANCH"
git fetch --no-tags origin "refs/heads/$DEFAULT_BRANCH:$branch_ref"
branch_sha=$(git rev-parse "$branch_ref")

# A retry may follow a successful atomic push and an interrupted release step.
if git ls-remote --exit-code --refs origin "refs/tags/$tag" >/dev/null; then
  git fetch --no-tags origin "refs/tags/$tag:refs/tags/$tag"
  commit_sha=$(git rev-parse "$tag^{commit}")
  if [[ $(git rev-list --parents -n 1 "$commit_sha") != "$commit_sha $SOURCE_SHA" ]] \
    || [[ $(git rev-parse "$commit_sha^{tree}") != "$expected_tree" ]] \
    || ! git merge-base --is-ancestor "$commit_sha" "$branch_sha"; then
    abort "Existing $tag does not match this release on $DEFAULT_BRANCH; reconcile it before retrying."
  fi
else
  status=$?
  [[ $status == 2 ]] || abort "Could not check the remote release tag."
  [[ $branch_sha == "$SOURCE_SHA" ]] \
    || abort "$DEFAULT_BRANCH advanced after this release started; reconcile the uploaded version before retrying."

  GIT_AUTHOR_NAME="Zachary Sexton" GIT_AUTHOR_EMAIL="zsexton@stanford.edu" \
    GIT_COMMITTER_NAME="Zachary Sexton" GIT_COMMITTER_EMAIL="zsexton@stanford.edu" \
    git -c commit.gpgSign=false commit -m "Release pysvzerod $RELEASE_VERSION"
  commit_sha=$(git rev-parse HEAD)
  git tag "$tag" "$commit_sha"
  if ! git push --atomic origin "HEAD:refs/heads/$DEFAULT_BRANCH" "refs/tags/$tag"; then
    abort "The release commit and tag could not be recorded; reconcile the remote branch and uploaded version before retrying."
  fi
fi

printf 'commit_sha=%s\n' "$commit_sha" >> "$GITHUB_OUTPUT"
