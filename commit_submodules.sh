#!/bin/bash

# if this does not run, run this:
# chmod +x commit_submodules.sh

# Ensure all submodules are initialized and updated
git submodule update --init --recursive

# Prompt for submodule commit message
echo "Enter commit message for submodules:"
read submodule_commit_message

# Iterate over each submodule
git submodule foreach ' 
  echo "Processing submodule $name in $path"
  git add .
  git commit -m "$submodule_commit_message"
  git push
'

# Prompt for parent repository commit message
echo "Enter commit message for parent repository:"
read parent_commit_message

# Commit and push changes in the parent repository
echo "Updating parent repository..."
git add .
git commit -m "$parent_commit_message"
git push

echo "All submodules and parent repository processed."