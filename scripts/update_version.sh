#!/bin/bash

# SAIGEQTL Version Update Script
# This script updates the version and triggers automated builds

set -e

if [ $# -eq 0 ]; then
    echo "Usage: $0 <new_version>"
    echo "Example: $0 0.3.5"
    exit 1
fi

NEW_VERSION="$1"

echo "🔄 Updating SAIGEQTL to version $NEW_VERSION"

# Validate version format
if [[ ! $NEW_VERSION =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "❌ Error: Version must be in format X.Y.Z (e.g., 0.3.5)"
    exit 1
fi

# Update DESCRIPTION file
echo "📝 Updating DESCRIPTION file..."
sed -i.bak "s/^Version: .*/Version: $NEW_VERSION/" DESCRIPTION

# Check if the change was made
NEW_VERSION_CHECK=$(grep "^Version:" DESCRIPTION | cut -d' ' -f2)
if [ "$NEW_VERSION_CHECK" != "$NEW_VERSION" ]; then
    echo "❌ Error: Failed to update version in DESCRIPTION"
    mv DESCRIPTION.bak DESCRIPTION  # Restore backup
    exit 1
fi

echo "✅ Version updated in DESCRIPTION: $NEW_VERSION"

# Commit and push changes
echo "📤 Committing changes..."
git add DESCRIPTION
git commit -m "Bump version to $NEW_VERSION

🤖 This will automatically:
- Create GitHub release
- Build Docker images: weizhou0/saigeqtl:$NEW_VERSION and weizhou0/saigeqtl:latest  
- Update all installation methods

📖 See INSTALLATION.md for usage"

echo "🚀 Pushing to GitHub..."
git push origin main

echo "
✅ Version update complete! 

🤖 GitHub Actions will now automatically:
1. 📦 Create release v$NEW_VERSION
2. 🐳 Build Docker image weizhou0/saigeqtl:$NEW_VERSION  
3. 🔄 Update weizhou0/saigeqtl:latest
4. 📋 Create release notes with installation instructions

⏱️  Check progress at: https://github.com/$(git config --get remote.origin.url | sed 's/.*github.com[:/]\(.*\)\.git/\1/')/actions

🎉 Users can immediately install with:
   docker run --rm weizhou0/saigeqtl:$NEW_VERSION step1_fitNULLGLMM_qtl.R --help
"

# Cleanup
rm -f DESCRIPTION.bak