#!/bin/bash

# SAIGEQTL Documentation Cleanup Script
# This removes development docs from GitHub and keeps only user-facing files

echo "🧹 Cleaning up SAIGEQTL documentation..."

# Files to remove from repository (development/internal docs)
FILES_TO_REMOVE=(
    "SAIGEQTL_DEVELOPMENT_WORKFLOW.md"
    "BUILD_BINARY_GUIDE.md"
    "BINARY_TEST_STEPS.md" 
    "RELEASE_AUTOMATION.md"
    "OPTIMIZATION_README.md"
    "README_TESTING.md"
    "INSTALLATION_OLD.md"
)

# Files to keep for users
FILES_TO_KEEP=(
    "README.md"
    "INSTALLATION.md" 
)

# Files to hide (rename with .)
FILES_TO_HIDE=(
    "CLAUDE.md:.claude.md"
)

echo "📋 Current documentation files:"
ls -la *.md

echo ""
echo "🗑️  Removing development documentation files..."

for file in "${FILES_TO_REMOVE[@]}"; do
    if [[ -f "$file" ]]; then
        echo "  Removing: $file"
        rm "$file"
    else
        echo "  Not found: $file"
    fi
done

echo ""
echo "🔒 Hiding internal files..."

for pair in "${FILES_TO_HIDE[@]}"; do
    old_name="${pair%%:*}"
    new_name="${pair##*:}"
    
    if [[ -f "$old_name" ]]; then
        echo "  Hiding: $old_name → $new_name"
        mv "$old_name" "$new_name"
    else
        echo "  Not found: $old_name"
    fi
done

echo ""
echo "✅ Documentation cleanup complete!"
echo ""
echo "📚 Remaining user documentation:"
for file in "${FILES_TO_KEEP[@]}"; do
    if [[ -f "$file" ]]; then
        echo "  ✅ $file"
    else
        echo "  ❌ $file (missing)"
    fi
done

echo ""
echo "🔐 Private files (not in git):"
echo "  📝 DEVELOPER_GUIDE_PRIVATE.md (your private guide)"
echo "  🤖 .claude.md (hidden Claude instructions)"

echo ""
echo "🎯 Next steps:"
echo "1. Review the changes: git status"
echo "2. Commit cleanup: git add . && git commit -m 'Clean up documentation - keep only user guides'"
echo "3. Add DEVELOPER_GUIDE_PRIVATE.md to .gitignore"