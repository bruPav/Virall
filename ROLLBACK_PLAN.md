# Rollback Plan for Single-Cell Feature Development

## Current Status
- **Backup Branch**: `backup-before-single-cell` (committed and pushed)
- **Working Branch**: `feature/single-cell-support`
- **File Backup**: `/home/Chema/Program/vas_backup_20251018_102617/`

## Rollback Options

### Option 1: Git Rollback (Recommended)
```bash
# If you need to completely revert to the backup state:
git checkout backup-before-single-cell
git checkout -b main-backup
git push origin main-backup

# Or merge backup into main if needed:
git checkout main
git merge backup-before-single-cell
git push origin main
```

### Option 2: File System Rollback
```bash
# If git rollback doesn't work:
cd /home/Chema/Program
rm -rf vas
mv vas_backup_20251018_102617 vas
```

### Option 3: Selective Rollback
```bash
# If only some files are problematic:
git checkout backup-before-single-cell -- path/to/problematic/file.py
```

## Testing Before Rollback

### 1. Run Test Suite
```bash
cd /home/Chema/Program/vas
python tests/test_environment.py
```

### 2. Test Basic Functionality
```bash
# Test CLI help
python -m virall.cli --help

# Test basic assembly (if dependencies available)
python -m virall.cli assemble --help
```

### 3. Check Git Status
```bash
git status
git log --oneline -10
```

## Emergency Contacts
- **Backup Location**: `/home/Chema/Program/vas_backup_20251018_102617/`
- **Git Remote**: `origin/backup-before-single-cell`
- **Test Script**: `tests/test_environment.py`

## Development Strategy
1. **Incremental Changes**: Make small changes and test frequently
2. **Frequent Commits**: Commit after each working feature
3. **Test After Each Change**: Run test suite after modifications
4. **Document Changes**: Keep track of what was modified

## Safe Development Practices
- Always work in the feature branch
- Never modify the backup branch
- Test before committing
- Keep the backup branch as a safety net
