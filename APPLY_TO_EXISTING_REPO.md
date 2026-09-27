# Updating the existing GitHub repository

This bundle is intended to replace/update the files in the existing `GenPept-Curated-2025` repository while preserving Git history.

Recommended procedure:

1. Clone/pull the existing repository.
2. Create a backup branch/tag before replacement.
3. Copy the contents of this bundle into the repository root.
4. Remove obsolete public files that contradict the current manuscript, especially the previous CTD-only Step 05 workflow and old convenience dataset tables.
5. Do **not** delete repository history.
6. Run the validation commands in `docs/REPRODUCIBILITY.md`.
7. Review `git diff` before committing.
8. Commit the update with a message such as `Align public repository with manuscript and frozen dataset v1.1`.

Do not copy the internal source-audit report into the public repository.
