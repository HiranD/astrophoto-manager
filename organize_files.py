#!/usr/bin/env python3
"""
Astrophotography File Organizer

This script organizes FITS files in a directory by moving them into filter-specific subdirectories.
It detects the filter from filenames and creates appropriate folder structure.

Usage:
    python3 organize_files.py /path/to/mixed/files

Example:
    python3 organize_files.py "/Users/user/Pictures/ASTRO/2025 July/SMC Panel 6 copy"

The script will:
1. Scan for FITS files in the specified directory
2. Detect filter names from filenames
3. Create filter subdirectories (L, R, G, B, Ha, OIII, SII, ALPT, etc.)
4. Move files to appropriate subdirectories
5. Show a summary of what was moved
"""

import os
import glob
import shutil
import sys
from collections import defaultdict
from typing import Dict, List, Optional

# Import filter detection from the main exposures script
try:
    from exposures import detect_filter, DEFAULT_ORDER
except ImportError:
    # Fallback if exposures.py is not available
    DEFAULT_ORDER = ["L", "R", "G", "B", "Ha", "OIII", "SII", "ALPT", "ALPT2", "UVIR"]

    def detect_filter(path: str) -> Optional[str]:
        """Detect filter from filename or directory path."""
        # Directory-based detection
        for part in reversed(os.path.dirname(path).split(os.sep)):
            if part in DEFAULT_ORDER:
                return part

        # Filename-based detection
        basename = os.path.basename(path)
        for token in basename.replace("-", "_").split("_"):
            if token in DEFAULT_ORDER:
                return token
        return None


def organize_fits_files(directory: str, dry_run: bool = False) -> Dict[str, List[str]]:
    """Organize FITS files in a directory by filter.

    Parameters
    ----------
    directory : str
        Path to directory containing mixed FITS files
    dry_run : bool
        If True, show what would be done without actually moving files

    Returns
    -------
    dict
        Dictionary of {filter: [list of files moved]}
    """
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory not found: {directory}")

    # Find all FITS files in the directory (non-recursive to avoid moving already organized files)
    fits_pattern = os.path.join(directory, "*.fits")
    fits_files = glob.glob(fits_pattern)

    if not fits_files:
        print(f"No FITS files found in {directory}")
        return {}

    print(f"Found {len(fits_files)} FITS files to organize")

    # Group files by detected filter
    files_by_filter = defaultdict(list)
    unknown_files = []

    for file_path in fits_files:
        filter_name = detect_filter(file_path)
        if filter_name:
            files_by_filter[filter_name].append(file_path)
        else:
            unknown_files.append(file_path)

    # Show summary of what was detected
    print(f"\nDetected filters:")
    for filter_name in sorted(files_by_filter.keys()):
        count = len(files_by_filter[filter_name])
        print(f"  {filter_name}: {count} files")

    if unknown_files:
        print(f"  Unknown: {len(unknown_files)} files")

    if dry_run:
        print(f"\n{'='*50}")
        print("DRY RUN - No files will be moved")
        print(f"{'='*50}")
        return dict(files_by_filter)

    # Create filter directories and move files
    moved_files = {}

    for filter_name, file_list in files_by_filter.items():
        # Create filter subdirectory
        filter_dir = os.path.join(directory, filter_name)
        os.makedirs(filter_dir, exist_ok=True)

        moved_files[filter_name] = []

        for file_path in file_list:
            filename = os.path.basename(file_path)
            destination = os.path.join(filter_dir, filename)

            # Check if destination already exists
            if os.path.exists(destination):
                print(f"  Warning: {filename} already exists in {filter_name}/ - skipping")
                continue

            try:
                shutil.move(file_path, destination)
                moved_files[filter_name].append(filename)
                print(f"  Moved {filename} → {filter_name}/")
            except Exception as e:
                print(f"  Error moving {filename}: {e}")

    # Handle unknown files
    if unknown_files:
        print(f"\nUnknown files (not moved):")
        for file_path in unknown_files:
            print(f"  {os.path.basename(file_path)}")

    return moved_files


def show_summary(moved_files: Dict[str, List[str]]):
    """Show summary of files moved."""
    if not moved_files:
        print("\nNo files were moved.")
        return

    print(f"\n{'='*50}")
    print("ORGANIZATION COMPLETE")
    print(f"{'='*50}")

    total_moved = sum(len(files) for files in moved_files.values())
    print(f"Successfully moved {total_moved} files into {len(moved_files)} filter directories:")

    for filter_name in sorted(moved_files.keys()):
        count = len(moved_files[filter_name])
        print(f"  {filter_name}/: {count} files")


def main():
    """Main function for command-line usage."""
    if len(sys.argv) < 2:
        print("Usage: python3 organize_files.py <directory> [--dry-run]")
        print("\nExample:")
        print('  python3 organize_files.py "/Users/user/Pictures/ASTRO/2025 July/SMC Panel 6 copy"')
        print("\nOptions:")
        print("  --dry-run    Show what would be done without actually moving files")
        sys.exit(1)

    directory = sys.argv[1]
    dry_run = "--dry-run" in sys.argv

    try:
        # Get user confirmation unless it's a dry run
        if not dry_run:
            print(f"This will organize FITS files in: {directory}")
            print("Files will be moved into filter subdirectories (L/, R/, G/, B/, etc.)")
            response = input("\nContinue? (y/N): ").strip().lower()
            if response != 'y':
                print("Cancelled.")
                sys.exit(0)

        moved_files = organize_fits_files(directory, dry_run=dry_run)
        show_summary(moved_files)

    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nCancelled by user.")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()