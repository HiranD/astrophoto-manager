#!/usr/bin/env python3
"""
Astrophotography Session Organizer

This script reorganizes filter directories into session directories based on sensor angles.
Each session groups files with similar telescope angles, making it easier to apply
appropriate flat frames for calibration.

Usage:
    python3 organize_sessions.py /path/to/filter/directories

Example:
    python3 organize_sessions.py "/Users/user/Pictures/ASTRO/2025 July/The War and Peace Nebula copy"

The script will:
1. Scan filter subdirectories (L/, R/, G/, B/, Ha/, OIII/, SII/, ALPT/, ALPT2/, etc.)
2. Extract sensor angles from FITS filenames
3. Group files by angle (±2.5° tolerance) into sessions
4. Create session directories (Session_70deg/, Session_71deg/, etc.)
5. Move filter directories into appropriate session directories
6. Show a summary of sessions created

Directory structure before:
    target_dir/
    ├── ALPT/
    │   ├── file_70.5_angle.fits
    │   └── file_71.2_angle.fits
    └── ALPT2/
        ├── file_70.8_angle.fits
        └── file_71.1_angle.fits

Directory structure after:
    target_dir/
    ├── Session_70deg/
    │   ├── ALPT/
    │   │   └── file_70.5_angle.fits
    │   └── ALPT2/
    │       └── file_70.8_angle.fits
    └── Session_71deg/
        ├── ALPT/
        │   └── file_71.2_angle.fits
        └── ALPT2/
            └── file_71.1_angle.fits
"""

import os
import glob
import shutil
import sys
import json
from collections import defaultdict
from typing import Dict, List, Set, Optional, Tuple

# Import angle extraction from the main exposures script
try:
    from exposures import extract_angle, DEFAULT_ORDER, FLAT_ANGLE_TOLERANCE
except ImportError:
    # Fallback if exposures.py is not available
    DEFAULT_ORDER = ["L", "R", "G", "B", "Ha", "OIII", "SII", "ALPT", "ALPT2", "UVIR"]
    FLAT_ANGLE_TOLERANCE = 2.5

    def extract_angle(path: str) -> Optional[float]:
        """Extract telescope angle from FITS filename."""
        basename = os.path.basename(path)
        parts = basename.split('_')

        # Look for angle in format: filter_angle_exposure
        for i, part in enumerate(parts):
            if part in DEFAULT_ORDER and i + 1 < len(parts):
                try:
                    angle_str = parts[i + 1]
                    return float(angle_str)
                except (ValueError, IndexError):
                    continue
        return None


def find_filter_directories(directory: str) -> List[str]:
    """Find all filter directories in the given directory."""
    filter_dirs = []

    if not os.path.exists(directory):
        return filter_dirs

    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)
        if os.path.isdir(item_path) and item in DEFAULT_ORDER:
            filter_dirs.append(item_path)

    return sorted(filter_dirs)


def analyze_filter_angles_for_session(filter_dir: str, session_angle: int, tolerance: float = FLAT_ANGLE_TOLERANCE) -> List[str]:
    """Find files in a filter directory that match a session angle within tolerance."""
    matching_files = []

    fits_files = glob.glob(os.path.join(filter_dir, "*.fits"))

    for fits_file in fits_files:
        angle = extract_angle(fits_file)
        if angle is not None and abs(angle - session_angle) <= tolerance:
            matching_files.append(fits_file)

    return matching_files


def load_available_flat_angles() -> Dict[str, List[int]]:
    """Load available flat frame angles from cache."""
    flat_cache_path = os.path.expanduser("~/.astrophotography_flats.json")

    if not os.path.exists(flat_cache_path):
        # Default flat angles if no cache found
        return {
            "ALPT": [70, 75, 82],
            "ALPT2": [10, 70, 75, 82]
        }

    try:
        with open(flat_cache_path, 'r') as f:
            flat_data = json.load(f)

        # Convert string angles to integers
        flat_angles = {}
        for filter_name, angles_dict in flat_data.items():
            if isinstance(angles_dict, dict):
                flat_angles[filter_name] = [int(angle) for angle in angles_dict.keys()]

        return flat_angles

    except (json.JSONDecodeError, ValueError):
        # Fallback to defaults
        return {
            "ALPT": [70, 75, 82],
            "ALPT2": [10, 70, 75, 82]
        }


def find_all_sessions(directory: str, tolerance: float = FLAT_ANGLE_TOLERANCE) -> Dict[int, Dict[str, bool]]:
    """Find ALL session angles from actual data and check flat availability.

    Returns:
        Dict mapping session_angle to {filter: has_flats}
    """
    available_flats = load_available_flat_angles()
    filter_dirs = find_filter_directories(directory)

    if not filter_dirs:
        return {}

    # Collect all actual angles from image data
    all_image_angles = []

    for filter_dir in filter_dirs:
        fits_files = glob.glob(os.path.join(filter_dir, "*.fits"))
        for fits_file in fits_files:
            angle = extract_angle(fits_file)
            if angle is not None:
                all_image_angles.append(angle)

    if not all_image_angles:
        return {}

    # Get all available flat angles across all filters
    all_flat_angles = set()
    for angles in available_flats.values():
        all_flat_angles.update(angles)

    # New approach: find best flat angle for each image angle range
    all_image_angles = sorted(set(all_image_angles))
    sessions = {}
    assigned_angles = set()

    # Calculate coverage for each flat angle and prioritize by coverage + filter availability
    flat_coverage = []
    for flat_angle in all_flat_angles:
        compatible_angles = [img_angle for img_angle in all_image_angles
                           if abs(img_angle - flat_angle) <= tolerance]

        if compatible_angles:
            # Count how many filters have flats at this angle
            filter_count = sum(1 for filter_dir in filter_dirs
                             if (os.path.basename(filter_dir) in available_flats and
                                 flat_angle in available_flats[os.path.basename(filter_dir)]))

            # Calculate distance from center of compatible angles
            center_angle = sum(compatible_angles) / len(compatible_angles)
            distance_from_center = abs(flat_angle - center_angle)

            # Prioritize by: 1) filter coverage, 2) closeness to center, 3) number of angles, 4) prefer 70° over 71°
            # We want flats that are available for our filters and close to data center
            prefer_70_bonus = 1000 if flat_angle == 70 else 0
            priority = filter_count * 10000 - distance_from_center * 100 + len(compatible_angles) + prefer_70_bonus
            flat_coverage.append((priority, flat_angle, compatible_angles))

    # Sort by priority (descending) to process best coverage first
    flat_coverage.sort(key=lambda x: x[0], reverse=True)


    # Assign angles to flat angles based on priority
    for priority, flat_angle, compatible_angles in flat_coverage:
        # Only assign angles that haven't been assigned yet
        available_angles = [angle for angle in compatible_angles if angle not in assigned_angles]

        if available_angles:
            sessions[flat_angle] = {}
            for filter_dir in filter_dirs:
                filter_name = os.path.basename(filter_dir)
                has_flats = (filter_name in available_flats and
                           flat_angle in available_flats[filter_name])
                sessions[flat_angle][filter_name] = has_flats

            # Mark these angles as assigned
            for angle in available_angles:
                assigned_angles.add(angle)

    # Handle remaining angles that don't match any flat angles
    remaining_angles = [a for a in all_image_angles if a not in assigned_angles]

    if remaining_angles:
        # Group remaining angles by proximity
        remaining_angles.sort()
        clusters = []

        for angle in remaining_angles:
            # Find if this angle belongs to an existing cluster
            placed = False
            for cluster in clusters:
                if any(abs(angle - c) <= tolerance for c in cluster):
                    cluster.append(angle)
                    placed = True
                    break

            if not placed:
                clusters.append([angle])

        # Create sessions for remaining clusters
        for cluster in clusters:
            center_angle = round(sum(cluster) / len(cluster))

            sessions[center_angle] = {}
            for filter_dir in filter_dirs:
                filter_name = os.path.basename(filter_dir)
                has_flats = (filter_name in available_flats and
                           center_angle in available_flats[filter_name])
                sessions[center_angle][filter_name] = has_flats

    return sessions


def organize_into_sessions(directory: str, dry_run: bool = False, tolerance: float = FLAT_ANGLE_TOLERANCE) -> Dict[int, Dict[str, int]]:
    """Organize filter directories into session directories by angle based on available flats.

    Parameters
    ----------
    directory : str
        Path to directory containing filter subdirectories
    dry_run : bool
        If True, show what would be done without actually moving files
    tolerance : float
        Angle tolerance in degrees (default: ±2.5°)

    Returns
    -------
    dict
        Dictionary of {session_angle: {filter: file_count}}
    """
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory not found: {directory}")

    # Load available flat frame angles
    available_flats = load_available_flat_angles()
    print("Available flat frame angles:")
    for filter_name, angles in available_flats.items():
        print(f"  {filter_name}: {angles}°")

    # Find filter directories
    filter_dirs = find_filter_directories(directory)
    if not filter_dirs:
        print(f"No filter directories found in {directory}")
        print("Run organize_files.py first to create filter directories!")
        return {}

    print(f"\nFound {len(filter_dirs)} filter directories: {[os.path.basename(d) for d in filter_dirs]}")

    # Find ALL session angles from actual data and check flat availability
    session_info = find_all_sessions(directory, tolerance)
    if not session_info:
        print("No sessions found in image data")
        return {}

    print(f"\nAll sessions found in data (±{tolerance}°):")

    # Show flat availability summary
    sessions_with_flats = []
    sessions_without_flats = []

    for session_angle, filter_flats in session_info.items():
        has_any_flats = any(filter_flats.values())
        if has_any_flats:
            sessions_with_flats.append(session_angle)
        else:
            sessions_without_flats.append(session_angle)

    # Analyze each session
    sessions_data = {}

    for session_angle in sorted(session_info.keys()):
        sessions_data[session_angle] = {}
        session_range = f"{session_angle-tolerance:.0f}° to {session_angle+tolerance:.0f}°"
        print(f"\nSession_{session_angle}deg ({session_range}):")

        for filter_dir in filter_dirs:
            filter_name = os.path.basename(filter_dir)
            has_flats = session_info[session_angle].get(filter_name, False)

            matching_files = analyze_filter_angles_for_session(filter_dir, session_angle, tolerance)

            if matching_files:
                file_count = len(matching_files)
                sessions_data[session_angle][filter_name] = file_count

                # Show actual angle range and flat status
                angles = [extract_angle(f) for f in matching_files if extract_angle(f) is not None]
                if angles:
                    angle_range = f"{min(angles):.2f}° to {max(angles):.2f}°"
                    flat_status = "✓ Flats available" if has_flats else "⚠ No flats - need to capture"
                    print(f"  {filter_name}: {file_count} files ({angle_range}) {flat_status}")

    if dry_run:
        print(f"\n{'='*50}")
        print("DRY RUN - No directories will be moved")
        print(f"{'='*50}")
        return sessions_data

    # Create session directories and reorganize ALL files
    for session_angle in sorted(session_info.keys()):
        if session_angle not in sessions_data or not sessions_data[session_angle]:
            continue

        session_dir = os.path.join(directory, f"Session_{session_angle}deg")
        os.makedirs(session_dir, exist_ok=True)
        print(f"\nCreated session directory: Session_{session_angle}deg/")

        for filter_dir in filter_dirs:
            filter_name = os.path.basename(filter_dir)

            matching_files = analyze_filter_angles_for_session(filter_dir, session_angle, tolerance)

            if matching_files:
                # Create filter subdirectory in session
                session_filter_dir = os.path.join(session_dir, filter_name)
                os.makedirs(session_filter_dir, exist_ok=True)

                # Move files for this session
                for file_path in matching_files:
                    filename = os.path.basename(file_path)
                    destination = os.path.join(session_filter_dir, filename)

                    if os.path.exists(destination):
                        print(f"    Warning: {filename} already exists in Session_{session_angle}deg/{filter_name}/ - skipping")
                        continue

                    try:
                        shutil.move(file_path, destination)
                        print(f"    Moved {filename} → Session_{session_angle}deg/{filter_name}/")
                    except Exception as e:
                        print(f"    Error moving {filename}: {e}")

    # Clean up empty filter directories
    for filter_dir in filter_dirs:
        try:
            if not os.listdir(filter_dir):  # Directory is empty
                os.rmdir(filter_dir)
                filter_name = os.path.basename(filter_dir)
                print(f"Removed empty directory: {filter_name}/")
        except OSError:
            pass  # Directory not empty or other error

    # Show flat frame recommendations
    show_flat_recommendations(session_info)

    return sessions_data


def show_flat_recommendations(session_info: Dict[int, Dict[str, bool]]):
    """Show recommendations for missing flat frames."""
    print(f"\n{'='*60}")
    print("FLAT FRAME RECOMMENDATIONS")
    print(f"{'='*60}")

    missing_flats = {}

    for session_angle, filter_flats in session_info.items():
        for filter_name, has_flats in filter_flats.items():
            if not has_flats:
                if filter_name not in missing_flats:
                    missing_flats[filter_name] = []
                missing_flats[filter_name].append(session_angle)

    if not missing_flats:
        print("✅ All sessions have flat frames available!")
        return

    print("The following flat frames need to be captured:")
    print()

    for filter_name in sorted(missing_flats.keys()):
        angles = sorted(missing_flats[filter_name])
        print(f"📷 {filter_name} filter:")
        for angle in angles:
            print(f"   • Capture flats at {angle}° telescope angle")
        print()

    print("💡 Tip: Capture flat frames at these angles to enable full calibration")
    print("   of all your image sessions.")



def show_session_summary(sessions_data: Dict[int, Dict[str, int]]):
    """Show summary of sessions created."""
    if not sessions_data:
        print("\nNo sessions were created.")
        return

    print(f"\n{'='*50}")
    print("SESSION ORGANIZATION COMPLETE")
    print(f"{'='*50}")

    total_sessions = len(sessions_data)
    total_files = sum(sum(filters.values()) for filters in sessions_data.values())

    print(f"Created {total_sessions} session directories with {total_files} total files:")

    for session_angle in sorted(sessions_data.keys()):
        filters = sessions_data[session_angle]
        session_file_count = sum(filters.values())
        print(f"\n  Session_{session_angle}deg/ ({session_file_count} files):")

        for filter_name in sorted(filters.keys()):
            file_count = filters[filter_name]
            print(f"    {filter_name}/: {file_count} files")


def main():
    """Main function for command-line usage."""
    if len(sys.argv) < 2:
        print("Usage: python3 organize_sessions.py <directory> [--dry-run]")
        print("\nExample:")
        print('  python3 organize_sessions.py "/Users/user/Pictures/ASTRO/2025 July/The War and Peace Nebula copy"')
        print("\nOptions:")
        print("  --dry-run    Show what would be done without actually moving files")
        sys.exit(1)

    directory = sys.argv[1]
    dry_run = "--dry-run" in sys.argv

    try:
        # Get user confirmation unless it's a dry run
        if not dry_run:
            print(f"This will organize filter directories into sessions based on available flat frames: {directory}")
            print("Only files with matching flat frame angles will be moved into session subdirectories.")
            response = input("\nContinue? (y/N): ").strip().lower()
            if response != 'y':
                print("Cancelled.")
                sys.exit(0)

        sessions_data = organize_into_sessions(directory, dry_run=dry_run, tolerance=FLAT_ANGLE_TOLERANCE)
        show_session_summary(sessions_data)

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