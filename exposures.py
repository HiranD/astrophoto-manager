import os
import glob
import json
from collections import defaultdict
from tabulate import tabulate
from typing import List, Tuple, Dict, Any, Optional
from datetime import datetime

try:
    from astropy.io import fits as astropy_fits
    ASTROPY_AVAILABLE = True
except ImportError:
    ASTROPY_AVAILABLE = False
    print("Warning: astropy not installed. FITS header reading disabled.")

from tqdm import tqdm

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Only these seven canonical filter folders are considered valid
DEFAULT_ORDER = ["L", "R", "G", "B", "Ha", "OIII", "SII", "ALPT", "ALPT2", "UVIR"]

# Extra exposure cushion (percentage) to compensate for light-frame rejection
# during quality selection / calibration.  Example: 15.0 means *capture 15 % more
# data* than strictly needed so the final stack still meets the requested
# exposure-time goal after discarding poor frames.
REJECTION_BUFFER_PCT: float = 15.0

# Flat frame directory paths to scan (camera folders)
# Add your flat directories here - will be scanned for FITS with IMAGETYP='FLAT'
FLAT_SCAN_PATHS: List[str] = [
    # "/path/to/camera1/Flats",
    # "/path/to/camera2/Flats",
]

# Base ASTRO folder for automatic date-based folder discovery
# Will scan subfolders matching "YYYY Month" pattern (e.g., "2025 Sep", "2025 July")
ASTRO_BASE_PATH: str = ""  # e.g., "/Users/you/Pictures/ASTRO"

# Year(s) to scan for date-based folders (e.g., ["2025"] or ["2024", "2025"])
ASTRO_SCAN_YEARS: List[str] = ["2025"]

# Tolerance for flat frame angle matching (±degrees)
FLAT_ANGLE_TOLERANCE: float = 2.5

# Default sensor symmetry for angle normalization (degrees)
# 180 for rectangular sensors (181° ≈ 1°)
# 90 for square sensors (91° ≈ 1°)
# Auto-detected from camera name when possible
DEFAULT_SENSOR_SYMMETRY: int = 180

# Known square sensor cameras (90° symmetry)
# Add camera model identifiers here (matched case-insensitive in camera name)
SQUARE_SENSOR_CAMERAS: List[str] = [
    "533",    # ASI533MM, ASI533MC (square 1:1)
]
# Rectangular: 585, 2600, 571 use default 180° symmetry

# Maximum age for flat frames before warning (days)
FLAT_MAX_AGE_DAYS: int = 90

# Path to cache flat scan results (configurable)
FLAT_CACHE_FILE: str = os.path.expanduser("~/.astro_flats_cache.json")

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

def extract_exposure_time(filename: str) -> float:
    """Return exposure time in **seconds** parsed from a FITS filename."""
    try:
        return float(filename.split("_")[-2].rstrip("s"))
    except (IndexError, ValueError):
        return 0.0


def extract_angle(filename: str) -> Optional[float]:
    """Extract angle from FITS filename (appears after filter name).

    Example: 2025-09-14_03-06-28_ALPT_56.26_300.00s_0114.fits
    Returns 56.26 for the above example.
    """
    try:
        parts = os.path.basename(filename).replace("-", "_").split("_")
        # Find filter position first
        filter_idx = -1
        for i, part in enumerate(parts):
            if part in DEFAULT_ORDER:
                filter_idx = i
                break

        # Angle should be the next part after filter
        if filter_idx >= 0 and filter_idx + 1 < len(parts):
            # Try to parse as float
            angle_str = parts[filter_idx + 1]
            return float(angle_str)
    except (IndexError, ValueError):
        pass
    return None


def minutes_to_hm(total_minutes: float) -> Tuple[int, int]:
    """Convert minutes to (hours, minutes) integer tuple."""
    hours = int(total_minutes // 60)
    minutes = int(round(total_minutes % 60))
    return hours, minutes


def hm_string(total_minutes: float) -> str:
    """Return a nicely formatted "H h M m" string from minutes."""
    h, m = minutes_to_hm(total_minutes)
    return f"{h} h {m} m"


def get_sensor_symmetry(camera_name: Optional[str] = None) -> int:
    """Get sensor symmetry based on camera name.

    Returns 90 for square sensors (e.g., ASI533), 180 for rectangular.
    """
    if camera_name:
        camera_upper = camera_name.upper()
        for model in SQUARE_SENSOR_CAMERAS:
            if model in camera_upper:
                return 90
    return DEFAULT_SENSOR_SYMMETRY


def normalize_angle(angle: float, symmetry: int = None) -> float:
    """Normalize angle based on sensor symmetry.

    For rectangular sensors (180° symmetry): 181° → 1°, 270° → 90°
    For square sensors (90° symmetry): 91° → 1°, 181° → 1°

    This allows flats at 1° to match lights at 181° (rectangular)
    or lights at 91°, 181°, 271° (square).
    """
    if symmetry is None:
        symmetry = DEFAULT_SENSOR_SYMMETRY
    return angle % symmetry


def angles_match(angle1: float, angle2: float, tolerance: float = None, symmetry: int = None) -> bool:
    """Check if two angles match within tolerance, considering sensor symmetry.

    Normalizes both angles before comparison to handle wraparound.
    """
    if tolerance is None:
        tolerance = FLAT_ANGLE_TOLERANCE
    if symmetry is None:
        symmetry = DEFAULT_SENSOR_SYMMETRY

    norm1 = normalize_angle(angle1, symmetry)
    norm2 = normalize_angle(angle2, symmetry)

    # Direct comparison
    if abs(norm1 - norm2) <= tolerance:
        return True

    # Handle wraparound near 0/symmetry boundary
    # e.g., norm1=1°, norm2=179° should match if tolerance allows
    if abs(norm1 - norm2 - symmetry) <= tolerance:
        return True
    if abs(norm2 - norm1 - symmetry) <= tolerance:
        return True

    return False


def detect_filter(path: str) -> Optional[str]:
    """Infer canonical filter name from **directory** or **filename**.

    Search order:
        1. Walk upward from the file's parent directory; the first directory
           whose name matches a canonical filter wins.
        2. If no directory matches, scan the filename tokens (split on `_` or
           `-`).  The first token that exactly matches a canonical filter wins.

    Returns ``None`` if no filter keyword is found.
    """
    # 1) Directory-based detection (nearest wins)
    for part in reversed(os.path.dirname(path).split(os.sep)):
        if part in DEFAULT_ORDER:
            return part

    # 2) Filename-based detection
    basename = os.path.basename(path)
    for token in basename.replace("-", "_").split("_"):
        if token in DEFAULT_ORDER:
            return token
    return None


def extract_date_from_filename(filename: str) -> Optional[str]:
    """Extract date from FITS filename in YYYY-MM-DD format.

    Expected format: 2025-09-14_03-06-28_FILTER_ANGLE_EXPOSURE.fits
    """
    try:
        basename = os.path.basename(filename)
        # First part should be date: 2025-09-14
        date_part = basename.split("_")[0]
        # Validate it looks like a date
        datetime.strptime(date_part, "%Y-%m-%d")
        return date_part
    except (IndexError, ValueError):
        return None


def extract_setup_from_header(fits_file: str) -> Optional[Tuple[str, str]]:
    """Extract camera and telescope from FITS header.

    Returns (camera, telescope) tuple or None if headers not found.
    Uses INSTRUME for camera and TELESCOP for telescope.
    """
    if not ASTROPY_AVAILABLE:
        return None

    try:
        with astropy_fits.open(fits_file) as hdul:
            header = hdul[0].header
            camera = header.get("INSTRUME", "").strip()
            telescope = header.get("TELESCOP", "").strip()

            if camera and telescope:
                return (camera, telescope)
    except Exception:
        pass

    return None


def make_setup_key(camera: str, telescope: str) -> str:
    """Create a setup key from camera and telescope names."""
    return f"{camera}|{telescope}"


def parse_setup_key(setup_key: str) -> Tuple[str, str]:
    """Parse a setup key back into camera and telescope."""
    parts = setup_key.split("|", 1)
    if len(parts) == 2:
        return (parts[0], parts[1])
    return (setup_key, "Unknown")


def discover_astro_folders() -> List[str]:
    """Discover folders in ASTRO_BASE_PATH matching patterns.

    Matches:
    - "YYYY Month" pattern (e.g., "2025 Sep", "2025 July")
    - Folders containing "Mosaic" anywhere in name

    Returns list of folder paths to scan for flats.
    """
    discovered = []

    if not ASTRO_BASE_PATH or not os.path.exists(ASTRO_BASE_PATH):
        return discovered

    # Month names to match
    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
              "January", "February", "March", "April", "May", "June",
              "July", "August", "September", "October", "November", "December"]

    for entry in os.listdir(ASTRO_BASE_PATH):
        entry_path = os.path.join(ASTRO_BASE_PATH, entry)
        if not os.path.isdir(entry_path):
            continue

        # Check if folder contains "Mosaic" (case-insensitive)
        if "mosaic" in entry.lower():
            discovered.append(entry_path)
            continue

        # Check if folder matches "YYYY Month" pattern
        for year in ASTRO_SCAN_YEARS:
            if entry.startswith(year):
                for month in months:
                    if month in entry:
                        discovered.append(entry_path)
                        break
                break

    return discovered


def get_all_scan_paths() -> List[str]:
    """Get all paths to scan for flats (manual + auto-discovered)."""
    all_paths = list(FLAT_SCAN_PATHS)  # Start with manual paths
    all_paths.extend(discover_astro_folders())  # Add discovered paths
    return all_paths


# Type alias for flat info: {"date": str, "path": str}
FlatInfoType = Dict[str, str]
# Type alias for the cache structure: {setup: {filter: {angle: {"date": ..., "path": ...}}}}
FlatCacheType = Dict[str, Dict[str, Dict[float, FlatInfoType]]]


def load_flat_cache() -> Optional[FlatCacheType]:
    """Load cached flat scan results from file.

    Returns dict: {setup: {filter: {angle: {"date": ..., "path": ...}}}} or None if cache doesn't exist.
    """
    if not os.path.exists(FLAT_CACHE_FILE):
        return None

    try:
        with open(FLAT_CACHE_FILE, 'r') as f:
            data = json.load(f)
            # Convert string keys back to floats for angles
            result: FlatCacheType = {}
            for setup, filters in data.items():
                result[setup] = {}
                for filt, angles in filters.items():
                    result[setup][filt] = {}
                    for angle, info in angles.items():
                        # Handle old cache format (string date) vs new format (dict)
                        if isinstance(info, str):
                            result[setup][filt][float(angle)] = {"date": info, "path": "Unknown"}
                        else:
                            result[setup][filt][float(angle)] = info
            return result
    except (json.JSONDecodeError, KeyError, ValueError):
        return None


def save_flat_cache(data: FlatCacheType):
    """Save flat scan results to cache file."""
    # Convert float keys to strings for JSON serialization
    json_data = {}
    for setup, filters in data.items():
        json_data[setup] = {}
        for filt, angles in filters.items():
            json_data[setup][filt] = {str(angle): info for angle, info in angles.items()}

    with open(FLAT_CACHE_FILE, 'w') as f:
        json.dump(json_data, f, indent=2)


def scan_flat_directories(use_cache: bool = True, show_progress: bool = False) -> FlatCacheType:
    """Scan FLAT_SCAN_PATHS for flat frames efficiently.

    Optimization: Since flats are always in separate folders:
    1. For each folder, verify FIRST file is IMAGETYP='FLAT'
    2. Extract camera+telescope setup from header
    3. Extract angle/date from ALL filenames (no more header reads)
    4. Cache results for future use

    Parameters
    ----------
    use_cache : bool
        If True, return cached data if available. If False, force rescan.
    show_progress : bool
        If True, show progress bar during scanning.

    Returns dict: {setup: {filter: {angle: date}}}
    """
    # Check cache first
    if use_cache:
        cached = load_flat_cache()
        if cached is not None:
            return cached

    result: FlatCacheType = {}

    # Get all paths to scan (manual + auto-discovered)
    all_scan_paths = get_all_scan_paths()

    if not all_scan_paths:
        return result

    if not ASTROPY_AVAILABLE:
        print("Warning: astropy required for FITS header scanning.")
        return result

    # Collect all folders to process
    all_folders: Dict[str, List[str]] = {}

    for scan_path in all_scan_paths:
        if not os.path.exists(scan_path):
            continue

        # Find all FITS files recursively
        fits_files = glob.glob(os.path.join(scan_path, "**", "*.fits"), recursive=True)

        if not fits_files:
            continue

        # Group files by parent folder
        for fits_file in fits_files:
            folder = os.path.dirname(fits_file)
            if folder not in all_folders:
                all_folders[folder] = []
            all_folders[folder].append(fits_file)

    # Create iterator with optional progress bar
    folder_items = list(all_folders.items())
    if show_progress and folder_items:
        folder_iterator = tqdm(folder_items, desc="Scanning folders", unit="folder")
    else:
        folder_iterator = folder_items

    # Process each folder
    for folder, folder_files in folder_iterator:
        # Verify FIRST file is a FLAT and extract setup (one header read per folder)
        first_file = folder_files[0]
        setup_key = None
        try:
            with astropy_fits.open(first_file) as hdul:
                header = hdul[0].header
                image_type = header.get("IMAGETYP", "").strip().upper()

                if image_type != "FLAT":
                    # Not a flat folder, skip all files in this folder
                    continue

                # Extract camera and telescope for setup key
                camera = header.get("INSTRUME", "").strip()
                telescope = header.get("TELESCOP", "").strip()

                if camera and telescope:
                    setup_key = make_setup_key(camera, telescope)
                else:
                    # Missing setup info, use folder name as fallback
                    setup_key = f"Unknown|{os.path.basename(folder)}"

        except Exception:
            # Can't read header, skip this folder
            continue

        # Initialize setup in result if needed
        if setup_key not in result:
            result[setup_key] = {}

        # Folder verified as flats - extract from ALL filenames (no more header reads)
        for fits_file in folder_files:
            filt = detect_filter(fits_file)
            angle = extract_angle(fits_file)

            if filt is None or angle is None:
                continue

            # Extract date from filename or use file modification time
            date_captured = extract_date_from_filename(fits_file)
            if date_captured is None:
                mtime = os.path.getmtime(fits_file)
                date_captured = datetime.fromtimestamp(mtime).strftime("%Y-%m-%d")

            # Store in result[setup][filter][angle]
            if filt not in result[setup_key]:
                result[setup_key][filt] = {}

            # Keep the most recent date for each angle
            existing = result[setup_key][filt].get(angle)
            if existing is None or date_captured > existing["date"]:
                result[setup_key][filt][angle] = {"date": date_captured, "path": folder}

    # Save to cache for future use
    if result:
        save_flat_cache(result)

    return result


def refresh_flat_cache(delete_existing: bool = False) -> FlatCacheType:
    """Force rescan of flat directories and update cache.

    Use this when you've added new flat frames and want to update the cache.

    Parameters
    ----------
    delete_existing : bool
        If True, delete existing cache file before scanning (hard refresh).
    """
    # Delete existing cache if requested
    if delete_existing and os.path.exists(FLAT_CACHE_FILE):
        os.remove(FLAT_CACHE_FILE)
        print(f"Deleted existing cache: {FLAT_CACHE_FILE}")

    # Show paths being scanned
    discovered = discover_astro_folders()

    print("Scanning flat directories...")
    if FLAT_SCAN_PATHS:
        print(f"  Manual paths: {len(FLAT_SCAN_PATHS)}")
    if discovered:
        print(f"  Auto-discovered ({ASTRO_SCAN_YEARS}): {len(discovered)} folders")
        for path in discovered[:5]:  # Show first 5
            print(f"    - {os.path.basename(path)}")
        if len(discovered) > 5:
            print(f"    ... and {len(discovered) - 5} more")

    data = scan_flat_directories(use_cache=False, show_progress=True)

    # Count totals per setup
    total_angles = 0
    for setup, filters in data.items():
        setup_angles = sum(len(angles) for angles in filters.values())
        total_angles += setup_angles
        camera, telescope = parse_setup_key(setup)
        print(f"  {camera} + {telescope}: {setup_angles} flat angles")

    print(f"\nCache refreshed: {total_angles} flat angles across {len(data)} setups")
    print(f"Cache saved to: {FLAT_CACHE_FILE}")

    return data


def get_setup_from_lights(fits_files: List[str]) -> Optional[str]:
    """Extract camera|telescope setup from light frames.

    Reads header from first file to determine setup.
    Returns setup key or None if not determinable.
    """
    if not fits_files or not ASTROPY_AVAILABLE:
        return None

    # Try first file
    setup = extract_setup_from_header(fits_files[0])
    if setup:
        return make_setup_key(setup[0], setup[1])

    return None


def calculate_flats_needed(angles: List[float], symmetry: int = None) -> List[int]:
    """Calculate optimal flat angles based on light frame distribution.

    Normalizes angles based on sensor symmetry before processing.
    For rectangular sensors: 181° → 1°, so lights at 181° need flats at ~1°

    Logic:
    - If normalized angle range <= 5°: Return ONE angle (weighted mean)
    - If range > 5°: Divide into 5° bins, compute weighted mean for each

    Parameters
    ----------
    angles : List[float]
        List of angles from light frames (may contain duplicates for weighting)
    symmetry : int, optional
        Sensor symmetry in degrees (90 for square, 180 for rectangular)

    Returns
    -------
    List[int]
        List of rounded optimal angles for flat frames (normalized)
    """
    if not angles:
        return []

    if symmetry is None:
        symmetry = DEFAULT_SENSOR_SYMMETRY

    # Normalize all angles based on sensor symmetry
    normalized_angles = [normalize_angle(a, symmetry) for a in angles]

    unique_angles = sorted(set(normalized_angles))
    min_angle = min(unique_angles)
    max_angle = max(unique_angles)
    angle_range = max_angle - min_angle

    # Count frequency of each normalized angle for weighting
    angle_counts: Dict[float, int] = {}
    for angle in normalized_angles:
        angle_counts[angle] = angle_counts.get(angle, 0) + 1

    if angle_range <= 5.0:
        # Single flat: weighted mean favoring angles with more frames
        weighted_sum = sum(angle * count for angle, count in angle_counts.items())
        total_count = sum(angle_counts.values())
        optimal_angle = round(weighted_sum / total_count)
        return [optimal_angle]
    else:
        # Divide into 5° bins
        suggested = []
        bin_size = 5.0

        # Create bins from min to max angle
        current_bin_start = min_angle
        while current_bin_start <= max_angle:
            bin_end = current_bin_start + bin_size

            # Find angles in this bin
            bin_angles = [a for a in normalized_angles if current_bin_start <= a < bin_end]

            if bin_angles:
                # Compute weighted mean for this bin
                bin_counts = {}
                for a in bin_angles:
                    bin_counts[a] = bin_counts.get(a, 0) + 1

                weighted_sum = sum(a * c for a, c in bin_counts.items())
                total_count = sum(bin_counts.values())
                bin_optimal = round(weighted_sum / total_count)
                suggested.append(bin_optimal)

            current_bin_start = bin_end

        return sorted(list(set(suggested)))


def is_flat_old(date_captured: str) -> bool:
    """Check if a flat frame is older than the maximum age threshold.

    Parameters
    ----------
    date_captured : str
        Date in YYYY-MM-DD format

    Returns
    -------
    bool
        True if flat is older than FLAT_MAX_AGE_DAYS
    """
    try:
        captured = datetime.strptime(date_captured, "%Y-%m-%d")
        age_days = (datetime.now() - captured).days
        return age_days > FLAT_MAX_AGE_DAYS
    except ValueError:
        # If date parsing fails, assume it's old
        return True


def list_existing_flats():
    """Display all existing flat frames from scanned directories, organized by setup."""
    all_flats = scan_flat_directories()

    if not all_flats:
        print("No flat frames found in configured directories.")
        print(f"Configure FLAT_SCAN_PATHS in exposures.py to scan for flats.")
        return

    print("\nEXISTING FLAT FRAMES")
    print("="*60)

    total_all = 0
    old_all = 0

    # Display flats organized by setup
    for setup in sorted(all_flats.keys()):
        camera, telescope = parse_setup_key(setup)
        filters = all_flats[setup]

        print(f"\n{camera} | {telescope}")
        print("-" * 50)

        # Build table rows for this setup
        table_rows = []
        for filt in sorted(filters.keys()):
            for angle in sorted(filters[filt].keys()):
                info = filters[filt][angle]
                date = info["date"]
                path = info.get("path", "Unknown")
                # Format path for display:
                # 1. If under ASTRO_BASE_PATH, show relative path
                # 2. Otherwise, replace home dir with ~
                if ASTRO_BASE_PATH and path.startswith(ASTRO_BASE_PATH):
                    display_path = path[len(ASTRO_BASE_PATH):].lstrip(os.sep)
                else:
                    home = os.path.expanduser("~")
                    if path.startswith(home):
                        display_path = "~" + path[len(home):]
                    else:
                        display_path = path

                age_days = (datetime.now() - datetime.strptime(date, "%Y-%m-%d")).days

                if is_flat_old(date):
                    status = f"⚠ OLD ({age_days}d)"
                    old_all += 1
                else:
                    status = f"✓ Recent ({age_days}d)"

                table_rows.append([filt, f"{angle}°", date, display_path, status])
                total_all += 1

        if table_rows:
            print(tabulate(
                table_rows,
                headers=["Filter", "Angle", "Date", "Location", "Status"],
                tablefmt="grid"
            ))

    # Summary
    print(f"\n{'='*60}")
    print(f"Total: {total_all} flat angles across {len(all_flats)} setups")
    if old_all > 0:
        print(f"Warning: {old_all} flats need updating (>{FLAT_MAX_AGE_DAYS} days old)")


# -----------------------------------------------------------------------------
# Main routine
# -----------------------------------------------------------------------------

def calculate_total_exposure(
    directory: str,
    total_target_hours: Optional[float] = None,
    filter_pct: Optional[Dict[str, float]] = None,
    *,
    sort_by: str = "folder",
    include_unknown: bool = False,
):
    """Recursively compute exposure statistics for *directory*.

    Parameters
    ----------
    directory : str
        Root directory that will be searched recursively for ``*.fits`` files.
    total_target_hours : float, optional
        Desired **overall exposure time** for the imaging target, expressed *as
        decimal hours* (e.g. ``7.5`` for 7 h 30 m).  The code inflates this
        value by ``REJECTION_BUFFER_PCT`` so that more photons are captured
        than the final stack strictly needs:

        .. code:: python

            effective_minutes = total_target_hours * 60 * (1 + REJECTION_BUFFER_PCT/100)

        Combined with ``filter_pct`` every filter ``f`` receives a per-filter
        target of ``effective_minutes * (filter_pct[f] / 100)``.
    filter_pct : Mapping[str, float], optional
        Dictionary whose **keys are canonical filters** (``"L"``, ``"R"`` …)
        and **values are the desired percentage** of the overall exposure time
        to allocate to that filter.  Percentages need not sum to exactly
        100 % – any remainder is ignored.
    sort_by : str, default "folder"
        Column key to sort by.  Recognised values are identical to the earlier
        implementation (case-insensitive): ``folder``, ``exposure_time``,
        ``num_subs``, ``total_exposure``, ``remaining_time``,
        ``remaining_subs``.
    include_unknown : bool, default ``False``
        Whether to include files whose filter cannot be detected.
    """

    # ------------------------------------------------------------------
    # Pass 1 – gather statistics at (filter, exposure) granularity
    # ------------------------------------------------------------------
    stats: Dict[str, Dict[float, Dict[str, Any]]] = defaultdict(
        lambda: defaultdict(lambda: {"total_minutes": 0.0, "num_subs": 0, "angles": []})
    )
    grand_total: Dict[float, Dict[str, float]] = defaultdict(
        lambda: {"total_minutes": 0.0, "num_subs": 0}
    )

    fits_files = glob.glob(os.path.join(directory, "**", "*.fits"), recursive=True)

    for file in fits_files:
        exposure_sec = extract_exposure_time(file)
        if exposure_sec <= 0:
            continue

        filt = detect_filter(file)
        if filt is None:
            if not include_unknown:
                continue
            filt = "Unknown"

        minutes = exposure_sec / 60.0
        angle = extract_angle(file)

        stats[filt][exposure_sec]["total_minutes"] += minutes
        stats[filt][exposure_sec]["num_subs"] += 1
        if angle is not None:
            stats[filt][exposure_sec]["angles"].append(angle)

        grand_total[exposure_sec]["total_minutes"] += minutes
        grand_total[exposure_sec]["num_subs"] += 1

    if not stats:
        print(f"No matching FITS files found under '{directory}'.")
        return

    # ------------------------------------------------------------------
    # Compute per-filter aggregates and target minutes (if requested)
    # ------------------------------------------------------------------
    current_minutes_per_filter: Dict[str, float] = {
        f: sum(d["total_minutes"] for d in stats_f.values()) for f, stats_f in stats.items()
    }

    target_minutes_per_filter: Dict[str, float] = {}
    if total_target_hours is not None and filter_pct is not None:
        effective_total_minutes = total_target_hours * 60.0 * (1 + REJECTION_BUFFER_PCT / 100.0)

        # Warn if supplied percentages sum to > 100 (not necessarily wrong but suspicious)
        pct_sum = sum(filter_pct.values())
        if pct_sum > 100.0:
            print(f"Warning: filter_pct values sum to {pct_sum:.1f} %, exceeding 100 %.")

        for f in stats:
            pct = filter_pct.get(f, 0.0)
            target_minutes_per_filter[f] = effective_total_minutes * (pct / 100.0)

    # ------------------------------------------------------------------
    # Build numeric data rows
    # ------------------------------------------------------------------
    data_rows: List[List[Any]] = []
    overall_minutes = 0.0

    ordered_filters = [f for f in DEFAULT_ORDER if f in stats] + [
        f for f in sorted(stats) if f not in DEFAULT_ORDER
    ]

    # Collect all angle data for flat frame analysis
    flat_data: Dict[str, List[float]] = {}

    for filt in ordered_filters:
        for exposure_sec in sorted(stats[filt]):
            d = stats[filt][exposure_sec]
            total_min = d["total_minutes"]
            if d["angles"]:
                if filt not in flat_data:
                    flat_data[filt] = []
                flat_data[filt].extend(d["angles"])
            row = [filt, exposure_sec, d["num_subs"], total_min]

            # Target time and remaining time / subs (buffer-aware)
            if target_minutes_per_filter:
                target = target_minutes_per_filter.get(filt)
                if target is not None:
                    remaining = max(0.0, target - current_minutes_per_filter[filt])
                    additional = (remaining * 60.0 / exposure_sec) if exposure_sec else 0.0
                    row.extend([target, remaining, additional])
                else:
                    row.extend([None, None, None])
            data_rows.append(row)
            overall_minutes += total_min

    # ------------------------------------------------------------------
    # Headers & sort-index resolution
    # ------------------------------------------------------------------
    headers = ["Filter", "Exposure (s)", "# Subs", "Total Exposure"]
    if target_minutes_per_filter:
        headers.extend(["Target Time", "Remaining Time", "Remaining Subs"])

    header_lut = {h.lower(): i for i, h in enumerate(headers)}

    legacy_map = {
        "folder": 0,
        "exposure_time": 1,
        "num_subs": 2,
        "total_exposure": 3,
        "target_time": 4 if len(headers) > 4 else None,
        "remaining_time": 5 if len(headers) > 5 else None,
        "remaining_subs": 6 if len(headers) > 6 else None,
    }

    sort_by_lc = sort_by.lower()
    sort_idx = legacy_map.get(sort_by_lc, header_lut.get(sort_by_lc))

    # ------------------------------------------------------------------
    # Sorting – DESCENDING for numeric columns; canonical for filters
    # ------------------------------------------------------------------
    def _filter_key(r):
        return (
            DEFAULT_ORDER.index(r[0]) if r[0] in DEFAULT_ORDER else len(DEFAULT_ORDER),
            -r[1],  # secondary: exposure seconds descending
        )

    if sort_idx == 0:
        data_rows.sort(key=_filter_key)
    elif sort_idx is not None:
        data_rows.sort(key=lambda r: (r[sort_idx] is None, -r[sort_idx] if r[sort_idx] is not None else 0))
    else:
        print(f"Warning: unsupported sort_by '{sort_by}' – using default order.")
        data_rows.sort(key=_filter_key)

    # ------------------------------------------------------------------
    # Convert to display rows
    # ------------------------------------------------------------------
    display_rows: List[List[Any]] = []
    for r in data_rows:
        display = [r[0], r[1], r[2], hm_string(r[3])]
        if len(r) > 4:
            target_min = r[4]
            rem_min = r[5]
            display.extend([
                hm_string(target_min) if target_min is not None else "—",
                hm_string(rem_min) if rem_min is not None else "—",
                int(round(r[6])) if r[6] is not None else "—",
            ])
        display_rows.append(display)

    print(tabulate(display_rows, headers=headers, tablefmt="grid"))

    # ------------------------------------------------------------------
    # Per-Filter Summary table (only when targets are specified)
    # ------------------------------------------------------------------
    if target_minutes_per_filter:
        summary_rows = []
        for filt in ordered_filters:
            captured = current_minutes_per_filter.get(filt, 0.0)
            target = target_minutes_per_filter.get(filt, 0.0)
            remaining = max(0.0, target - captured)

            if target > 0:
                progress = (captured / target) * 100.0
                progress_str = f"{progress:.1f}%"
            else:
                progress_str = "—"

            summary_rows.append([
                filt,
                hm_string(captured),
                hm_string(target) if target > 0 else "—",
                hm_string(remaining) if target > 0 else "—",
                progress_str
            ])

        print("\nPer-Filter Summary:")
        print(tabulate(
            summary_rows,
            headers=["Filter", "Captured", "Target", "Remaining", "Progress"],
            tablefmt="grid"
        ))

    # ------------------------------------------------------------------
    # Grand-total table
    # ------------------------------------------------------------------
    grand_display: List[List[Any]] = []
    for exposure_sec in sorted(grand_total, reverse=True):  # largest sec first
        g = grand_total[exposure_sec]
        grand_display.append([
            exposure_sec,
            g["num_subs"],
            hm_string(g["total_minutes"]),
        ])

    print("\nGrand Total Exposure for All Filters:")
    print(
        tabulate(
            grand_display,
            headers=["Exposure Time (seconds)", "Total Number of Subs", "Total Exposure"],
            tablefmt="grid",
        )
    )

    oh, om = minutes_to_hm(overall_minutes)
    print(f"\nOverall total exposure (all filters) → {oh} h {om} m")

    # ------------------------------------------------------------------
    # Flat Frame Recommendations Table
    # ------------------------------------------------------------------
    if flat_data:
        print(f"\n{'='*60}")
        print(f"FLAT FRAME RECOMMENDATIONS")

        # Detect setup from light frames
        light_setup = get_setup_from_lights(fits_files)
        camera = None
        if light_setup:
            camera, telescope = parse_setup_key(light_setup)
            print(f"Setup: {camera} | {telescope}")
        else:
            print("Setup: Unknown (could not read FITS headers)")

        # Get sensor symmetry based on camera
        symmetry = get_sensor_symmetry(camera)
        symmetry_type = "square" if symmetry == 90 else "rectangular"
        print(f"Sensor: {symmetry_type} ({symmetry}° symmetry)")

        print(f"{'='*60}")

        # Scan flat directories for existing flats
        all_flats = scan_flat_directories()

        # Get flats only for matching setup
        if light_setup and light_setup in all_flats:
            flat_history = all_flats[light_setup]
        else:
            flat_history = {}
            if light_setup and all_flats:
                print(f"Note: No flats found for this setup. Available setups:")
                for setup in all_flats.keys():
                    c, t = parse_setup_key(setup)
                    print(f"  - {c} | {t}")

        flat_rows = []
        for filt in sorted(flat_data.keys()):
            # Get existing flats for this filter (from matching setup only)
            existing = flat_history.get(filt, {})

            # Get detected angles from image data
            detected_angles = sorted(set(flat_data[filt]))

            # Calculate optimal flat angles needed based on distribution
            optimal_angles = calculate_flats_needed(flat_data[filt], symmetry)

            # Determine which optimal angles are already covered by existing flats
            needed_angles = []
            old_angles = []

            for opt_angle in optimal_angles:
                # Check if any existing flat covers this angle (within tolerance)
                # Uses angles_match() which normalizes for sensor symmetry
                is_covered = False
                covering_flat = None
                for existing_angle, info in existing.items():
                    if angles_match(opt_angle, existing_angle, symmetry=symmetry):
                        is_covered = True
                        covering_flat = (existing_angle, info["date"])
                        break

                if is_covered and covering_flat:
                    # Check if the covering flat is old
                    if is_flat_old(covering_flat[1]):
                        old_angles.append(opt_angle)
                else:
                    needed_angles.append(opt_angle)

            # Format detected angles display (show normalized range)
            if detected_angles:
                # Normalize angles for display using camera-specific symmetry
                norm_angles = sorted(set(normalize_angle(a, symmetry) for a in detected_angles))
                min_norm = min(norm_angles)
                max_norm = max(norm_angles)

                if len(norm_angles) == 1:
                    detected_str = f"{norm_angles[0]:.1f}°"
                elif len(norm_angles) <= 3:
                    detected_str = ", ".join(f"{a:.1f}°" for a in norm_angles)
                else:
                    detected_str = f"{min_norm:.1f}°-{max_norm:.1f}° ({len(norm_angles)} angles)"
            else:
                detected_str = "None"

            # Format matching flats display (only flats that cover this target's angles)
            matching_flats = []
            if existing and optimal_angles:
                for opt_angle in optimal_angles:
                    for flat_angle, info in existing.items():
                        if angles_match(opt_angle, flat_angle, symmetry=symmetry):
                            is_old = is_flat_old(info["date"])
                            matching_flats.append((flat_angle, is_old))
                            break  # Found a match for this optimal angle

            if matching_flats:
                # Remove duplicates and sort
                unique_matches = sorted(set(matching_flats), key=lambda x: x[0])
                matching_display = []
                for angle, is_old in unique_matches:
                    if is_old:
                        matching_display.append(f"{angle}°⚠")
                    else:
                        matching_display.append(f"{angle}°")
                matching_str = ", ".join(matching_display)
            else:
                matching_str = "None"

            # Format flats needed display
            all_needed = []
            if needed_angles:
                all_needed.extend([f"{a}°" for a in needed_angles])
            if old_angles:
                all_needed.extend([f"{a}°(old)" for a in old_angles])

            flats_needed_str = ", ".join(all_needed) if all_needed else "—"

            # Determine status
            if not needed_angles and not old_angles:
                status = "✓ Have"
            elif old_angles and not needed_angles:
                status = "⚠ Old"
            elif not matching_flats:
                status = "✗ Need All"
            else:
                status = "⚠ Need More"

            flat_rows.append([filt, detected_str, matching_str, flats_needed_str, status])

        print(tabulate(
            flat_rows,
            headers=["Filter", "Detected Angles", "Matching Flats", "Flats Needed", "Status"],
            tablefmt="grid"
        ))

        # Show cache info
        if FLAT_SCAN_PATHS:
            if os.path.exists(FLAT_CACHE_FILE):
                cache_mtime = datetime.fromtimestamp(os.path.getmtime(FLAT_CACHE_FILE))
                print(f"\nUsing cached flat data from: {FLAT_CACHE_FILE}")
                print(f"Cache last updated: {cache_mtime.strftime('%Y-%m-%d %H:%M')}")
                print("To refresh: refresh_flat_cache()")
            else:
                print(f"\nFlat data scanned fresh (no cache existed)")
                print(f"Cache saved to: {FLAT_CACHE_FILE}")
        else:
            print("\nNote: Configure FLAT_SCAN_PATHS in exposures.py to auto-detect existing flats.")


# -----------------------------------------------------------------------------
# HOW TO USE THIS SCRIPT
# -----------------------------------------------------------------------------
#
# STEP 1: Configure paths and sensor (edit the variables above)
# --------------------------------------------------------------
#   FLAT_SCAN_PATHS = ["/path/to/camera1/Flats", "/path/to/camera2/Flats"]
#   ASTRO_BASE_PATH = "/Users/you/Pictures/ASTRO"
#   ASTRO_SCAN_YEARS = ["2025"]  # or ["2024", "2025"]
#   Sensor symmetry is auto-detected from camera name:
#     - Square (90°): ASI533 cameras
#     - Rectangular (180°): ASI585, ASI2600, etc.
#   Angle normalization: 181° matches flats at 1° for rectangular sensors
#
#   Auto-discovers folders matching:
#     - "YYYY Month" (e.g., "2025 Sep", "2025 July")
#     - Contains "Mosaic" anywhere in name
#
# STEP 2: Scan for existing flats (first time or after adding new flats)
# -----------------------------------------------------------------------
#   from exposures import refresh_flat_cache
#   refresh_flat_cache()                      # normal refresh
#   refresh_flat_cache(delete_existing=True)  # hard refresh (deletes cache first)
#
# STEP 3: View all existing flats (organized by camera+telescope setup)
# ----------------------------------------------------------------------
#   from exposures import list_existing_flats
#   list_existing_flats()
#
#   Shows: Filter, Angle, Date, Location (folder path), Status
#
# STEP 4: Analyze light frames and get flat recommendations
# ----------------------------------------------------------
#   from exposures import calculate_total_exposure
#
#   calculate_total_exposure(
#       "/path/to/your/light/frames",
#       total_target_hours=10.0,                    # optional: target hours
#       filter_pct={"L": 40, "Ha": 30, "OIII": 30}, # optional: filter %
#   )
#
#   Output shows:
#     - Exposure stats per filter
#     - Remaining time/subs needed
#     - Flat recommendations (setup-aware: camera + telescope)
#     - Which flat angles you need to capture
#
# REQUIREMENTS
# ------------
#   pip install astropy tabulate tqdm
#
# -----------------------------------------------------------------------------