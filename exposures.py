import os
import glob
import json
from collections import defaultdict
from tabulate import tabulate
from typing import List, Tuple, Dict, Any, Optional, Set
from datetime import datetime

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

# Path to store flat frame history
FLAT_HISTORY_FILE: str = os.path.expanduser("~/.astrophotography_flats.json")

# Tolerance for flat frame angle matching (±degrees)
FLAT_ANGLE_TOLERANCE: float = 2.5

# Maximum age for flat frames before warning (days)
FLAT_MAX_AGE_DAYS: int = 90

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


def load_flat_history() -> Dict[str, Dict[float, str]]:
    """Load flat frame history from persistent storage.

    Returns dict: {filter: {angle: date_captured}}
    Filter is the primary key since flats are specific to optical setup, not object.
    Supports both old format (list of angles) and new format (dict with dates).
    Supports both integer and decimal angles.
    """
    if os.path.exists(FLAT_HISTORY_FILE):
        try:
            with open(FLAT_HISTORY_FILE, 'r') as f:
                data = json.load(f)

                # Convert old format to new format with today's date
                result = {}
                today = datetime.now().strftime("%Y-%m-%d")

                for filt, angles in data.items():
                    if isinstance(angles, list):
                        # Old format: list of angles -> convert to dict with today's date
                        result[filt] = {float(angle): today for angle in angles}
                    elif isinstance(angles, dict):
                        # New format: already has dates, convert string keys to float
                        result[filt] = {float(angle): date for angle, date in angles.items()}

                return result
        except (json.JSONDecodeError, KeyError, ValueError):
            return {}
    return {}


def save_flat_history(history: Dict[str, Dict[float, str]]):
    """Save flat frame history to persistent storage."""
    # Convert float keys to strings for JSON serialization
    data = {filt: {str(angle): date for angle, date in angles.items()}
           for filt, angles in history.items()}

    with open(FLAT_HISTORY_FILE, 'w') as f:
        json.dump(data, f, indent=2)


def find_matching_flat(suggested_angle: float, existing_flats: Dict[float, str]) -> Optional[Tuple[float, str]]:
    """Find an existing flat within tolerance of the suggested angle.

    Parameters
    ----------
    suggested_angle : float
        The angle we need a flat for
    existing_flats : dict of {float: str}
        Available flat frame angles and their capture dates

    Returns
    -------
    tuple of (float, str) or None
        The (angle, date) of closest existing flat within tolerance, or None if no match
    """
    for existing_angle, date in existing_flats.items():
        if abs(existing_angle - suggested_angle) <= FLAT_ANGLE_TOLERANCE:
            return (existing_angle, date)
    return None


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


def suggest_flat_angles(angles: List[float], existing_flats: Optional[Dict[float, str]] = None) -> List[int]:
    """Suggest optimal flat frame angles based on light frame angles.
    Now considers existing flats to avoid suggesting redundant angles.

    Parameters
    ----------
    angles : List[float]
        Light frame angles from image data
    existing_flats : Dict[int, str], optional
        Existing flat frame angles and their capture dates

    Returns
    -------
    List[int]
        List of rounded integer angles that need new flat frames

    Strategy:
    1. Find which angles are already covered by existing flats (±2° tolerance)
    2. For uncovered angles, suggest optimal flat positions
    3. Prefer existing flat angles when they can cover the data
    """
    if not angles:
        return []

    if existing_flats is None:
        existing_flats = {}

    unique_angles = sorted(set(angles))
    min_angle = min(unique_angles)
    max_angle = max(unique_angles)
    angle_range = max_angle - min_angle

    # Count frequency of each angle
    angle_counts = {}
    for angle in angles:
        angle_counts[angle] = angle_counts.get(angle, 0) + 1

    # Step 1: Find which angles are already covered by existing flats
    uncovered_angles = []
    covered_by_existing = []

    for angle in unique_angles:
        is_covered = False
        for flat_angle in existing_flats.keys():
            if abs(angle - flat_angle) <= FLAT_ANGLE_TOLERANCE:
                is_covered = True
                if flat_angle not in covered_by_existing:
                    covered_by_existing.append(flat_angle)
                break

        if not is_covered:
            uncovered_angles.append(angle)

    # If all angles are covered by existing flats, return empty list
    if not uncovered_angles:
        return []

    # Step 2: For uncovered angles, suggest optimal positions
    # If we have existing flats, try to extend coverage efficiently
    if existing_flats and uncovered_angles:
        suggested = []

        # Group uncovered angles into clusters
        clusters = []
        for angle in uncovered_angles:
            # Find if this angle belongs to an existing cluster
            placed = False
            for cluster in clusters:
                if any(abs(angle - c) <= FLAT_ANGLE_TOLERANCE * 2 for c in cluster):
                    cluster.append(angle)
                    placed = True
                    break

            if not placed:
                clusters.append([angle])

        # For each cluster, suggest one flat angle
        for cluster in clusters:
            # Use weighted mean of cluster
            cluster_counts = {a: angle_counts[a] for a in cluster}
            weighted_sum = sum(angle * count for angle, count in cluster_counts.items())
            total_count = sum(cluster_counts.values())
            cluster_mean = weighted_sum / total_count
            suggested.append(round(cluster_mean))

        return sorted(list(set(suggested)))

    # Step 3: Fallback to original algorithm for uncovered angles only
    if not uncovered_angles:
        return []

    uncovered_range = max(uncovered_angles) - min(uncovered_angles)
    uncovered_counts = {a: angle_counts[a] for a in uncovered_angles}

    # Small range - single flat angle at weighted mean
    if uncovered_range < 10:
        weighted_sum = sum(angle * count for angle, count in uncovered_counts.items())
        total_count = sum(uncovered_counts.values())
        mean_angle = round(weighted_sum / total_count)
        return [mean_angle]

    # Medium range - suggest 2-3 angles
    elif uncovered_range <= 30:
        if len(uncovered_angles) <= 2:
            return [round(a) for a in uncovered_angles]
        else:
            suggested = []
            # Include min and max of uncovered range
            suggested.append(round(min(uncovered_angles)))
            suggested.append(round(max(uncovered_angles)))

            # Add middle if range is large enough
            if uncovered_range > 15:
                median_angle = uncovered_angles[len(uncovered_angles) // 2]
                suggested.append(round(median_angle))

            return sorted(list(set(suggested)))

    # Large range - suggest key angles based on distribution
    else:
        # Group into bins and find peaks
        bin_size = 10  # 10-degree bins
        bins = {}
        for angle, count in uncovered_counts.items():
            bin_idx = int(angle // bin_size)
            if bin_idx not in bins:
                bins[bin_idx] = []
            bins[bin_idx].extend([angle] * count)

        # Find top bins by frame count
        sorted_bins = sorted(bins.items(), key=lambda x: len(x[1]), reverse=True)

        suggested = []
        for bin_idx, bin_angles in sorted_bins[:3]:  # Top 3 bins for uncovered
            # Use weighted mean of bin
            mean_bin_angle = sum(bin_angles) / len(bin_angles)
            suggested.append(round(mean_bin_angle))

        return sorted(list(set(suggested)))


def update_flat_history(flats: Dict[str, List[float]], replace: bool = False, date_captured: Optional[str] = None):
    """Update flat frame history by filter.

    Parameters
    ----------
    flats : dict
        Dictionary of {filter: [angles]} that have been captured (supports decimal angles)
    replace : bool
        If True, replace existing angles. If False, add to existing.
    date_captured : str, optional
        Date when flats were captured (YYYY-MM-DD). If None, uses today's date.
    """
    history = load_flat_history()

    if date_captured is None:
        date_captured = datetime.now().strftime("%Y-%m-%d")

    for filt, angles in flats.items():
        if replace:
            history[filt] = {float(angle): date_captured for angle in angles}
        else:
            if filt not in history:
                history[filt] = {}
            for angle in angles:
                history[filt][float(angle)] = date_captured

    save_flat_history(history)
    print(f"Updated flat history")
    print(f"Filters updated: {', '.join(flats.keys())}")
    print(f"Date: {date_captured}")


def list_flat_history():
    """Display all flat frame history with dates and age warnings in table format."""
    history = load_flat_history()

    if not history:
        print("No flat frame history found.")
        return

    print("\nFLAT FRAME HISTORY")
    print("="*60)

    # Build table rows
    table_rows = []
    for filt in sorted(history.keys()):
        for angle in sorted(history[filt].keys()):
            date = history[filt][angle]
            age_days = (datetime.now() - datetime.strptime(date, "%Y-%m-%d")).days

            if is_flat_old(date):
                status = f"⚠ OLD ({age_days}d)"
            else:
                status = f"✓ Recent ({age_days}d)"

            table_rows.append([filt, f"{angle}°", date, status])

    if table_rows:
        print(tabulate(
            table_rows,
            headers=["Filter", "Angle", "Date Captured", "Status"],
            tablefmt="grid"
        ))

        # Summary
        total_flats = len(table_rows)
        old_flats = sum(1 for row in table_rows if "OLD" in row[3])

        print(f"\nSummary: {total_flats} flat frames total, {old_flats} need updating (>{FLAT_MAX_AGE_DAYS} days old)")


def clear_flat_history(filters: Optional[List[str]] = None):
    """Clear flat frame history for filters.

    Parameters
    ----------
    filters : list of str, optional
        Specific filters to clear. If None, clears all filters.
    """
    history = load_flat_history()

    if not history:
        print("No flat history found")
        return

    if filters is None:
        # Clear all filters
        history.clear()
        print("Cleared all flat history")
    else:
        # Clear specific filters
        cleared = []
        for filt in filters:
            if filt in history:
                del history[filt]
                cleared.append(filt)

        if cleared:
            print(f"Cleared flat history for filters: {', '.join(cleared)}")
        else:
            print(f"No matching filters found: {', '.join(filters)}")

    save_flat_history(history)


def remove_flat_angles(filter_name: str, angles: List[int]):
    """Remove specific angles from flat frame history.

    Parameters
    ----------
    filter_name : str
        Filter to remove angles from
    angles : list of int
        Specific angles to remove
    """
    history = load_flat_history()

    if filter_name not in history:
        print(f"No flat history found for filter {filter_name}")
        return

    # Remove specified angles
    original = history[filter_name].copy()
    for angle in angles:
        history[filter_name].discard(angle)

    # Remove filter if no angles left
    if not history[filter_name]:
        del history[filter_name]

    removed = original - history.get(filter_name, set())
    if removed:
        save_flat_history(history)
        print(f"Removed angles from {filter_name}: {', '.join(f'{a}°' for a in sorted(removed))}")
        remaining = history.get(filter_name, set())
        if remaining:
            print(f"Remaining angles: {', '.join(f'{a}°' for a in sorted(remaining))}")
        else:
            print(f"No angles remaining for {filter_name}")
    else:
        print(f"No matching angles found to remove")


def reset_all_flat_history(confirm: bool = False):
    """Reset all flat frame history.

    Parameters
    ----------
    confirm : bool
        Must be True to actually reset. Safety mechanism to prevent accidental deletion.
    """
    if not confirm:
        print("WARNING: This will delete ALL flat frame history for ALL objects!")
        print("To confirm, run: reset_all_flat_history(confirm=True)")
        return

    if os.path.exists(FLAT_HISTORY_FILE):
        os.remove(FLAT_HISTORY_FILE)
        print("All flat frame history has been reset.")
    else:
        print("No flat frame history file found.")


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
        print(f"{'='*60}")

        # Load flat history
        flat_history = load_flat_history()

        # Update history with current session if user confirms
        current_flats = {}

        flat_rows = []
        for filt in sorted(flat_data.keys()):
            # Get existing flats for this filter first
            existing = flat_history.get(filt, {})

            # Get detected angles from image data
            detected_angles = sorted(set(flat_data[filt]))

            # Pass existing flats to suggestion algorithm
            suggested = suggest_flat_angles(flat_data[filt], existing)
            current_flats[filt] = set(suggested)

            # Determine which angles are covered by existing flats (within tolerance)
            covered_angles = []
            new_angles = []
            old_flats = []
            recent_flats = []

            for suggested_angle in suggested:
                matching_flat = find_matching_flat(suggested_angle, existing)
                if matching_flat is not None:
                    angle, date = matching_flat
                    covered_angles.append(suggested_angle)
                    if is_flat_old(date):
                        old_flats.append((angle, date))
                    else:
                        recent_flats.append((angle, date))
                else:
                    new_angles.append(suggested_angle)

            # Analyze coverage groups
            coverage_groups = {}
            uncovered_angles = []

            if existing and detected_angles:
                # Group detected angles by which flat covers them
                for det_angle in detected_angles:
                    covered_by = None
                    for flat_angle in existing.keys():
                        if abs(det_angle - flat_angle) <= FLAT_ANGLE_TOLERANCE:
                            covered_by = flat_angle
                            break

                    if covered_by is not None:
                        if covered_by not in coverage_groups:
                            coverage_groups[covered_by] = []
                        coverage_groups[covered_by].append(det_angle)
                    else:
                        uncovered_angles.append(det_angle)

            # Format coverage groups display
            if coverage_groups:
                group_strs = []
                for flat_angle in sorted(coverage_groups.keys()):
                    covered = coverage_groups[flat_angle]
                    if len(covered) == 1:
                        range_str = f"{covered[0]:.1f}°"
                    else:
                        range_str = f"{min(covered):.1f}°-{max(covered):.1f}°"
                    group_strs.append(f"{flat_angle:.1f}°→{range_str}")

                if uncovered_angles:
                    if len(uncovered_angles) == 1:
                        uncov_str = f"{uncovered_angles[0]:.1f}°"
                    else:
                        uncov_str = f"{min(uncovered_angles):.1f}°-{max(uncovered_angles):.1f}°"
                    group_strs.append(f"?→{uncov_str}")

                coverage_str = "; ".join(group_strs)
            elif uncovered_angles:
                if len(uncovered_angles) == 1:
                    coverage_str = f"?→{uncovered_angles[0]:.1f}°"
                else:
                    coverage_str = f"?→{min(uncovered_angles):.1f}°-{max(uncovered_angles):.1f}°"
            else:
                coverage_str = "No data"

            # Format detected angles display with range
            if detected_angles:
                min_angle = min(detected_angles)
                max_angle = max(detected_angles)
                if len(detected_angles) == 1:
                    detected_str = f"{detected_angles[0]:.1f}°"
                elif len(detected_angles) <= 3:
                    detected_str = ", ".join(f"{a:.1f}°" for a in detected_angles)
                else:
                    detected_str = f"{min_angle:.1f}°-{max_angle:.1f}° ({len(detected_angles)} angles)"
            else:
                detected_str = "None"

            if existing:
                # Format existing flats with age indicators
                existing_display = []
                for angle in sorted(existing.keys()):
                    date = existing[angle]
                    if is_flat_old(date):
                        existing_display.append(f"{angle}°⚠")
                    else:
                        existing_display.append(f"{angle}°")

                existing_str = ", ".join(existing_display)

                # Determine status considering age
                if not new_angles and not old_flats:
                    status = "✓ Have"
                elif old_flats and not new_angles:
                    status = "⚠ Old"
                elif new_angles:
                    status = "⚠ Need More"
                else:
                    status = "⚠ Need More"

                # Show what's still needed
                needed = []
                if new_angles:
                    needed.extend([f"{a}°" for a in new_angles])
                if old_flats:
                    needed.extend([f"{a}°(old)" for a, _ in old_flats])

                new_str = ", ".join(needed) if needed else "—"
            else:
                existing_str = "None"
                status = "✗ Need All"
                new_str = ", ".join(f"{a}°" for a in suggested) if suggested else "—"

            flat_rows.append([filt, detected_str, existing_str, new_str, status, coverage_str])

        print(tabulate(
            flat_rows,
            headers=["Filter", "Detected Angles", "Existing Flats", "New Flats Needed", "Status", "Coverage Groups"],
            tablefmt="grid"
        ))

        # Ask user if they want to update the flat history
        print(f"\nFlat frame history stored in: {FLAT_HISTORY_FILE}")

        # Provide option to mark flats as completed
        print("\nTo mark flat frames as captured, run:")
        print(f"  python3 -c \"from exposures import update_flat_history; update_flat_history({dict(current_flats)})\"")


# -----------------------------------------------------------------------------
# Example Usage
# -----------------------------------------------------------------------------
# if __name__ == "__main__":
#     pct = {"L": 40, "R": 20, "G": 20, "B": 20}  # percentages
#     calculate_total_exposure(
#         "/path/to/data",
#         total_target_hours=10.0,  # decimal hours
#         filter_pct=pct,
#         sort_by="remaining_time",
#     )