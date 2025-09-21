# Astrophotography Exposure Tracker & File Organizer

A comprehensive toolkit for astrophotographers featuring exposure tracking, flat frame management, and file organization. This suite helps monitor imaging sessions, calculate total integration time, plan future captures, manage flat frame calibration files, and organize FITS files by filter.

## Features

### Exposure Analysis (`exposures.py`)
- **Automatic Filter Detection**: Intelligently identifies filters (L, R, G, B, Ha, OIII, SII, etc.) from directory structure or filenames
- **Exposure Statistics**: Calculates total exposure time per filter and exposure duration
- **Target Planning**: Set exposure goals and see remaining time needed per filter
- **Rejection Buffer**: Accounts for frame rejection during calibration/stacking
- **Flexible Sorting**: Sort results by various criteria (filter, exposure time, number of subs, etc.)
- **Multi-Exposure Support**: Handles mixed exposure times within the same filter
- **Flat Frame Management**:
  - Analyzes telescope angles from light frames
  - Suggests optimal flat frame angles (rounded to whole degrees)
  - Tracks flat frame history with date tracking
  - Shows which flats you have vs. need
  - Warns about old flats (>90 days)
  - ±2.5° tolerance for angle matching

### File Organization (`organize_files.py`)
- **Automatic Filter Detection**: Uses same logic as exposure analyzer
- **Smart Organization**: Moves FITS files into filter-specific subdirectories
- **Safety Features**: Dry-run mode, confirmation prompts, duplicate detection
- **Progress Tracking**: Shows each file being moved with clear feedback
- **Multi-Filter Support**: Handles all standard filters (L, R, G, B, Ha, OIII, SII, ALPT, etc.)

### Session Organization (`organize_sessions.py`)
- **Angle-Based Sessions**: Groups filter directories into session directories based on telescope angles
- **Flat Frame Integration**: Uses available flat frame angles to determine optimal session grouping
- **Smart Clustering**: Prioritizes flat angles that cover the most data with best filter availability
- **±2° Tolerance**: Groups files within ±2° of flat frame angles for practical calibration
- **Complete Workflow**: Works with organize_files.py for full file organization pipeline
- **Flat Recommendations**: Shows which additional flat frames need to be captured

## Installation

### Prerequisites

- Python 3.6 or higher
- `tabulate` library for formatted table output

### Setup

```bash
pip install tabulate
```

## Usage

### Basic Usage

```python
from exposures import calculate_total_exposure

# Analyze all FITS files in a directory
calculate_total_exposure("/path/to/your/data")
```

### Advanced Usage with Target Planning

```python
from exposures import calculate_total_exposure

# Define filter allocation percentages
filter_allocation = {
    "L": 40,    # 40% to Luminance
    "R": 20,    # 20% to Red
    "G": 20,    # 20% to Green
    "B": 20     # 20% to Blue
}

# Set 10 hours total target and calculate remaining exposures needed
calculate_total_exposure(
    "/path/to/your/data",
    total_target_hours=10.0,  # Target 10 hours total
    filter_pct=filter_allocation,
    sort_by="remaining_time"   # Sort by filters needing most time
)
```

## API Reference

### Main Functions

#### `calculate_total_exposure(directory, total_target_hours=None, filter_pct=None, sort_by="folder", include_unknown=False)`

Analyzes FITS files in a directory and generates exposure statistics with flat frame recommendations.

**Parameters:**
- **`directory`** (str): Root directory to recursively search for FITS files
- **`total_target_hours`** (float, optional): Total desired exposure time in decimal hours (e.g., 7.5 = 7h 30m)
- **`filter_pct`** (dict, optional): Dictionary mapping filter names to percentage allocations
- **`sort_by`** (str, default="folder"): Sort criteria - options include:
  - `"folder"`: Sort by filter in canonical order
  - `"exposure_time"`: Sort by exposure duration
  - `"num_subs"`: Sort by number of subframes
  - `"total_exposure"`: Sort by total integration time
  - `"remaining_time"`: Sort by time still needed (requires targets)
  - `"remaining_subs"`: Sort by subframes still needed
- **`include_unknown`** (bool, default=False): Include files where filter cannot be detected

**Example:**
```python
from exposures import calculate_total_exposure

# Basic analysis
calculate_total_exposure("/data/M31")

# With target planning
calculate_total_exposure(
    "/data/M31",
    total_target_hours=10.0,
    filter_pct={"L": 40, "R": 20, "G": 20, "B": 20}
)
```

### Flat Frame Management Functions

#### `update_flat_history(flats, replace=False, date_captured=None)`

Records which flat frames have been captured, organized by filter with date tracking.

**Parameters:**
- **`flats`** (dict): Dictionary of {filter: [angles]} that have been captured
- **`replace`** (bool, default=False): If True, replaces existing angles. If False, adds to existing.
- **`date_captured`** (str, optional): Date in YYYY-MM-DD format. Defaults to today if not specified.

**Examples:**
```python
from exposures import update_flat_history

# Add new flats by filter (uses today's date)
update_flat_history({'L': [45, 56], 'R': [33]})

# Add more L flats (keeps existing ones)
update_flat_history({'L': [67, 89]})

# Add flats with specific date
update_flat_history({'L': [45, 56]}, date_captured='2025-09-15')

# Replace all L flats (overwrites existing)
update_flat_history({'L': [45, 56, 67, 89]}, replace=True)

# Clear all flats for a filter
update_flat_history({'L': []}, replace=True)
```

#### `list_flat_history()`

Displays all flat frame history in a table format with dates and aging status.

**Example:**
```python
from exposures import list_flat_history

# Display all flat history
list_flat_history()

# Output:
# FLAT FRAME HISTORY
# ============================================================
# +----------+---------+-----------------+---------------+
# | Filter   | Angle   | Date Captured   | Status        |
# +==========+=========+=================+===============+
# | ALPT     | 70°     | 2025-09-21      | ✓ Recent (0d) |
# | ALPT     | 75°     | 2025-09-21      | ✓ Recent (0d) |
# | L        | 82°     | 2025-09-21      | ✓ Recent (0d) |
# | Ha       | 45°     | 2025-04-01      | ⚠ OLD (173d)  |
# +----------+---------+-----------------+---------------+
#
# Summary: 4 flat frames total, 1 need updating (>90 days old)
```

#### `clear_flat_history(filters=None)`

Clears flat frame history for specific filters.

**Parameters:**
- **`filters`** (list of str, optional): Specific filters to clear. If None, clears all filters.

**Examples:**
```python
from exposures import clear_flat_history

# Clear all flat history
clear_flat_history()

# Clear only L and R filters (keep others)
clear_flat_history(['L', 'R'])

# Clear just Ha filter
clear_flat_history(['Ha'])
```

#### `remove_flat_angles(filter_name, angles)`

Removes specific angles from flat frame history while keeping others.

**Parameters:**
- **`filter_name`** (str): Filter to remove angles from
- **`angles`** (list of int): Specific angles to remove

**Examples:**
```python
from exposures import remove_flat_angles

# Remove just 45° from L filter, keep others
remove_flat_angles('L', [45])

# Remove multiple angles from Ha filter
remove_flat_angles('Ha', [30, 60, 90])

# Output shows what was removed and what remains:
# Removed angles from L: 45°
# Remaining angles: 56°, 67°, 89°
```

#### `reset_all_flat_history(confirm=False)`

Resets all flat frame history for all filters. Requires confirmation to prevent accidental deletion.

**Parameters:**
- **`confirm`** (bool): Must be True to actually reset. Safety mechanism.

**Examples:**
```python
from exposures import reset_all_flat_history

# Safe check - shows warning without deleting
reset_all_flat_history()
# Output: WARNING: This will delete ALL flat frame history for ALL objects!
#         To confirm, run: reset_all_flat_history(confirm=True)

# Actually reset everything
reset_all_flat_history(confirm=True)
# Output: All flat frame history has been reset.
```

### Helper Functions

#### `extract_exposure_time(filename)`

Extracts exposure time in seconds from a FITS filename.

**Parameters:**
- **`filename`** (str): FITS filename

**Returns:**
- **float**: Exposure time in seconds (0.0 if not found)

**Example:**
```python
from exposures import extract_exposure_time

time = extract_exposure_time("light_M31_L_120s_001.fits")
# Returns: 120.0
```

#### `extract_angle(filename)`

Extracts telescope angle from a FITS filename (appears after filter name).

**Parameters:**
- **`filename`** (str): FITS filename

**Returns:**
- **float or None**: Angle in degrees, or None if not found

**Example:**
```python
from exposures import extract_angle

angle = extract_angle("2025-09-14_03-06-28_ALPT_56.26_300.00s_0114.fits")
# Returns: 56.26
```

#### `detect_filter(path)`

Detects filter name from file path or filename.

**Parameters:**
- **`path`** (str): Full path to FITS file

**Returns:**
- **str or None**: Filter name (L, R, G, B, Ha, etc.) or None

**Example:**
```python
from exposures import detect_filter

filter_name = detect_filter("/data/M31/L/light_001.fits")
# Returns: "L"
```

#### `detect_object_name(path)`

Extracts object name from directory path structure.

**Parameters:**
- **`path`** (str): Full path to FITS file

**Returns:**
- **str**: Object name or "Unknown"

**Example:**
```python
from exposures import detect_object_name

obj = detect_object_name("/data/2025/M31/lights/frame_001.fits")
# Returns: "M31"
```

#### `load_flat_history()`

Loads flat frame history from persistent storage.

**Returns:**
- **dict**: {object_name: {filter: set(angles)}}

**Example:**
```python
from exposures import load_flat_history

history = load_flat_history()
# Returns: {'M31': {'L': {45, 56}, 'R': {33}}, 'SMC Panel 6': {'ALPT': {1, 15}}}
```

#### `save_flat_history(history)`

Saves flat frame history to persistent storage.

**Parameters:**
- **`history`** (dict): {object_name: {filter: set(angles)}}

**Example:**
```python
from exposures import save_flat_history, load_flat_history

history = load_flat_history()
history['NGC1234'] = {'Ha': {30, 45, 60}}
save_flat_history(history)
```

### Command-Line Usage Examples

```bash
# Basic analysis
python3 -c "from exposures import calculate_total_exposure; calculate_total_exposure('/path/to/data')"

# With sorting
python3 -c "from exposures import calculate_total_exposure; calculate_total_exposure('/path/to/data', sort_by='total_exposure')"

# Mark flats as captured
python3 -c "from exposures import update_flat_history; update_flat_history({'L': [45, 56]})"

# View all flat history
python3 -c "from exposures import list_flat_history; list_flat_history()"

# Clear all flat history
python3 -c "from exposures import clear_flat_history; clear_flat_history()"

# Clear specific filters
python3 -c "from exposures import clear_flat_history; clear_flat_history(['L', 'R'])"

# Remove specific angles while keeping others
python3 -c "from exposures import remove_flat_angles; remove_flat_angles('L', [45])"

# Reset all flat history (with confirmation)
python3 -c "from exposures import reset_all_flat_history; reset_all_flat_history(confirm=True)"
```

## Configuration

### Supported Filters

The script recognizes these canonical filters in order of preference:
```python
["L", "R", "G", "B", "Ha", "OIII", "SII", "ALPT", "ALPT2", "UVIR"]
```

### Rejection Buffer

The script includes a 15% buffer (configurable via `REJECTION_BUFFER_PCT`) to account for frames that will be rejected during calibration and quality selection. For example, if you target 10 hours, the script will recommend capturing 11.5 hours to ensure 10 hours of usable data after rejection.

## Filter Detection Logic

The script uses a two-step process to identify filters:

1. **Directory-based**: Walks up the directory tree from each FITS file. The first directory name matching a canonical filter is used.

2. **Filename-based**: If no directory matches, it parses the filename (splitting on `_` and `-`) looking for filter tokens.

### Examples

```
/data/M31/L/light_M31_L_120s_001.fits     → Filter: L (from directory)
/data/light_M31_Ha_300s_001.fits          → Filter: Ha (from filename)
/data/NGC1234/OIII/subs/frame_300s.fits   → Filter: OIII (from directory)
```

## Expected File Naming Convention

FITS files should follow this naming pattern for proper exposure time and angle extraction:
```
*_<filter>_<angle>_<exposure>s_*.fits
```

The angle should appear immediately after the filter name.

Examples:
- `2025-09-14_03-06-28_ALPT_56.26_300.00s_0114.fits` (ALPT filter at 56.26°, 300 second exposure)

## Output

The script generates two tables:

### 1. Per-Filter Breakdown
Shows statistics for each filter/exposure combination:
- Filter name
- Exposure time in seconds
- Number of subframes
- Total exposure time (formatted as "H h M m")
- Target time (if specified)
- Remaining time needed
- Remaining subframes needed

### 2. Grand Total Summary
Aggregates all filters showing:
- Exposure times used
- Total number of subframes per exposure
- Combined integration time

### 3. Flat Frame Recommendations (Enhanced!)
Separate table showing:
- **Filter**: Filter name
- **Detected Angles**: Complete range of telescope angles found in your data
- **Existing Flats**: Angles you've already captured (from history)
- **New Flats Needed**: Angles you still need to capture
- **Status**: Visual indicator (✓ Have, ⚠ Need More, ✗ Need All)
- **Coverage Groups**: Shows which existing flats cover which angle ranges (format: `flat_angle→covered_range`)

## Example Output

```
+--------+---------------+---------+------------------+-------------+----------------+-----------------+
| Filter | Exposure (s)  | # Subs  | Total Exposure   | Target Time | Remaining Time | Remaining Subs  |
+========+===============+=========+==================+=============+================+=================+
| L      | 120           | 150     | 5 h 0 m          | 4 h 36 m    | 0 h 0 m        | 0               |
| R      | 180           | 40      | 2 h 0 m          | 2 h 18 m    | 0 h 18 m       | 6               |
| G      | 180           | 35      | 1 h 45 m         | 2 h 18 m    | 0 h 33 m       | 11              |
| B      | 180           | 38      | 1 h 54 m         | 2 h 18 m    | 0 h 24 m       | 8               |
| Ha     | 300           | 25      | 2 h 5 m          | 0 h 0 m     | —              | —               |
+--------+---------------+---------+------------------+-------------+----------------+-----------------+

Grand Total Exposure for All Filters:
+------------------------+---------------------+------------------+
| Exposure Time (seconds)| Total Number of Subs| Total Exposure   |
+========================+=====================+==================+
| 300                    | 25                  | 2 h 5 m          |
| 180                    | 113                 | 5 h 39 m         |
| 120                    | 150                 | 5 h 0 m          |
+------------------------+---------------------+------------------+

Overall total exposure (all filters) → 12 h 44 m

============================================================
FLAT FRAME RECOMMENDATIONS
============================================================
+----------+-------------------------+----------------------+--------------------+-------------+-----------------------------------------------------+
| Filter   | Detected Angles         | Existing Flats       | New Flats Needed   | Status      | Coverage Groups                                     |
+==========+=========================+======================+====================+=============+=====================================================+
| ALPT     | 70.4°-72.5° (16 angles) | 69.6°, 75.0°, 128.0° | —                  | ✓ Have      | 69.6°→70.4°-71.7°; 75.0°→72.5°-72.5°; ?→72.4°-72.5° |
| ALPT2    | 10.0°-72.5° (10 angles) | 10.0°, 69.6°         | 73°                | ⚠ Need More | 10.0°→10.0°; 69.6°→70.5°-71.7°; ?→72.5°-72.5°       |
+----------+-------------------------+----------------------+--------------------+-------------+-----------------------------------------------------+
```

## Flat Frame Management

### Automatic Angle Detection
The script extracts telescope angles from your FITS filenames and suggests optimal angles for flat frames.

### Persistent History
Flat frame history is stored in `~/.astrophotography_flats.json`, tracking:
- Which flats you've captured organized by filter
- The date each flat was captured
- Automatic aging detection (warns for flats >90 days old)

This makes sense since flats are specific to your optical setup (telescope, camera, filter) rather than the astronomical object being imaged. The date tracking helps you maintain fresh calibration files.

### Managing Flat History

#### Adding Flat Frame Records
```python
# Mark flats as captured (organized by filter)
from exposures import update_flat_history
update_flat_history({'L': [45, 56], 'R': [33]})

# Add more flats (keeps existing ones)
update_flat_history({'L': [67, 89]})

# Replace all flats for specific filters
update_flat_history({'L': [45, 56, 67, 89]}, replace=True)
```

#### Viewing Flat History
```python
# View all flat history
from exposures import list_flat_history
list_flat_history()
```

#### Clearing Flat Frame Records

**Clear all filters:**
```python
from exposures import clear_flat_history

# Clear all flat history
clear_flat_history()
```

**Clear specific filters:**
```python
# Clear only L and R filters (keeps other filters)
clear_flat_history(['L', 'R'])

# Clear just Ha filter
clear_flat_history(['Ha'])
```

**Remove specific angles:**
```python
from exposures import remove_flat_angles

# Remove just 45° from L filter, keep others
remove_flat_angles('L', [45])

# Remove multiple angles from Ha filter
remove_flat_angles('Ha', [30, 60, 90])
```

**Reset everything:**
```python
from exposures import reset_all_flat_history

# Safe check (shows warning)
reset_all_flat_history()

# Actually reset all flat history
reset_all_flat_history(confirm=True)
```

#### Common Clearing Scenarios

**Scenario 1: Made a mistake entering angles**
```python
# Remove wrong angle and add correct one
remove_flat_angles('L', [45])  # Remove wrong angle
update_flat_history({'L': [46]})  # Add correct angle
```

**Scenario 2: Re-did flats for a filter**
```python
# Replace all L flats with new ones
update_flat_history({'L': [30, 45, 60]}, replace=True)
```

**Scenario 3: Starting fresh for a filter**
```python
# Clear all flats for the filter
clear_flat_history(['L'])
```

**Scenario 4: Complete reset**
```python
# Nuclear option - clear everything
reset_all_flat_history(confirm=True)
```

### Angle Strategy

The script analyzes telescope angles from your light frames and suggests optimal angles for flat frames using the following strategy:

- **Small angle range (< 10°)**: Single angle at the weighted mean position
- **Medium range (10-30°)**: 3-5 angles distributed across the range
- **Large range (> 30°)**: Up to 5 angles based on frame distribution
- **All angles rounded**: Angles are rounded to whole degrees for practicality
- **Tolerance matching**: Existing flats within ±2.5° of suggested angles are considered adequate

This ensures your flat frames accurately correct vignetting and dust motes at the actual telescope positions used during imaging. The ±2.5° tolerance makes the system practical - you don't need exact angle matches for effective flat field correction.

### Smart Flat Frame Integration (Enhanced!)

The flat frame recommendation system now intelligently considers your existing flat frame inventory:

1. **Coverage Analysis**: First determines which light frame angles are already covered by existing flats (within ±2.5° tolerance)
2. **Gap Identification**: Only suggests new flats for uncovered angle ranges
3. **Efficient Clustering**: Groups uncovered angles to minimize the number of new flats needed
4. **Priority-Based Suggestions**: Recommends angles that provide optimal coverage for the remaining data

**Example**: If you have 70° flats and image data from 70.4° to 73.0°:
- ✓ **Covered**: 70.4°-72.5° data (covered by existing 70° ± 2.5° flats)
- ⚠ **Needs flats**: 72.6°-73.0° data (suggests new flat at ~73°)

This prevents redundant flat frame suggestions and ensures you only capture the flats you actually need.

### Coverage Groups Analysis (New!)

The enhanced flat frame recommendations now include a "Coverage Groups" column that shows exactly which existing flats protect which angle ranges in your data:

**Format**: `flat_angle→covered_range; flat_angle→covered_range; ?→uncovered_range`

**Example Interpretation**:
- `69.6°→70.4°-71.7°` means your 69.6° flat covers data from 70.4° to 71.7°
- `75.0°→72.5°-72.5°` means your 75.0° flat covers the 72.5° data point
- `?→72.4°-72.5°` means the 72.4°-72.5° range needs new flats (uncovered)

This analysis uses the same ±2.5° tolerance applied throughout the toolkit, so you can see at a glance:
- Which of your existing flats are actually protecting your data
- Which angle ranges are left unprotected
- How efficiently your current flat inventory covers your imaging sessions

**Benefits**:
- **Visual Coverage Map**: Instantly see which flats are working for you
- **Gap Identification**: Quickly spot uncovered angle ranges that need attention
- **Efficiency Analysis**: Understand if you have redundant flats or coverage gaps
- **Planning Tool**: Make informed decisions about which new flats to capture

## File Organization Script

### Overview

The `organize_files.py` script helps organize mixed FITS files into filter-specific subdirectories. This is especially useful after imaging sessions where all light frames are saved to a single directory.

### Usage

#### Basic Usage
```bash
python3 organize_files.py "/path/to/mixed/files"
```

#### Dry Run (Preview Mode)
```bash
python3 organize_files.py "/path/to/mixed/files" --dry-run
```

### Examples

#### Organize a night's imaging session
```bash
python3 organize_files.py "/Users/user/Pictures/ASTRO/2025 July/SMC Panel 6 copy"
```

This will:
1. Scan for all FITS files in the directory
2. Detect filter names from each filename
3. Create subdirectories for each filter (ALPT/, L/, Ha/, etc.)
4. Move files to appropriate subdirectories
5. Show a summary of what was moved

#### Preview without moving files
```bash
python3 organize_files.py "/path/to/data" --dry-run

# Output:
# Found 116 FITS files to organize
#
# Detected filters:
#   ALPT: 116 files
#
# ==================================================
# DRY RUN - No files will be moved
# ==================================================
```

### How It Works

1. **Filter Detection**: The script examines each FITS filename to identify the filter
2. **Directory Creation**: Creates subdirectories for each detected filter
3. **Safe Moving**: Moves files with duplicate checking and error handling
4. **Summary Report**: Shows how many files were moved to each filter directory

### Before and After

**Before:**
```
SMC Panel 6 copy/
├── 2025-09-08_03-51-55_ALPT_117.41_300.00s_0081.fits
├── 2025-09-13_20-57-38_ALPT_0.99_300.00s_0096.fits
├── 2025-07-12_01-54-25_ALPT_14.68_300.00s_0021.fits
└── ... (113 more FITS files)
```

**After:**
```
SMC Panel 6 copy/
└── ALPT/
    ├── 2025-09-08_03-51-55_ALPT_117.41_300.00s_0081.fits
    ├── 2025-09-13_20-57-38_ALPT_0.99_300.00s_0096.fits
    ├── 2025-07-12_01-54-25_ALPT_14.68_300.00s_0021.fits
    └── ... (113 more FITS files)
```

### Safety Features

- **Confirmation Prompt**: Asks for confirmation before moving files
- **Dry Run Mode**: Preview what will happen without actually moving files
- **Duplicate Detection**: Won't overwrite existing files in destination
- **Unknown Files**: Files with undetected filters are not moved
- **Error Handling**: Graceful handling of permission issues or other errors

## Session Organization Script

### Overview

The `organize_sessions.py` script takes filter-organized directories and groups them into session directories based on telescope angles. This creates a calibration-ready structure where each session contains all files that can use the same flat frames.

### Usage

#### Basic Usage
```bash
python3 organize_sessions.py "/path/to/filter/directories"
```

#### Dry Run (Preview Mode)
```bash
python3 organize_sessions.py "/path/to/filter/directories" --dry-run
```

### Complete Workflow Example

Here's the complete two-step workflow using "The War and Peace Nebula" dataset:

#### Step 1: Organize by Filter
```bash
python3 organize_files.py "/Users/user/Pictures/ASTRO/2025 July/The War and Peace Nebula copy"

# Output:
# Found 324 FITS files to organize
#
# Detected filters:
#   ALPT: 165 files
#   ALPT2: 159 files
#
# Successfully moved 324 files into 2 filter directories:
#   ALPT/: 165 files
#   ALPT2/: 159 files
```

#### Step 2: Organize by Session
```bash
python3 organize_sessions.py "/Users/user/Pictures/ASTRO/2025 July/The War and Peace Nebula copy"

# Output:
# Available flat frame angles:
#   ALPT2: [70, 75, 82, 10, 71]°
#   ALPT: [70, 75, 82, 71]°
#
# All sessions found in data (±2.0°):
#
# Session_10deg (8° to 12°):
#   ALPT2: 38 files (10.00° to 10.00°) ✓ Flats available
#
# Session_70deg (68° to 72°):
#   ALPT: 141 files (70.40° to 71.67°) ✓ Flats available
#   ALPT2: 104 files (70.47° to 71.67°) ✓ Flats available
#
# Session_71deg (69° to 73°):
#   ALPT: 24 files (72.41° to 72.53°) ✓ Flats available
#   ALPT2: 17 files (72.49° to 72.53°) ✓ Flats available
#
# Created 3 session directories with 324 total files
```

### Final Directory Structure

**Before Session Organization:**
```
The War and Peace Nebula copy/
├── ALPT/
│   ├── 2025-09-08_20-41-17_ALPT_71.53_300.00s_0152.fits
│   ├── 2025-08-01_01-45-43_ALPT_70.40_300.00s_0024.fits
│   └── ... (163 more files)
└── ALPT2/
    ├── 2025-09-05_21-19-04_ALPT2_71.53_300.00s_0071.fits
    ├── 2025-09-12_21-27-30_ALPT2_10.00_300.00s_0145.fits
    └── ... (157 more files)
```

**After Session Organization:**
```
The War and Peace Nebula copy/
├── Session_10deg/
│   └── ALPT2/
│       ├── 2025-09-12_21-27-30_ALPT2_10.00_300.00s_0145.fits
│       └── ... (37 more files at 10.00°)
├── Session_70deg/
│   ├── ALPT/
│   │   ├── 2025-08-01_01-45-43_ALPT_70.40_300.00s_0024.fits
│   │   └── ... (140 more files 70.40-71.67°)
│   └── ALPT2/
│       ├── 2025-09-05_21-19-04_ALPT2_71.53_300.00s_0071.fits
│       └── ... (103 more files 70.47-71.67°)
└── Session_71deg/
    ├── ALPT/
    │   ├── 2025-07-27_18-44-30_ALPT_72.41_300.00s_0001.fits
    │   └── ... (23 more files 72.41-72.53°)
    └── ALPT2/
        ├── 2025-08-01_20-12-22_ALPT2_72.49_300.00s_0017.fits
        └── ... (16 more files 72.49-72.53°)
```

### How It Works

1. **Flat Frame Analysis**: Reads your available flat frame angles from `~/.astrophotography_flats.json`
2. **Angle Extraction**: Extracts telescope angles from FITS filenames
3. **Smart Clustering**: Groups files using available flat frame angles with ±2.5° tolerance
4. **Priority System**: Prefers flat angles that:
   - Cover the most image data
   - Are available for the most filters
   - Are closest to the center of the data range
5. **Session Creation**: Creates session directories and moves filter subdirectories into them
6. **Flat Recommendations**: Shows which additional flats need to be captured

### Key Features

#### Intelligent Flat Frame Integration
- Uses your existing flat frame inventory to determine session grouping
- Only creates sessions for angles where you have (or need) flat frames
- Prioritizes angles that cover the most data efficiently

#### Smart Angle Clustering
- Groups files within ±2.5° tolerance for practical flat frame application
- Avoids creating redundant sessions for very similar angles
- Prefers existing flat frame angles over creating new angle requirements

#### Complete Calibration Workflow
Each session directory contains all files that can use the same flat frames:
- All filters are grouped by compatible telescope angles
- Each session shows flat frame availability status
- Missing flat frames are clearly identified for future capture

### Flat Frame Recommendations

The script provides specific guidance on which flat frames you still need:

```
============================================================
FLAT FRAME RECOMMENDATIONS
============================================================
The following flat frames need to be captured:

📷 ALPT filter:
   • Capture flats at 10° telescope angle

💡 Tip: Capture flat frames at these angles to enable full calibration
   of all your image sessions.
```

### Session Benefits

**For Image Processing:**
- Each session can be calibrated with the same flat frames
- Reduces the number of master flats needed
- Ensures consistent vignetting correction across filter sets

**For Data Management:**
- Logical grouping by telescope configuration
- Easy identification of calibration requirements
- Clear separation of data requiring different flat frames

### Integration with Exposure Tracker

The session organization integrates seamlessly with the flat frame management:
- Reads flat frame history from the same cache file
- Updates recommendations based on your current flat inventory
- Maintains the same ±2.5° tolerance used throughout the toolkit

## Tips for Best Results

1. **Organize by Filter**: Place FITS files in directories named after their filters for most reliable detection
2. **Consistent Naming**: Use a consistent naming scheme with exposure time in seconds (e.g., `_300s_`)
3. **Plan Ahead**: Use the target planning feature to ensure balanced filter coverage
4. **Account for Weather**: The rejection buffer helps, but consider adding extra margin for weather losses
5. **Regular Monitoring**: Run the script periodically during imaging sessions to track progress

## Troubleshooting

- **No files found**: Ensure your FITS files have the `.fits` extension (case-sensitive)
- **Missing exposure times**: Check that filenames include `_<number>s_` pattern
- **Unknown filters**: Either organize files in filter-named directories or ensure filter names appear in filenames
- **Percentages > 100%**: While allowed, this will generate a warning. Check your `filter_pct` values

## License

This script is provided as-is for the astrophotography community. Feel free to modify and distribute as needed.

## Contributing

Suggestions and improvements are welcome! Common enhancement areas include:
- Support for additional filter types
- Integration with FITS header reading for more accurate metadata
- Export to CSV/Excel formats
- GUI interface for easier use