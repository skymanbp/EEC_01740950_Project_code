# Growth Curve Data Format

This folder holds the raw growth curve data: 101 tab-separated plate files covering 18 of the
21 strains (85, 88, 100, 186, 322, 333, 350, 353, 374, 380, 390, 398, 436, 442, 448, 487, 527,
565). AUC for the remaining three strains (74, 331, 371) comes pre-computed from
`data/spline-fits.csv`.

## File Naming Convention

```
{strain}_p_{library}_r_{replica}.txt
```

**Examples:**
- `88_p_1_r_1.txt` - Strain 88, Library 1, Replica 1
- `100_p_2_r_3.txt` - Strain 100, Library 2, Replica 3

## File Format

- **Format**: Tab-separated values (.txt)
- **Header row**: Column names
- **Columns**: `Time`, `T° 600`, A1, A2, ..., H12 (96 wells; 98 columns in total)
- **Encoding**: Latin-1 — the degree sign in the `T° 600` header is the single byte `0xB0`

### Sample Data Structure

First and last rows of `88_p_1_r_1.txt`, middle columns elided:

```
Time	T° 600	A1	A2	A3	...	H12
0:00:06	25.0	0.043	0.042	0.054	...	0.046
1:00:03	24.9	0.043	0.043	0.055	...	0.047
2:00:05	24.7	0.045	0.044	0.056	...	0.048
...
72:01:10	30.5	0.590	0.616	0.678	...	0.573
```

## Notes

- **Time format**: HH:MM:SS (converted to hours in analysis)
- **`T° 600`**: Temperature logged alongside every OD600 row, recorded and then dropped by
  `remove.temp()` (`Project_Main_Code.R:47-51`), which selects columns `c(1, 3:98)` by
  position, not by name
- **Wells**: 96 wells in standard microplate format (A-H rows, 1-12 columns)
- **OD600**: Optical density readings at 600nm
- **Duration**: 72 hours at 1-hour intervals — 73 readings per plate, except strains 353 and 527 which have 65 (a 9-hour gap between the 6 h and 15 h readings, i.e. 8 missing readings)

## Required Strains

**Main strains:**
88, 100, 186, 322, 333, 350, 353, 374, 380, 390, 442, 448, 487, 527, 565

**Special-format strains:** 436, 398, 85

Their 11 files sit flat in this folder alongside the rest, named
`_{strain}_p_{library}_r_{replica}_2.txt` (leading underscore, `_2` suffix): 4 files for 85,
4 for 398, 3 for 436. `Project_Main_Code.R:126` reads them from a separate directory
(`setwd("D:/Masters_project/data/weird")`), so running the script means copying just these 11
files somewhere and pointing that line at it — there is no `weird/` directory in the repo.

**That destination path must contain exactly one underscore before the filename.** Line 127
builds absolute paths (`full.names=TRUE`), line 131 splits each one on `_`, and lines 138-140
take strain, library and replica from split positions 3, 5 and 7 — an alignment that holds
only because `Masters_project` contributes exactly one underscore ahead of the filename. Any
other count shifts those positions and labels strains 85, 398 and 436 from the wrong fields,
silently rather than erroring.

**Extra strains (from spline-fits.csv):**
331, 371, 74
