# Growth Curve Data Format

This folder should contain growth curve data files.

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
- **Columns**: Time, Temperature, A1, A2, ..., H12 (96 wells)

### Sample Data Structure

```
Time	Temperature	A1	A2	A3	...	H12
0:00:00	37.0	0.05	0.06	0.05	...	0.05
0:30:00	37.0	0.08	0.09	0.08	...	0.07
1:00:00	37.0	0.12	0.14	0.11	...	0.10
...
72:00:00	37.0	1.25	1.32	1.18	...	0.95
```

## Notes

- **Time format**: HH:MM:SS (converted to hours in analysis)
- **Temperature**: Recorded but removed during analysis
- **Wells**: 96 wells in standard microplate format (A-H rows, 1-12 columns)
- **OD600**: Optical density readings at 600nm
- **Duration**: Typically 72 hours with 30-minute intervals

## Required Strains

**Main strains:**
88, 100, 186, 322, 333, 350, 353, 374, 380, 390, 442, 448, 487, 527, 565

**Special strains (in 'weird' subfolder):**
436, 398, 85

**Extra strains (from spline-fits.csv):**
331, 371, 74

## Missing Data Notice

The original growth curve data files need to be added to this folder for the analysis to run successfully.
