# Metadata File Format

## Metadata format

**Header Format**
The top three rows must contain:
1) Machine-readable keywords, eg. "probe_area", that will become dictionary keys
2) Unit string, or blank for dimensionless units, eg. "mm^2"
3) Human-readable title or note, eg. "Bdot Area"

**Types of Metadata files**
Metadata files are sorted into four types automatically based on whether the csv file includes the "run" or "probe" keywords.

Type |Has "run" key?|Has "probe" key?|Explanation
-----|--------------|----------------|-----------
Experiment Data|No|No|Data that pertains to the whole experiment, eg. experiment name, vacuum chamber used
Run Data|Yes|No|Data that pertains to a given run number, eg. background field, fill pressure
Probe Data|No|Yes|Data about a particular probe for the whole experiment, eg. calibration constants
RunProbe Data|Yes|Yes|Data about a probe on a particular run, eg. positon or attenuation



## Key Dictionary

**Experiment Data**
No standard keys exist for this type of file yet

**Run Data**
Keys in this section are all optional metadata (not used in processing probe data).

**Probe Data**
