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

"datafile" (str): Datafile name. HDF/datafile.hdf5 should be the datafile for each "datafile" in this column.

**Probe Data**

All probes:

"probe_type" (str) -> Specifies the probe type, used to decide which analysis routines to run. Example: "bdot", "tdiode"

Bdot:

"nturns" (int) -> Number of bdot turns (assume to be the same for all axes)

"{xyz}area" (float) -> Area of the probe tip from calibration

"{xyz}tau" -> High-frequency calibration constant for each axis

Langmuir:
"area" (float) -> Area of the probe tip (used when calculating isat density)


**Run-Probe Data**

All Probes:

"{xyz}pol" (1 or -1): This factor is multiplied by the data that is read in to potentially reverse it if a probe was in upside down.

"gain" (float): Amplifier gain prior to the digitizer.


Probes Using LAPD digitizer:

"digitizer" (str): Name of the digitizer used. Example "SIS crate"

"adc" (str): Name of the analog-to-digital converter used. Example "SIS 3305", "SIS 3302"

"brd{i}" (int): Digitizer board used, where {i} is the number of the channel. There should be one of these columns for each channel.

"chan{i}" (int): Channel on the digitizer used, where {i} is the number of the channel. There should be a coorresponding "brd{i}" for each one.

Probes Using HRR Digitizer:

"resource{i}" (int): Resource number for each channel {i}

"chan{i}" (int): Channel number. The number of these columns should match the number of resource columns.


All Probes with Position Information:

"probe_origin_{xyz}" (float): Position of the probe origin relative to the experiment coordinate system. These will be added to all the probe positions.

"{xyz}pos" (float): Position of the probe. Overridden by motor drive information if a probe is being scanned. These positions are relative to the probe origin.

"roll" (float): Angle, in degrees, that the probe was rotate about its central (x) axis. This is included during the rotation correction phase, and can be used to correct a probe that was misaligned.


Probes Using LAPD Motor Drives:

"motion_controller" (str): Which LAPD drive was associated with the probe. Example "6K Compumotor", "NI_XYZ"

"motion_receptacle" (int): Which instance of that motor drive was used? Starts at 1. For example, 6K Compumotor can control four XY drives, labeled 1,2,3,4

"rot_center_{xyz}" (float): Required only for probes that rotate on a ball valve (LAPD probe drives). This specifies the center position of the ball valve, for use in angle corrections.


Bdot Probes:

"{xyz}atten" (float): Attenuation on each axis of the probe. Required to be in dB currently?

Langmuir Probes:

"atten" (float): Attenuation on the digizer.

"ramp_atten" (float): Attenuation for ramp, only required for vsweep.

"ramp_gain" (float): Gain on ramp signal: only required for vsweep. Regular "gain" keyword is used for the other signal channel.

"sweep_type" (str): Currently "langmuir_isat" or "langmuir_vsweep" are used. This keyword is just useful for deciding which type of analysis routine to call on this dataset.

"resistor" (float): Measurement resistor for vsweep runs.

"bias" (float): Probe bias for isat runs.



