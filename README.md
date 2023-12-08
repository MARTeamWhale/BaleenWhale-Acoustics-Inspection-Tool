# BaleenWhale Acoustics Inspection Tool (BAIT)

## Overview
This project is primarily a MATLAB tool for reviewing baleen whale call detections, notably those output by the *Low Frequency Detection and Classification System* (LFDCS)<sup>[1]</sup>. It was originally written for validating North Atlantic right whale upcalls, but can be used for other baleen whale calls as well.


## Installation
1) Install the *MUCA* library from the MARTeamWhale GitHub account (https://github.com/MARTeamWhale/MATLAB-Utilities-for-Cetacean-Acoustics) 
2) Clone or download the BAIT repository in its entirety and save it somewhere on your computer
3) Add the location of the cloned repository folder to your MATLAB path


## Usage
The project consists of three scripts:

- ***BAIT_validate*** - This is the main tool. Launches a graphical user interface (GUI) application from which users can review and validate baleen whale call autodetections, mark possible and definite missed calls, or general annotations for other signals of interest.

- ***BAIT_convertLFDCSTable*** - Converts LFDCS autodetection CSV tables into XLSX files in a format that can be read by BAIT.

- ***BAIT_isolateLFDCSDetections*** - Finds the locations of LFDCS detections within the original audio files and generates short clips and/or spectrogram images from those detections.

The scripts are designed such that **there is no need to edit any code** to run them. Certain input arguments may be specified via parameter files (see the *PARAMS* folder).

More information on usage is available through the following resources:

- Script file headers (i.e. comment block at the top of the file)

- *SUMMARY_OF_IMPORTANT_CHANGES.pdf*

- *MERIDIAN_Detection_Browser_ReferenceManual_20211028.pdf* (OUT OF DATE, but still relevant for understanding how to use the GUI)


## References
1) Baumgartner, M. F., & Mussoline, S. E. (2011). A generalized baleen whale call detection and classification system. *J. Acoust. Soc. Am.* **129**, 2889-2902.