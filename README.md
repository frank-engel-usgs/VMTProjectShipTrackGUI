VMT
===

The Velocity Mapping Toolbox: Project ShipTrack Mod

A Matlab-based program for visualization and analysis of moving boat ADCP data. This is a modified version of the full VMT v4.09 that enables processing of TRDI ADCP data without GPS. See the full instructions and caveats in the file [https://github.com/frank-engel-usgs/VMTProjectShipTrackGUI/blob/master/VMT%20Project%20ShipTracks%20Instructions.pdf](VMT Project ShipTracks Instructions). Absolutely *no support* is given. This is a custom one-off tool I made for a specific project. It has been shared with other people, with some sucess, but I unfortunately cannot support bug or feature requests for this modification at this time. 

Notes
-----
I have compiled a tool that may help you to project the ADCP bottom track data into real-world (UTM) coordinates. I offer no guarantees that the method will work for you, or that the method contains less than an acceptable level of error. Basically, I just don't condone it. That said, it may still be of use you to you and others on this forum. Below I explain my method:

```
% To locate the position of each measured ensemble within a transect, a
% method based on ADCP bottom track data that accounts for positional
% errors has been developed. The positions of the ADCP probe at the
% beginning and end of each pass is determined relative to the surveyed
% endpoints of each cross section, accounting for the offset between the
% endpoint and the probe. The position of each ensemble along the transect
% is then computed by distributing the probe’s reported values of linear
% distance traveled (dead reckoning), or “distance made good” (DMG),
% between the beginning and ending locations of the pass. In high flows, an
% actively mobile bed may exist, producing error in the ADCP bottom track
% positions.  In these cases, length scaling adjustments are made to
% ensemble locations by distributing the difference between the measured
% distance and probe-reported distance equally among all ensembles.
```
The file "VMT_ProjectShipTrackGUI.exe" file runs using the same Matlab MCR that VMT uses. I highly recommend running the file from a command/bash window, as the code will output important information to this window as it runs. If you do not run the code from the command line, you will not see these data. When VMT_ProjectShipTrackGUI.exe opens, it looks like this:

IMAGE HERE

1. Enter the Lat/Lon coordinates of the endpoints for your cross section in the 4 boxes at the top
2. Use the "Select Files" button to choose the ASCII data from WinRiver containing BOTTOM TRACK referenced data. You can load multiple ASCII files corresponding to ONE cross section
3. In the table, insert the distance from the ADCP center to your known endpoints in the Left/Right Offset columns
4. Enter the starting bank
5. Press the "Process Transect" button

This action will create a VMT v4.xx compliant MAT file that can be loaded into the program normally. As the program runs, it will send some text info to the command window:

IMAGE 2 HERE

Be sure to look at the data below "Transect pass X". The code compares the dead reckoning distance (Distance made good, or DMG) reported by bottom track against the estimation method. Large errors can result from this process, especially if the bed is mobile, or more than 5% of the bottom track data is bad. I would not trust this method in cases where the Percent error reported are greater than 5% without a firm understanding of the impacts on your data. Note that:

1. This method creates and corrects for ensemble coordinate positions, it does nothing to the velocity or discharge data
2. This means that if there is a moving bed, you will bias your velocities and discharges low (water velocity = what the ADCP measures as velocity minus the speed of movement of the ADCP).

If you can't tell, I am wary of this method. I am posting this because there are several other users who have requested support on this issue. I realize and acknowledge that real world conditions often lead to situations like this. PLEASE USE THIS METHOD WITH DISCRETION.

Disclaimer
----------
This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey  (USGS), an agency of the United States Department of Interior. For more information, see the official USGS copyright policy at [http://www.usgs.gov/visual-id/credit_usgs.html#copyright](http://www.usgs.gov/visual-id/credit_usgs.html#copyright)

Although this software program has been used by the USGS, no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.

This software is provided "AS IS."

This information is preliminary and is subject to revision. It is being provided to meet the need for timely
best science. The information is provided on the condition that neither the U.S. Geological Survey nor the
U.S. Government may be held liable for any damages resulting from the authorized or unauthorized use of
the information.
