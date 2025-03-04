# Ghandour2025_NatureCommunications

# 1. Install MATLAB.

# 2. Download the "Khaled_et_al_coactivity - without ca-data" folders.
   Aim:

Using custom-made Matlab codes to generete the following:

First_code: Extract, filter, & plot the calcium transients over time
Second_code: calculate the cooccurence of neurons from different categories.

Prerequisites:
Software: "MatlabR2018a" available http//www.mathworks.com
Instalation guide: https://jp.mathworks.com/videos/how-to-install-matlab-1525083586145.html
Installation time: approx. 30 min
Operating system: Windows 10 pro
Non-standard hardware: https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/support/sysreq/files/SystemRequirements-Release2018a_Windows.pdf
Data files: 1. Calcium raw data over time ("finalv_MouseID.CSV" file produced by "Hotaru" Automatic cell detection) 
	    2. "Frames for session.xlsx" file containing the timing of each session.

Steps:

Load 'v109entobe.mat'manually.

First_code: "extract_ca_transients_plot.m"
a. Open the software and set the directory to the folder containing all the codes and files.
b. Open the finalv_mouseID.CSV file to view the fluorescence data table then import and save as a numeric matrix file (finalv_mouseID.mat) in the same folder.
c. Write the name of the imported data file before loading the matlab numeric matrix file (finalv_mouseID.mat) by double clicking.
d. Open the Matlab code file "first_code_extract_ca_transients_plot.m" and Input the range (sec) of all sessions' time (from the start of the first session till the end of the last session) provided in the "Frames for session.xlsx" files.
e. From the above tool bar press "Editor" then "Run" to remove low frequency flactuations and noise by using a 0.01 Hz high-pass filtering and calculating the z-scores to sxclude those below 3 SD. "Running time" is approximately less than 30 sec. 
f. expected run time (couple of seconds)


Second_code:"cooccurence":
a. After runing all the steps in the first code, open the third code that is used to calculate the synchronized neuronal activity.
b. Input the classified neurons (results of the second_code) into Cell_type_A, Cell_type_B, or Cell_type_C, respectively.
c. Specify the value for the bin_Frame_Num as "4" to calculate the co-active neurons for each 250 ms time bin.
d. From the above tool bar press "Editor" then "Run"."Running time" is approximately less than 30 sec.
e. Co-activity data will be displayed in the command window in the folowing order: (Engram-to-be, Common engram cell, Specific engram, Engram-to-be + Common engram, Engram-to-be+ Specific engram, Engram-to-be+Common engram+Specific engram) in addition to the generated plots.
f. expected run time (couple of seconds)


# 3. Download the "Khaled_et_al_PCAICA - without ca-data" folders.
Aim:

Using custom-made Matlab codes to generate the following:

Extract, filter, & plot the PCAICA-detected patterns in a certain session and backproject them over the whole session to determine their activity in other sessions.


   Steps:

Load 'calcium_135.mat'manually.

a. Open the software and set the directory to the folder containing all the codes and files.
b. run the code and press y (as a yes for the question that arises)
c. all figures and data plots will appear automatically.
d. expected run time (couple of seconds) 
# Ca2+ data for running the above codes are uploaded at the following link. Please download the files from the link below.
https://doi.org/10.5281/zenodo.14963087


