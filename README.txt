Manuscript Title :
Efficient analytical partial derivatives of modal phase velocity with respect to layer parameters. 

Authors:
Prabir Das and Tarun Naskar


MATLAB Code : Input.m

Input: Soil parameters (S-wave velocity, P-wave velocity, density, and thickness) of a layered profile need to be read from an Excel file (Line No. 19).
       The minimum and maximum frequency, frequency and velocity resolution, and the mode for the expected Jacobian need to be provided. If "mode = 1" (Line No. 29), it indicates the fundamental mode.
       The user can also choose a higher number for higher-mode's Jacobians. For example, "mode = 2" gives the Jacobian of first higher mode.	


How to Run the Code: (1) Open the script "Input.m" in MATLAB.
                     (2) Specify the values of input parameters in the "INPUT" section.
                     (3) Now run the script.  

Output: The Jacobian of phase velocity will be calculated and plotted in the figures.