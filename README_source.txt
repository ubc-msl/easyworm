SOURCE CODE for EASYWORM


1. Copy all the files in Matlab workspace to get it working, 
	+ the folder 'Hist_gfit -ADD TO PATH' (THx to Dr. Thomas L. Gaussiran IIgauss@arlut.utexas.edu for his script)

2. Add the aboved-mentioned folder to the working directory (ie Matlab workspace)

3. Make sure you have following toolboxes installed or some Easyworm functionalities may not be working:
	-Curve Fit toolbox
	-Image Processing toolbox
	-Statistics toolbox
 

4. Supplementary scripts:

The scripts 'printo.m' ; 'savetop.m' ; 'synseries.m'  have been included 
as they complete the analysis tool available from the main executables.  
However they can be ran from Matlab workspace only (not from Easyworm interfaces). 


5. Outputs:

To get (some) variables generated in the Matlab workspace instead of or in addition to the .txt outputs, 
look in the code for all save() functions, at the bottom of which you will find the appropriate code.
In the current version it is commented out, so simply uncomment it, re-compile and it should work.