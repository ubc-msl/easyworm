# easyworm
## Easyworm - Windows 10 64-bit and Linux 64-bit

Easyworm is a very cool and easy-to-handle software tool that analyzes the shape fluctuations of semi-flexible polymers and provides accurate measurements of their persistence lengths. A complete set of additional tools is provided, including the derivation of the Young’s modulus.

### INSTALLATION NOTES

To use Easyworm, it is necessary that you install the Matlab Compiler Runtime (MCR) first (except if you have Matlab already installed and plan to work with the source code only). Make sure that you install the version that corresponds to your operating system, currently supported operating system is Windows 10 with Matlab 2020a or Matlab 2012b for Linux64.

#### STEP-BY-STEP INSTRUCTIONS (INSTALL AND EXECUTE)

1. Download the source or binaries for Windows 10 64-bit, requires both the MCR and the Easyworm software suite.

Windows & Linux 64-bit:

Easyworm recompiled with Matlab 2020a 64bit under windows 10 and Matlab 2012b 64bit for Linux.

* MCR 2020a 64 Bit: [MCRInstaller-win64-2020a(https://ssd.mathworks.com/supportfiles/downloads/R2020a/Release/5/deployment_files/installer/complete/win64/MATLAB_Runtime_R2020a_Update_5_win64.zip)]
* MCR 2012b 64 Bit: [MCRInstaller-linux64-2012b (https://ssd.mathworks.com/supportfiles/MCR_Runtime/R2012b/MCR_R2012b_glnxa64_installer.zip)]

Downloading the MCR directly on Mathworks’s website. However, since backwards compatibility issues are highly probable, it is very strongly advised to use only the versions as indicated with which Easyworm have been compiled (i.e. those available here, or equivalent versions on Mathworks website, and not the most recent ones).

2. Clone or download the source for Github.  How to tutorial: https://help.github.com/en/articles/cloning-a-repository

3. Install the MCR first.

Under Windows, all you need to do is run the MCR installer (you need administrative rights).
Under Linux, you should not need to login as root or be in super user mode, just cd to the uncompressed MCR_installer folder and type ./install in the command line, that should work. Then follow the instructions that come to screen. If you get the message ‘permission denied’ after typing ./install, that probably means you downloaded the file or copied it from a Windows machine. In this case, re-download the MCR_installer and try to install it, but doing all steps under Linux system only.
4. Run Easyworm.

Under Windows, all you have to do is double-click on the executable file (e.g. Easyworm1.exe).
5. Both the source code for Easyworm and some example files are included in the Easyworm software package:

The source code will work only under the Matlab environment.
Example files include images of 3 different sets of polymers (2 experimental and 1 synthetic) with 3 distinct persistence lengths: 52 nm, 1.5 um, and 100 um. These can be used for testing and training purposes; that is, if you analyze the same sets and do not find a similar persistence length, it means you made a mistake during the analysis.
 

### CONTACT

If you encounter any trouble to download the files, install Easyworm, or just to use it, please send an email to: lamour99 [at] hotmail.com and I will try to help.
