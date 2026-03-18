Mean Field Approximation SPC
----------------------------
Information about this method can be found here: http://www.weizmann.ac.il/complex/compphys/group_papers1/omer_mscthesis_final.pdf

You can either use the code or Stand-alone software for Windows:

For support, please contact: Assif.Yitzhaky@weizmann.ac.il
__________________________________________________________
Code: the SPC C code is suitable for Windows, Linux and Mac.

Makefile can be used to compile on Linux)

Note: You should run the MF version like the Monte-Carlo version here: https://github.com/eytandomany/SPC
You can find documentation in the following files: README.PDF and document.pdf.
In the following folders you can find examples for input and output files: 3Conc and ImageDist .
_________________________________________________________________________________________________

Alternatively, you can use a friendly standalone GUI software instead of the code above. 

Stand-alone software:

One-way SPC software for Windows (64-bit) can be downloaded from:

http://www.weizmann.ac.il/complex/compphys/software/ctwc/spc_app_v16_030313.zip
Help menu is included.

Instructions:
1. Extract to some folder.
2. Enter into SPC folder and launch run.bat

IMPORTANT NOTE:  If Matlab R2011b (7.13) is not installed on your computer, you will have to download and install MCR (Matlab Component Runtime) for this version from
http://www.weizmann.ac.il/complex/compphys/software/mcr/mcr713/MCRInstaller.exe
Note on SPC:
This software can utilize either:

1) The classical Monte-Carlo SPC (default).
2) The Mean Field Approximation SPC which is faster and deterministic. The larger is the number of points to cluster, the better is the approximation.

Note: The input files must be tab-delimited (see also help menu and example).

 
