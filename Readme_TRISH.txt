This readme contains files needed by ReconAnalog as well as sample data and files useful for understanding ReconAnalog and its use by TRISH
        
=========== FILES

ReconAnalog.R -- the main script, which read json file and calls many functions
*.R --- functions ReconAnalog depends on
PackagesNeeded.txt --- R packages use must install 

Input data files used by Meko running ReconAnalog outside of TRISH
        TreeData.txt
        TreeMeta.txt
        HydroData.txt
        
inputTRISH.xlsx ---spreadsheet we have talked about as form in which user could submits tree-ring time series and metadata to UNH. Data in TreeData.txt and TreeMeta.txt are in separate sheets of this file.. 

init01.json -- file with the initial settings of file names, folder names, specifications, etc., for initial default run of ReconAnalog. User would typically change the settings after that initial, and perhaps build a new file (say, init02.json) to hold settings for current run of ReconAnalog.
        
ScreenInterfaceWithReconAnalog.pdf -- describes the contents of the above json file and how the settings link to screen 1 and screen 2 of TRISH

Files describing the four reconstruction methods available with ReconAnalog. When ReconAnalog is run, the appropriate pdf is one of the files written to the system test_out folder. Meko makes these four pdfs available to TRISH from folder designated in the json init file 
        TrishOutputDescribeAnalog.pdf
        TrishOutputDescribeMLR1-noPCA.pdf
        TRISHvisual/TrishOutputDescribeMLR1-PCA.pdf
        TrishOutputDescribeSLR1.pdf

TRISHvisual/AbbreviationsTRISH.pdf -- list of abbreviations frequently used in the above pdfs 

TRISHoverview.pdf:  Brief overview of the online reconstruction tool TRISH (Tree-Ring Integrated System for Hydrology)
        
Two mock-up TRISH screens designed by Meko. Suggested screens for user input to TRISH for using ReconAnalog
        Screen1_TRISH.pdf
        Screen3_TRISH.pdf


ListFilesTRISH01.txt: (not important to UNH) --- list that I use to copy all the files to be zipped together to send to UHN 

=========== Meko laptop set up for running ReconAnalog
        
1) ReconAnalog, TreeData.txt, TreeMeta.txt, HydroData.txt are in an R current working directory
2) User-written function ReconAnalog depends on are in another folder (designated in init01.json)
3) the four "Describe" pdfs (see above), which must be accessed by ReconAnalog are in a separate 	system folder
4) I can run ReconAnalog in Rstudio, or alternatively from  using "Rscript" from the Linux comman prompt, as long as I am in the current working directory
5) All output goes to a system folder "test_out"