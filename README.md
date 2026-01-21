Contained herein are data and R Code for carrying out several types of Global Sensitivity Analysis (GSA) on deterministic simulators. GSA methods include Morris Elementary Effects, Sobol' sensitivities, VARS-TO variogram method, and several regression-based methods.

Currently two deterministic simulators, GR6J and SimplyP are included.  

Note: The entire R Markdown files need to be knit, not separate knitting of individual chunks.

Easy way to install: The easiest way to run the code is to download and unzip the two compressed files, GR6J.zip and SimplyP.zip, and unzip them into a folder on your computer. The default location is drive D, but this can be changed by editing the R Code files, changing the value of the "root" object (in SimplyP_GSA_generate.Rmd and SimplyP_GSA_evaluate.Rmd, or GR6J_GSA_generate.Rmd and GR6J_GSA_evaluate.Rmd). 

More tedious way: Alternatively, the individual files and folders within GR6J and SimplyP can be downloaded.For each simulator user needs to create a folder, GR6J, then 3 subfolders: R Code, Data, and Output. Likewise for SimplyP.  

In both GR6J/Output and SimplyP/Output, create 6 more subfolders: 1_Morris, 2_Sobol, 3_VARS.TO, 4_Regression, 5_RegTree_RF, 6_GPReg. For SimplyP/Output add a 7th folder: 0_Input_Data_plots.

Put the R Code in he R Code Folder, one for generating output from the simulator (...generate.Rmd) and one for calculating the GSA measures.

Ine Data file put the data files for GR6J this is Coull_GR6J_R.xlsx. For SimplyP they are Tarland_data_matrix.RData, TarlandInputs.dat, and TarlandParameters_v0-4.dat.

Then use RStudio to open up the R code files. Edit the lines in the two R code files, e.g., GR6J_GSA_generate.Rmd (~ line 60) and GR6J_GSA_evaluate.Rmd (~ line 43), where the object root has your GR6J folder location.  To run Morris method, simply knit GR6J_GSA_generate.Rmd first, then knit GR6J_GSA_evaluate.Rmd.  Results will appear in the Output subfolder. Read the text at the top of the R code files for information about specifying the GSA method.


