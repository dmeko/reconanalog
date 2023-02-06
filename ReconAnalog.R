################################################################################
#
# ReconAnalog.R
# Runoff reconstruction by analog and non-analog methods, all two-stage
# D Meko
# Last revised 2022-12-15
#
# This script, as written, can be run from the R prompt on a linux system. The 
# script, with modification, can be run from from the TRISH (Tree-Ring Integrated 
# System for Hydrology) developed by a team of researchers from the University
# of New Hampshire (UNH) and University of Arizona (UA)
#
# ReconAnalog requires a JSON input file (e.g., init01.json) with reconstruction
# specifications that be tailored by the user. Other reconstruction specifications
# are hard-coded into the script, and are not meant to be accessible for setting
# to TRISH users.
#
# ReconAnalog can do reconstructions alternative by four different methods, which
# are described in detail in pdf files copied by ReconAnalog to a system folder
# when running ReconAnalog. Please refer to the appropriate pdf for details. 
#
#================ TREE RING INPUT
#
# ReconAnalog requires input time series of tree-ring chronologies and a predictand
# seasonal hyrdrologic or climatological variable. ReconAnalog also requires metadata
# for the chronologies and specifications for the reconstruciton. How ReconAnalog gets
# all this data depends on whether ReconAnalog is run within TRISH or from a laptop's
# R prompt. TRISH culls the tree-ring data from its internally stored time series and
# metadata, and prepares the predictand seasonal hydroclimatic predictand from monthly
# data via the UNM water balance model. Outside of TRISH, the user supplies the data
# files of tree-ring chronologies and metadata, and the time series of predictand 
# hydroclimate data in tab-separated files with format described below. For example, 
# the tree-ring time series matrix and metadata read by statements starting:
#   "U <- read.table('treeData.txt',sep="\t",header=TRUE)
#   "Tmeta <- read.table('treeMeta.txt',sep="\t",header=TRUE) 
#
# The tree-ring time series data provided by TRISH will already have been screened for
# time coverage and geographic domain. Outside of TRISH, all chronologies in 
# the provided file of chronologies will be used. The time series and metadata used
# outside of TRISH should have the following form.
#
# 1) Time series. Tab-sep matrix of chronologies. Data should have one header
#   row of site codes (following format rules, and with "Year" as first-col header). 
# 2) Metadata. Rows should match columns (after "Year" column) of the time series
#   matrix and should be tab-separated. Columns should be in this order:
#   1 sequential number (1 to however many are in the time series matrix)
#   2 site number corresponding to the column of the chronology in the network the
#     the user uploaded to UNH
#   3 site code (1-12 characters, no spaces, and maybe some other rules we should specify)
#   4 longitude east (decimal degrees; negative if west)
#   5 latitude north (decimal degrees; negatve if southern hemispher)
#   6 elevation (m above msl)
#   7 species code (4-letter code, following ITRDB convention)
#   8 data-type (1 letter code): R=total ring width, E=earlywood width, L=latewood width,
#     X=maximum density. 
#   9 First and last year that chronology had valid data.
#   Here is and example line:
#       7	 45	BUT    	  93.367000	  64.283000	  113	LAGM	R 1723   1999
# 
#
#
############# JSON INPUT FILE
# 
# ReconAnalog was revised in Nov 2022 so that user-changeable inputs are no longer modified
# by changing lines of code, but by changing setting of an input JSON file (e.g., init01.json).
# The JSON file has 23 inputs, described in some detail here. Default settings from Meko's 
# trial run outside of TRISH on a Linux laptop are included.
#
# "code_dir" :  "/home/dave/Data/RlibraryMeko/",
#   Directory with many user-written R function that ReconAnalog must be able to access
# "pdf_dir" :  "/home/dave/Data/RlibraryMeko/",
#   Directory with the four pdfs describing the alternative reconstruction methods.
#   Pdf files are copied from this directory to the output directory on the system
# "tr_file" : "treeData.txt",
#   Name of the tab-separated file with time series of tree-ring chronologies
# "trM_file" : "treeMeta.txt",
#   Name of file with tree-ring metadata
# "cl_file" : "hydroData.txt",
#   Name of file with time series (year and value) of predictand climatic or hydroclimatic predictand
# "outputDir" :  "/home/dave/AAAtrish2/test_out/",
#   The system folder to which the output will be written
# "NameNetwork":  "Kyzyl",
#   The name of the tree-ring network. In TRISH, use picks this from a dropdown menu
#   at screen 1.
# "PrewhitenOrder" : 0,
#   Whether or not to prewhiten chronologies with AR model before use in regression modeling.
#   If "0", do not prewhiten. I "1", "2", or "3," prewhiten with that order of AR model.
# "LallowLags" : true
#   Whether to allow lags in the SSR models. If "false," only lag-0 models are allowed.
#   If "true," lags t-2 to t+2 from the predictand are allowed in pool of potential predictors
# "NsitesUserNetwork" : 274,
#   The number of sites in the network provided to TRISH by the user.
# "YearScreen" : [1700,1997],
#   Start and end year of mandatory common period of coverage by all chronologies
#   to be used in the reconstruction. Others will be ignored. All must have complete
#   data for all years bracketed. TRISH lets user input these two years and does the 
#   screening using the last two columns of the metadata file. Outside TRISH this input
#   has no effect.
# "NafterYearScreen" : 36,
#   Number of chronologies remaining after screening for "YearScreen." Outside of TRISH
#   this setting has not effect.
# "NafterPolygon" : 36,
#   Number of chronologies remaining after screening for sites being in the specified
#   geographic domain (e.g., polygon drawn at screen 1 of TRISH). Outside of TRISH
#   this is moot.
# "HydroVariable" :  "RO",
#   Code for the hydrclimatic variable represented by the predictand. Must be a member of
#   a recognized set of codes.
# "ClimDatSet" :  "CRU",
#   The source of the climate data used by UNH water balance model to compute the seasonalized
#   predictand.Must be a member of a recognized set of codes.
# "HydroSeason" : [9,12],
#   Ending month (1=jan, 12=dec) and number of months in season of predictand. In TRISH,
#   this tells TRISH how to seasonalized from the original monthly climate data stored. Outside
#   of TRISH, this just defines the season of the input hydrologic time series predictand.
# "yrgoc" : -99999, 
#   Desired start of y data to be used for calibrating SSR models; If -99999, use earliest
#   possible. Here, -99999 is used instead of NA because JSON does not handle NA
#   if NA, the earliest year of flow is used; if yrgoc is incorrectly specified a before the
#   start of flow data, the first year is forced to the start year of flows
# "yrspc": -99999,
#   Desired start of y data to be used for calibrating SSR models; If -99999, use latest
#   possible. Analogous conditions regarding NA, etc., apply as for yrgoc 
#   Start and end year of flow data to be used for calibration of SSR models and final model.
#   If NA, all years of flow overlapping the tree-ring series are used as the calibration period
# "ktran" : 1,
#   Optional transformation of the predictand, y, before reconstruction.
#   None (1), square root(2) or log10 (3). Setting ktran allows the user to 
#   transform y before regression. The resulting reconstruction will  be in units 
#   of transformed flow. TRISH issues error messages if the selected transformation
#   is inconsistent with the data (e.g., log10 transformaton when there are flows 
#   of zero)
# "methMSR" :  2,
#   Method for recontruction. This and PCApredictors effectively specifies the 
#   method to be used. (1) Simple Linear regression (SLR), (2) MLR on SSRs or their PCs, 
#   (3) Analog Nearest neighbor PCA 
# "PCApredictors" : true,
#   Whether PCs fo the SSRs (true) or the SSRs themselves comprise the pool of
#   potential predictors for the multi-site reconstruction.
#  kHowPCA <-2  # if PCA, do it on the correlation (1) or covariance (2) matrix
#   Covariance matrx makes more sense if you attach importance to differences in the 
#   variances of indivual time series on which the PCA is run. Because ReconAnalog runs
#   PCA on single-site reconstructions (SSR), and because the variance of an SSR reflects
#   its accuracy of reconstruction, the best selection here is kHowPCA=2. However,
#   the option is available to do the PCA on the correlation matrix in case you are
#   interested in sensitivity of reconstructions to such a choice.
# "PCoption" : 2,
#   If reconstruction method is MLR of y on PCs of the SSRs (methMSR=2, PCApredictors=TRUE),
#   PCoption specifies how many and which PCs should be in the pool of potential
#   predictors. If a different reconstruction method, this setting has no effect. 
#   Options are 1) directly specify the number of PC, and  2) use the m<fN PCs whose 
#   scores are highest-correlated with y. Here, N is the number of years for calibration
#   of the MSR model, and f is a factor less than or equal to 0.2. By default, f=0.10, 
#   meaning that the pool of potential predictors must be less the 1/10 N, where 
#   N is the length of the calibration period.
# "nPCsKeep" :  1,
#   If PCoption=1, nPCsKeep means to keep the nPCskeep eigenvectors with highest
#   eigenvalues in the pool of potential predictors of the MSR model. If PCoption=2,
#   the setting for nPCsKeep has no effect.. 
# "f" :   0.10,
#   A decimal fraction limiting the sized of the pool of potential predictors of the
#   MSR model to less than this specified decimal fraction of number of years in the 
#   calibration period of the model. For example if 100 years in the calibration period,
#   f=0.10 would mean must have fewer than 10 potential predictors in that pool.
# "alphaR"  :  0.05,
#   Critical alpha-level for correlation (r) of y with SSRs or their PCs. In a TRISH 
#   window, the critical level of r is plotted on a bar chart of correlations. 
#   For the analog (methMSR=3) method, PCs whose r with the predictand, y, is not
#   significantly different than zero at alphaR are not used in identifying analogs. 
#   The only acceptable values for alphaR are 0.10, 0.05, and .01, corresponding to
#   90%, 95% and 99% significance by a two-tailed test of # the null hypothesis of 
#   zero population correlation. If methMSR=1 or methMSR=2, alphR has no effect.
# "Lcausal" : true
#   Whether SSR models should be rejected if the only predictors for y in year t are
#   tree-ring variable from earlier years (e.g., t-1 or t-2). It is sensible to set
#   Lcausal=true if the chronologies are standard chronologies, but may not be such
#   a good idea if chronologies are residual chronologies. 
#   Lcausal=TRUE instructs to reject any final lag model for which the only lags on the 
#   tree-ring series are negatve (implausible to predict current years annual flow, say,
#   from past tree-ring data only). Beware your might want to set this as FALSE if
#   you use residual rather than standard chronologies. Conversion of standard core
#   indices to residual core indices could conceivably severely distort the climate signal
#   in the residual index if that signal itself (climate, not tree-growth) has very strong 
#   autocorrelation If so, a negative lag might be critical for restoring the climate
#   signal removed in tree-ring standardization.

######## LOAD LIBRARIES
 
library(car) # for durbinWatsonTest()
library("rjson") # for reading json input files

################################################################################

# Remove all variables except functions
rm(list = setdiff(ls(), lsf.str()))

################################################################################
#
# JSON input. Input data as list from json; in for loop, store the list elements

naJSON <- -99999 # regard -99999 as NA
myData <- fromJSON(file="init01.json")
X<- myData
for (j in 1:length(X)){
  a <- names(X)[j]
  b <- paste(a,'<-X$',a,sep='')
  s <- paste(a,'<-',b)
  eval(parse(text=s))
}
rm(a,b,s)
if (yrgoc == naJSON){
  yrgoc <- NA
}
if (yrspc == naJSON){
  yrspc <- NA
}


######## SOURCE THE FUNCTIONS ReconAnalog() DEPENDS ON

source(paste(code_dir,"TranFlow.R",sep="")) # optional transformation of flows
source(paste(code_dir,"trimnan.R",sep="")) # get indices to allow trimming of leading and trailing NA
source(paste(code_dir,"emssgUNH.R",sep="")) # write error file to system, specified output folder
source(paste(code_dir,"reconsw4.R",sep=""))  # single-site reconstruction (SSR) by distributed-lag stepwise regres
source(paste(code_dir,"TrimTsm1.R",sep=""))  # trim tree-ring time series matrix in preparation for SSR
source(paste(code_dir,"trimRowNA.R",sep=""))  # row indices of a matrix after trimming off leading and trailing all-NA rows
source(paste(code_dir,"tsmExtend.R",sep="")) # extend time series matrix on recent end
source(paste(code_dir,"RecSLR1.R",sep="")) # Reconstruction by simple linear regression of y on mean of SSRs
source(paste(code_dir,"RecMLR1.R",sep="")) # Reconstruction by multiple linear regression or analog method
source(paste(code_dir,"SignalDrop1.R",sep="")) # drop in maximum accuracy as tree-ring network thins in recent years
source(paste(code_dir,"PrewhitenChrons.R",sep="")) # convert tsm of chronologies to prewhitened chronologies


########## READ FILES OF PREDICTAND, TREE-RING TIME SERIES, TREE-RING METADATA

V <- read.table('hydroData.txt',sep="\t",header=TRUE) # input flow; UNH will not read this in because
# usere picks the hydroData variable interactively in TRISH
U <- read.table('treeData.txt',sep="\t",header=TRUE) # input chronologies
Tmeta <- read.table('treeMeta.txt',sep="\t",header=TRUE) # input chronologies
# Status. U is tsm of chronologies and Tmeta is table of corresponding metadata. So far
# the chronologies in the tsm and metadata table are those satisfying (1) in polygon, and (2) complete


########## OPTIONAL PREWHITENING OF CHRONOLOGIES
#
# Option to prewhiten (remove autocorrelation from) chronologies using autoregressive model
# of order 1, 2, or 3 (AR(1),AR(2), AR(3)). Order 0 instructs to not prewhiten. 
# Prewhitening essentially allows user to check whether "residual" chronologies might have
# stronger signal than standard chronologies for y. If so, user may want to develop
# residual chronologies as separate network and submit those to TRISH for testing. 
# Prewhitened chronologies are not exactly the same as residual chronologies because with
# prehitening the autoregressive modeling is done on the site chronology, but with 
# residual chronologies (as defined in dendro literature) the modeling is done on the 
# standard core indices before averaging those into a site chronology. But experience
# has shown the two versions (prewhitened and residual) are very similar.

if (PrewhitenOrder==0){
  # No action needed
  PWtext <- 'Not prewhitened'
  } else {
    if (any(PrewhitenOrder==c(1,2,3))){
      # Call function to covert U to prewhitened chronologies
      # If order p, will end up converting first p values of chronology to NA (lose p leading years)
      PWtext <- paste('Prewhitened with AR(',as.character(PrewhitenOrder),') model',sep='')
      ResTemp <- PrewhitenChrons(U,PrewhitenOrder,outputDir)
      U <- ResTemp$Xwhitened
    } else {
      # Invalid order (must be 1,2,or 4); stop with error message
      # Write error file
      emssg <- 'Allowable AR(p) prewhitening models are p=0,1,2 or 3'
      ResTemp<-emssgUNH(emssg,outputDir)
      stop(emssg)
    }
  }

# Get first and last year of the tree-ring TSM. It is assumed that this matrix has 
# no all-NA rows
yrTree <- U[,1]
yrgoTree <- yrTree[1]; yrspTree <- yrTree[length(yrTree)]  

##########  CHECK NUMBER OF SERIES IN TREE-RING MATRIX

nms1<-colnames(U[-1]) # chronology ids, for all chronologies from polygon screening
nchron<-length(nms1) # number of chrons from polygon screening
if (NafterYearScreen != nchron){
  # Write error file
  emmsgUNH <- 'Number of data columns of input tree-ring matrix inconsistent with nAfterPolygon'
  ResTemp<-emssgUNH(emssg,outputDir)
  stop(emssg)
}

######### BUILD SOME STRINGS TO BE USED LATER IN LABELS

ClimDataSet <- c('CRU','Delaware','Reanalyis') # input climate data used by UNH WBM
# This variable is included here for information only. When ReconAnalog is used
# within TRISH, the user has a dropdown window that allows selection of the input
# data. But ReconAnalog does not care about the input data set and works only with
# whatever processed WBM output or seasonalized climate series is provided


Dtypes <- c('RO','P','T','Q')
LabsUnits <- c('(mm)','(mm)','(Deg C)','(cms)')
HydNames  <- c('Runoff','Precipitation','Temperature','Discharge')
ithis = which(Dtypes==HydroVariable)
Dtype<-Dtypes[ithis]; LabUnits<-LabsUnits[ithis]; HydName<-HydNames[ithis]


########## RENAME A VARIABLE; SET UP MENU OF METHODS

n1<- NsitesUserNetwork # number of sites in the full user-supplied network; 
# system must provide this; renamed here for convenience
msrMethods <- c('SLR1','MLR1-PCA','MLR1-noPCA','Analog')


########## BUILD FILENAME FOR PDF DESCRIBING MSR METHOD TO BE USED
#
# This pdf will we written to system output folder. User and refer to the pdf for
# details on the method used, including definitions of downloaded output files

if (methMSR==1){
  jMethod <- 1
} else if (methMSR==3) {
  jMethod <-4
} else {
  if (isTRUE(PCApredictors)) {
    jMethod <-2
  } else {
    jMethod <- 3
  }
}
PdfFile <- paste('TrishOutputDescribe',msrMethods[jMethod],'.pdf',sep='')
PdfDescribe <- paste('See',PdfFile,'for more',sep=' ')
rm(jMethod)

########## CLEAR OUTPUT DIRECTORY AND COPY PDF DESCRIBING MSR METHOD THERE

PathClean <-  paste(outputDir,'*',sep='')
unlink(PathClean,recursive=F,force=F)
file.copy(from=paste(pdf_dir,PdfFile,sep=''),to=paste(outputDir))
rm(PathClean)

############# HARD CODED SETTINGS
# 
# These setting not to be changeable by user. Can be used by developer to explore 
# possibilies for extending or modifying ReconAnalog. User should not be able to change
# these setting from TRISH menus.

nNeg<-2 #$ maximum negative lag allowed on chronologies for SSR models
nPos<-2 # maximum positive lag allowed on chronologies for SSR models
# ReconAnalog() was written to specifically apply lagged regression with maximum of
# 2 negative and positive lags on the chronology. Accordingly, leaving 9 out in 
# cross-validation guarantees that cross-validation estimated do not depend on
# any of the tree-ring data that are used also in calibrating the cross-validation
# model. 
# yrgo1 and yrsp1 best both set to NA. This lets the time coverage of the tree-ring data itself
# determine the time comverage of the reconstruction. In some early trialsl I played with varying
# yrgo1 and yrsp1, but decided makes more sense to go with NA
yrgo1<-yrgoTree+nNeg # desired start year of reconstruction; actual reconstruction coverage will depend on 
#   coverage of tree-ring chronologies in final model
yrsp1<-yrspTree-nPos  # desired end year of reconstuction (including through calibration period)
# yrgo1 and yrsp1 are the desired start and end years of reconstructed flow.
# The tree-ring matrix supplied by TRISH will be trimmed in ReconAnalog() to include only
# those chronologies with data in year yrgo1-2 (allows for two negative lags in 
# single-site regression).
# This will eliminate from consideration any chronologies that start at a later year. 
# The tree-ring matrix will be trimmed to end in the earlier of (yrsp1+2) and (the last year
# of data for any site passing the screening for yrgo1). If yrgo1=NA, the matrix is trimmed to
# start with the first year for which all sites in the basin domain have data. If yrsp1=NA,
# the matrix is trimmed to end with the last year of data at any one of the sites.
N1 <- 50 # in forward extension of SSR matrix, common period of all SSRs must be
# at least N1 year (e.g., 50
N2 <- 100 # in forward extension of SSR matrix, series needing an extension in year
# i must overlap by N>=N2 years with some other series that does have a value in 
# year i.
N3 <- 30 # minimum acceptable number of years for calibration on MSR. Allows 15/15, which is
# ridiculously low, for split-sample validation. 
incR2a<-0.01 # Critical increment in adjusted R-squared of MSR model. Stepwise model is assumed to reach
# "approximate maximum" when next step would yield increase in adjusted R-squared less than inR2a.
# Stepwise models in SSRs and MSR are not allowed to proceed beyond the approximate maximum of
# adjusted R-squared. Further, depending on "kstop" (see nex), the model may stop an an even
# earlier step to satisfy constraints on cross-validation. 
kstop <-2 # Stepwise forced to stop at either the maximum adjusted R-squared (kstop=1) or at some earlier
# step if cross-validation RE reaches a maximum at an earlier step (kstop=2)
ScreenAnalogPCs <- TRUE # Screen the PCs used in analog reconstruction (methMSr=3)
# by correlation with predictand.If TRUE, only those PCs whose correlations with y
# are significant at alpha-level alphaR (see earlier) are used to identify analogs
#
#
######### MAKE FLOWS MATRIX; TRIM OFF ANY LEADING OR TRAILING NA; STOP IF INTERNAL NA

V<-as.matrix(V)
v<-as.matrix(V[,2])
yrv <- V[,1,drop=FALSE]
i1<-trimnan(v) # row indices of v with leading and trailing NAs stripped
d = diff(i1)
L<-!all(d==1)
if (L) {
  stop('Internal missing value in V')
  #? ??UNH?? This stop() message will be changed so that instead of stop(), a
  # new function, emssgUNH() is called. The new function writes and error message
  # to the UNH unix system. My general approach will be NOT to use stop() 
  # messages such as above, but to instead write out a specific error message
  # to system when aborting.
  # ??UNH??end
}
V<-V[i1,,drop=FALSE] # trim the matrix of leading and trailing NAs
v <- v[i1,,drop=FALSE]; yrv<-yrv[i1,,drop=FALSE]
rm(L,d,i1)

############### TRANSFORM FLOWS (OPTIONAL)
# Allowed are square root or log10. Call a function to do the transform. If transform
# not reasonable physically, function called returns a flag that prompts this script
# to abort and also prints a message to outputDir

kBogus<-FALSE
Transformed<-FALSE
if (ktran==1){
  # If not call TranFlow, set ResTf to empty list
  ResTf<-vector(mode = "list", length = 0) 
} else if (ktran==2 || ktran==3) {
  ResTf <- TranFlow(v,ktran)
}  else {
  kbogus<-TRUE
}
#--- UNH HANDLING OF TRANSFORM MESSAGE

emssg<-'None' # initialize error message
if (kBogus){
  emssg<-'Invalid specified ktran: option must be 1, 2 or 3'
  # ??UNH?? this is the call to the new function mentioned previously ??UNHend??
}
if (length(ResTf)==0){
  sTran<-''
} else {
  if (ResTf$flag==0){
    sTran<-ResTf$sTran
    V[,2]<-ResTf$x
    
  } else if (ResTf$flag==1) {
    sTran<-''
    emssg<-'Sqrt transform invalid;<br />series has negative values'
  } else if (ResTf$flag==2) {
    emssg<- 'Log10 transform invalid;<br />series not all-positive'
  }
}

# Conditional error message to OutputDir
if (emssg=='None'){
} else {
  # Write error file
  ResTemp<-emssgUNH(emssg,outputDir)
  stop(emssg)
}

# Depending on transform, modify y units, and store with other information applicable to 
# all four MSR methods. First might have to modify label for units of y
if (ktran==1){
  # no transform, units OK
} else if (ktran==2) {
  LabUnits<-paste('[sqrt',LabUnits,']')
} else if (ktran==3){
  LabUnits<-paste('[log10',LabUnits,']')
}
txtSeas <- paste0(as.character(HydroSeason[2]),'-month season ending in month ',as.character(HydroSeason[1]))
RecListx<-c(Dtype,LabUnits,HydName,txtSeas) # general list for calls to MSR recon


###### CHECK FOR INTERNAL NA IN TREE-RING SERIES, ABORT WITH MESSAGE IF FOUND

emssg<-'None'
for (n in 1:nchron){
  j<-n+1
  u <- as.matrix(U[,j])
  ResTemp<-trimnan(u) # index to rows of nan-trimmed u
  u1<-u[ResTemp]
  L<-any(is.na(u1))
  if (L){
    emssg<-paste('Internal NA in tree-rng series ',as.character(n),': ',nms1[n])
    ResTemp<-emssgUNH(emssg,outputDir)
    stop(emssg)
  }
}

################################################################################
#
# TRIM TREE-RING MATRIX TO COVER JUST THE PERIOD NEEDED FOR SSR'S; COL-TRIM AS
# NEEDED SO THAT TSM INCLUDES ONLY THE CHRONS THAT COULD PROVIDE THE RECONSTRUCTION
# TIME COVERAGE; ALSO ROW-TRIM THE TREE-RING METADATA TABLE ACCORDINGLY
#
# See earlier comments for inputs yrgo1 and yrsp1
#
# The reconstruction method allows lags t-2 to t+2 on tree ring. Therefore, if you
# specify you want the reconstruction to start in year t0, all series must have 
# data in year t0-2. This next section may therefore reduce the number of colunns in 
# the tree ring matrix by removing those without data in year to-2. The matrix will
# be row-truncated to begin in year yrgo1-2 and to stop in the earlier of
# [year yrsp1+2 and the last year with data for any chronology]

mlead<-2; # this many leading tree-ring values are needed to produce first reconstructed value
ResTrim1 <- TrimTsm1(U,yrgo1,yrsp1,nNeg,nPos)
X<-ResTrim1$X # tsm of tree-ring chronologies to be converted by SSR
ix2 <- ResTrim1$ix # index of series in X to the columns of tsm of map-screened chronologies
# provided by system
nms2 <- nms1[ix2] # column names of the (possibly) column-reduced tree-ring matrix
# These ids are identically equal to the colnames (after Year) of X, as retured
# by TrimTsm1()

# Metadata trimming,and check for exact match of column headers of chronology tsm with
# Ids in the metadata
Tmeta <- Tmeta[ix2,]  # row-trim
IdMeta=Tmeta$Id
IdMeta <- gsub(" ", "", IdMeta, fixed = TRUE)
L<-identical(IdMeta,nms2)
#--- Bomb message if not exact match
emsg1 <- c('Ids in header row of time series matrix of chronologies do not exactly match Ids in metadata')
if (!L){
  eBomb<-emssgUNH(emsg1,outputDir)
  stop(eBomb)
}


################################################################################
#
#  CONVERT TREE-RING SERIES TO SINGLE-SITE RECONSTRUCTIONS (SSR's) OF FLOW
#
# All series in X will be converted. The regression statistics and the SSRs will be
# stored. Only a subset of those SSRs, dependkng on the calibration and validation
# statistics, will ultimately be used in the PCA analysis and final reconstruction.
# The SSR statistics will include a "reject" flag, which will be used to screen out
# bad tree-ring sites. A site will be rejected if any of the following are true:
# 1) p-value of overall F of regression model >=0.05
# 2) REcv (Reduction of error statistic from leave-9-out cross-validation) <=0
# 3) RE from either half of split-sample cross-validation <=0
#
# SSR modeling is done by reconsw4(), and comments there have more details on the method.
# Essentially, y(t) (flow in year t) is regressed stepwise forward on tree-ring
# {x(t-2), x(t-1), x(t), x(t+1),x(t+2)}. Entry of variables is stopped when an
# an additional step would lead to less than a c1 increment in adjusted R squared. The
# final model is truncated at the step beyond which RE of cross-validation does
# not ncrease by at least c1. 

# SSR time series will be in matrix Y1; these will be the full set of estimated SSRs
nchron2 <- dim(X)[2]-1 # number of chronologies after screening for time coverage
nY1<-nchron2;
mY1 <- dim(X)[1] # Y1 will store the SSRs, and will initially be same size as X
Y1 <- matrix(data=NA,nrow=mY1,ncol=nchron2) # to hold the SSRs
yrY1 <- matrix(X[,1])
yrgoY1<-yrY1[1,1]
yrspY1<-yrY1[mY1,1]


#--- Table of statistics of all models passing the spatial-temporal screening.
# The data frame SSRdf1 is the summary table, described below. A second table,
# SSRdf2, has a similar structure, but includes only those site passing statistical
# screening for a strong, validated, temporally stable signal for flow.

SSRdf1 <- data.frame(matrix(nrow = nchron2, ncol = 15))
names(SSRdf1) <- c("N1", "N2","Site","StartC","EndC","Model","R2adj","Fp",
"REcv","REa","REb","Refit","StartR","EndR","Reject")
# N1, N2, are sequential number in this table, and in the user's full network.
# Site is a site id; StartC & EndC are first and last years of calibration period;
#   LagsIn is a coded string indicated which of lags t-2 to t+2 are in final model;
# R2adj is regression adjusted R-square; Fp is the p-value of overall F of the
#   regression (pF<0.05 means significant model); 
# REcv is th reduction-of-error statistic from leave-9-out cross-validation; REa 
#   and REb are reduction-of-error for split-sample calibration/validation (e.g., 
#   REa is for validation on first half when model is fit to second half; 
# Refit is a logical (0 0r 1) indicating whether the final selectio of lags
#   allowed refitting of the model to a slightly longer calibration period than 
#   possible if all lags t-2 to t+2 are in model; 
# StartR and EndR are the start and end years of the reconstruction (SSR for this site)
# Reject is a logical (0 or 1) indicating if site was rejected for further use in 
#   in reconstruction. Note that only sites not rejected are used in the later
#   PCA step. Numbering n3 therefore goes with sites having a 0 in this colunn.
#   ) on last, validation on first half of 

#--- Generate SSRs
#
# Function reconsw4() is called in a loop over all chronologies in the drawn polygon
# that cover the minimum acceptable reconstruction interval input at the TRISH window.
# Each chronology is converted to an individual estimate of flow by lagged stepwise
# forward regression. Five lagged values, lagged t-2 to t+2 from the year of flow, are
# considered as potential predictors. The forward stepwise process is stopped when an
# an additional step would result in less than c1 increase in adjusted R squared.
# Leave-9-out cross-validation and split-sample calibration/validation are then applied
# to test the skill of prediction and to possibly simplify the model.
# 
# Status: V and X are data frames with prepared flows and tree-ring chronologies

timeGo<-proc.time() # start timer to check cpu time and clock time for running 

yrsCalWindow <- c(yrgoc,yrspc) # calibration of model will consider flows only within this window
#nchron<-22
SSRprelim<-vector(mode="list",nchron2) # to hold lists from preliminary SSR modeling for each chronology
SSRlags1 <- vector(); SSRlags2<-vector() # initialize empty vectors to hold concatenated lags
# in models (all SSRs for1 ; those passing final screening for 2)

# Fill list SSRdf1 members with network site numbers and site IDs as supplied in list treeMeta
SSRdf1[,2]<-Tmeta$N2
SSRdf1[,3]<-Tmeta$Id

for (n in 1:nchron2){ # loop over tree-ring chronologies
  x <- as.matrix(X[,c(1,(n+1))]) # 2-col data frame with year and chronology
  ResSSR <- reconsw4(x,V,nNeg,nPos,yrsCalWindow,incR2a,Lcausal,LallowLags)
  SSRprelim[[n]]<-ResSSR
  SSRdf1[n,1]<-n
  SSRdf1[n,4]<-ResSSR$yearsCal[1]
  SSRdf1[n,5]<-ResSSR$yearsCal[2]
  SSRdf1[n,6]<-ResSSR$ModelCoded
  SSRdf1[n,7]<-ResSSR$RsquaredAdj
  SSRdf1[n,8]<-ResSSR$Fp
  SSRdf1[n,9]<-ResSSR$REcv
  SSRdf1[n,10]<-ResSSR$REa
  SSRdf1[n,11]<-ResSSR$REb
  SSRdf1[n,12]<-ResSSR$Lrefit
  SSRdf1[n,13]<-ResSSR$yearsRec[1]
  SSRdf1[n,14]<-ResSSR$yearsRec[2]
  SSRdf1[n,15]<-ResSSR$Lreject
  
  # Build vector of lags in models (1=t-2... 5=t+2)
  # Make two versions: one for all SSRs, the other for those passing screeing.
  # Result (SSRlags1, SSRlags2) are vectors with the concatenated lags in models.
  # Values 1,2,3,4,5 correspond to lags t-2,t-2,t,t+1,t+2
  
  SSRlags1 <- c(SSRlags1,ResSSR$Model)  
  if (!ResSSR$Lreject){ 
    SSRlags2 <- c(SSRlags2,ResSSR$Model)
  }

  #---Store this SSR in matrix of SSRs
  irow <- ResSSR$yhat[,1]-yrgoY1[1]+1 # target rows in Y1
  Y1[irow,n]=ResSSR$yhat[,2,drop=FALSE]
}


#---Trim off any all-NA rows of SSR matrix Y1
i1<- trimRowNA(Y1)
if (is.vector(i1)){
 Y1<-Y1[i1,,drop=FALSE]
 yrY1<- yrY1[i1,,drop=FALSE]
}
rm(i1)

mSSRdf1<-nrow(SSRdf1) # number of rows in the table



#--- Write SSR results for all chronologies models to a file using fprintf
pf1<-file.path(outputDir,paste("Table1-SSR1",".txt",sep=""))
if (file.exists(pf1)){file.remove(pf1)} # must remove old versio of file
fmt1<-"%4s%4s%8s%7s%5s%8s%5s%9s%8s%6s%6s%7s%6s%6s%7s\n"
fmt2<-"%4d%4d%8s%7d%5d%8s%5.2f%9.2G%8.2f%6.2f%6.2f%7s%6d%6d%7s\n"

TitleLine <- 'Table1-SSR1 - Statistics of single site reconstruction (SSR) models'
fprintf('%s\n\n',TitleLine,file=pf1,append=FALSE)
fprintf(fmt1,"N1","N2","Site","Goc","Endc", "Model", "R2a",  "pF", "REcv", "REa",
        "REb", "Refit", "Gor", "Endr", "Reject",file=pf1,append=TRUE)
for (n in 1:mSSRdf1){
  fprintf(fmt2,SSRdf1[n,1],SSRdf1[n,2],SSRdf1[n,3],SSRdf1[n,4],SSRdf1[n,5],
          SSRdf1[n,6],SSRdf1[n,7],SSRdf1[n,8],SSRdf1[n,9],SSRdf1[n,10],
          SSRdf1[n,11],SSRdf1[n,12],SSRdf1[n,13],SSRdf1[n,14],SSRdf1[n,15],
          file=pf1,append=TRUE)
}
fprintf('%s\n\n','',file=pf1,append=TRUE)
fprintf('%s\n','This table applies to chronologies before screening for hydrologic signal',
        file=pf1,append=TRUE)
fprintf('%s\n',paste('Chronologies:',PWtext),file=pf1,append=TRUE)
fprintf('%s\n',PdfDescribe,file=pf1,append=TRUE)


#--- Make second table, SSRdf2, a data frame that is a subset of SSRdf1 with just those
# chronologies not rejected according to the Lreject criterion.
Lreject<-SSRdf1[,15]

# If no chronologies have stable signal for flow, punt
if (all(Lreject)){
  emssgThis<- paste('ReconAnalog aborted: no chronology has a stable sigmal for', HydName)
  eBomb<-emssgUNH(emssgThis,outputDir)
  stop(eBomb)
}

# If only one chronology has as stable signal and you have called for PCA, punt
nGood <- sum(!Lreject,na.rm=TRUE) 
if ((nGood == 1) && (methMSR==2 | methMSR==3) ){
  emssgThis<- paste('ReconAnalog aborted: Cannot use multivariate methods for reconstruction; only one chronology has a stable sigmal for', HydName)
  eBomb<-emssgUNH(emssgThis,outputDir)
  rm (nGood)
  stop(eBomb)
}


ix3 <- ix2[!Lreject] # col-pointer of "accepted" series to original tree-ring matrix
nms3 <- nms2[!Lreject] # ids of series passing screening
SSRdf2<-SSRdf1[!Lreject,]
mSSRdf2<-nrow(SSRdf2) # number of rows in the table
j1<- 1:mSSRdf2
SSRdf2[,1]=j1

pf2<-file.path(outputDir,paste("Table2-SSR2",".txt",sep=""))
if (file.exists(pf2)){file.remove(pf2)} # must remove old version of file
fmt1<-"%4s%4s%8s%7s%5s%8s%5s%9s%8s%6s%6s%7s%6s%6s%7s\n"
fmt2<-"%4d%4d%8s%7d%5d%8s%5.2f%9.2G%8.2f%6.2f%6.2f%7s%6d%6d%7s\n"

TitleLine <- 'Table1-SSR2 - Statistics of screened single site reconstruction (SSR) models'
fprintf('%s\n\n',TitleLine,file=pf2,append=FALSE)

fprintf(fmt1,"N1","N2","Site","Goc","Endc", "Model", "R2a",  "pF", "REcv", "REa",
        "REb", "Refit", "Gor", "Endr", "Reject",file=pf2,append=TRUE)
for (n in 1:mSSRdf2){
  fprintf(fmt2,SSRdf2[n,1],SSRdf2[n,2],SSRdf2[n,3],SSRdf2[n,4],SSRdf2[n,5],
          SSRdf2[n,6],SSRdf2[n,7],SSRdf2[n,8],SSRdf2[n,9],SSRdf2[n,10],
          SSRdf2[n,11],SSRdf2[n,12],SSRdf2[n,13],SSRdf2[n,14],SSRdf2[n,15],
          file=pf2,append=TRUE)
}
fprintf('%s\n\n','',file=pf2,append=TRUE)
fprintf('%s\n','This table applies to chronologies passing screening for hydrologic signal',
        file=pf2,append=TRUE)
fprintf('%s\n',paste('Chronologies:',PWtext),file=pf2,append=TRUE)
fprintf('%s\n',PdfDescribe,file=pf2,append=TRUE)

#--- Make time series matrix of SSRs passing signal for screening (n0n-rejects)
Y2 <- Y1[,!Lreject,drop=FALSE]
yrY2 <-yrY1

#---Trim off any all-NA rows of SSR matrix Y2
i1<- trimRowNA(Y2)
if (is.vector(i1)){
  Y2<-Y2[i1,,drop=FALSE]
  yrY2<- yrY2[i1,,drop=FALSE]
}
mY2<-dim(Y2)[1] # number of rows in Y2
jScreened=SSRdf2$N2 # pointer from cols of screened SSRs
# to columns in original users tree-ring network

Tcpu<- (proc.time()-timeGo)[1] #...... processing time
Tclock<- (proc.time()-timeGo)[3] #...... clock time


################################################################################
#
# ANALOG EXTENSION OF SCREENED SSR'S ON RECENT END
#
# Have, say, N SSRs. Have a target end year that you want all series to cover.
# Loop over all N SSRs, each time defining the current series as "key" and the
# rest as "others." Loop over key series: (1) Spearman r of key with others, and
# sorting of others from most correlated to least. (2) Proceed for next steps going
# from most to least correlated. (3) Pull full overlap of the two series. 
# (4) Loop over the years of key needing filled in. (5) Fill in all that are possible 
# from this member of others. (6) If still values to fill, proceed over more
# of the others. 
#
# Analog method used, assuming have a key series and a predictor is (1) find
# full overlap over two series. Sort the two series from smallest to largest for
# that overlap. Have the value, x, of predictor series for year with data missing
# at key. Compute the non-exceedance probability of that value in the predictor 
# serie for the overlap. Assign the estimate as the interpolated value of the 
# predictor series with the same non-exceedance probability in the overlap. 

# Truncate matrix Y2, yrY2 on early end so that first row has no NA. After truncation,
# Y2 should have no NAs on early end, but generally will have NAs in some cols on recent end
i1 <- which(complete.cases(Y2))
Y2 <- Y2[i1[1]:mY2,,drop=FALSE]
yrY2 <- yrY2[i1[1]:mY2,,drop=FALSE]

#--- Call function to extend tree-ring matrix on recent end

ResME <- tsmExtend(Y2,yrY2,yrsp1,N1,N2) # returns names list with Y, yrY

#--- Bomb out messages from tsmExtend()
emsgs2 <- c('No need to extend','OK, but yrsp later than last year of data in X',
            'No problem','tsmExtend aborted: common period of all series too short',
            'tsmExtend aborted: insufficient common period of predictor and predictan',
            'tsmExtend aborted: yrX not continuous','tsmExtend aborted: first year of X has some NA')
khow<-ResME$khow
if (khow>3){
  emsg2<-emsgs2[[khow]]
  eBomb<-emssgUNH(emsg2,outputDir)
  stop(eBomb)
}

Y3 <- ResME$Y; yrY3 <- ResME$yrY # forward-extend tree-ring matrix
if (any(is.na(Y3))){
  eBomb<-emssgUNH('ReconAnalog() aborted: matrix Y3 has a NA',outputDir)
  stop(eBomb)
} 


################################################################################
#
# FIGURES FOR SINGLE-SITE RECONSTRUCTIONS
#
# Idea is that in TRISH user can use radio buttonto choose from 1) bar chart of site-screeing,
# 2) boxplot summaries of model statistics, 3) line plot of time coverage of SSRs, and
# 4) z-score time-series plot (annual and smoothed on one set of axes; mean of the SSRs 
# after converting to z-scores)
n2 <-nchron # number of sites returned by polygon screening
n3 <- length(ix2) # number of chronologies after screening for reconstruction window
n4 <- dim(Y3)[2] # number of chronnologies passing final screening for hydrologic signal



#----FIGURE 01.  1x2 OF 1)BAR PLOT SUMMARIZING NUMBER OF CHRONOLOGIES, ORIGINAL AND SCREENED,
# AND 2) ADJUSTED R-SQUARED OF ALL FITTED SSR MODELS AND OF THOSE PASSING FINAL
# SCREENING FOR HYDROLOGIC SIGNAL


#--- Numbers of chronologies (bar chart)

# xtemp <- c(n1,n2,n3,n4)
# DeltaTemp<- 0.1*(max(xtemp)-min(xtemp))
# yhi <- DeltaTemp + max(xtemp)
# ylo <- 0
# ylims <- c(ylo,yhi)

# Preliminaries
xtemp <- c(n1,n2,n3,n4)
DeltaTemp<- 0.1*(max(xtemp)-min(xtemp))
DeltaTemp2 <- DeltaTemp/1.5
#yhi <- 0.1*(max(xtemp)-min(xtemp)) + max(xtemp)
yhi <- DeltaTemp + max(xtemp)
ylo <- 0
ylims <- c(ylo,yhi)
barnames <- c('N1','N2','N3','N4')

strAnnote1<-c('N1: source network','N2: polygon-selected',
              'N3: SSR-modeled','N4: Passed screening')

png(filename=paste(outputDir,"Figure01-SSR1.png",sep=""), width = 960, height = 480)
layout.matrix <- matrix(c(1,2), nrow = 1, ncol = 2)
layout(layout.matrix,heights=2,widths=c(2,1))
par(mar=c(5,4,4,1),cex.main=1.3)

bp <- barplot(xtemp,ylim=ylims,xlab='Screening Stage',col='Pink',border=TRUE,
              names.arg=barnames,main='Number of Chronologies',cex.lab=1.3)
text(bp,xtemp+DeltaTemp/2,labels=xtemp)
abline(h=0)

# Annotate meaning of N1, N2, N3, N4
xtemp=3
for (n in 1:4){
  if (n==1){
    ytemp = ylims[2]
  } else {
    ytemp <- ytemp-DeltaTemp2
  }
  text(xtemp,ytemp,strAnnote1[n],adj=c(0,1),cex=1.3)
}

rm(xtemp,barnames,ylims,bp,DeltaTemp,DeltaTemp2)

#--- Adjusted R squared

#screen(2)
par(mar=c(5,4,4,1),cex.main=1.3)

namesBP<-c('N3','N4')
boxplot(SSRdf1[,7],SSRdf2[,7],notch=FALSE,
        ylab = "Adj R-squared of SSR Models",
        main="Adjusted R-squared",names=namesBP,
        cex.lab=1.2)
dev.off()


#----FIGURE 02.  1x2 OF HISTOGRAMS OF WHICH LAGS ARE IN THE SSR MODELS. 1) ALL SSRS, AND
# 2) THOSE SSRS PASSING SCREENING FOR HYDROLOGIC SIGNAL. EACH SSR MODEL MAY HAVE
# FROM 1-5 LAGS, RANGING FROM t-2 to t+2 YEARS FROM THE YEAR OF FLOW. THE HISTOGRAMS
# SUM OVER MODELS, SO THA THE GRAND TOTAL NUMBER OF LAGS IN THE HISTOGRAM IS
# GREATER THAN THE NUMBER OF MODELS.


#--- Left: Histogram of lags, all N3 SSRs 

png(filename=paste(outputDir,"Figure02-SSR2.png",sep=""), width = 960, height = 480)
layout.matrix <- matrix(c(1,2), nrow = 1, ncol = 2)
layout(layout.matrix,heights=2,widths=c(1,1))

# Left histpgram
par(mar=c(5,4,4,1),cex.main=1.4)
hBreaks<-c(0.5,1.5,2.5,3.5,4.5,5.5)
# Changing x axis
xTicks<-seq(1, 5, by=1)
xTickLabs <- c('-2','-1','0','+1','+2')
n3a<-length(SSRlags1)
title1<-paste('Histogram of Lags (',n3,'N3 Models,',n3a,'Total Lags)')
hist(SSRlags1,breaks=hBreaks,xlim=c(0.5,5.5),xaxt='n',
     main=title1,xlab='',cex.lab=1.2)
vtemp<-par('usr')
mtext('Lag (yr)',side=1,line=1.5,cex=1.2)
for (k in 1:5){
  mtext(xTickLabs[k], side = 1, line = 0, outer = FALSE, 
      at = k, adj = NA, padj = NA, cex = NA, col = NA, 
      font = NA)
}
#text(xTicks,0,xTickLabs,adj=c(0,0))


#--- Right: Histogram of lags, SSRs passing screening for hydro signal

ylims<-c(vtemp[3],vtemp[4]) # same y limits as previuous
par(mar=c(5,4,4,1))
hBreaks<-c(0.5,1.5,2.5,3.5,4.5,5.5)
# Changing x axis
xTicks<-seq(1, 5, by=1)
xTickLabs <- c('-2','-1','0','+1','+2')
n4a<-length(SSRlags2)
title2<-paste('Histogram of Lags (',n4,'N4 Models,',n4a,'Total Lags)')
hist(SSRlags2,breaks=hBreaks,xlim=c(0.5,5.5),ylim=ylims,xaxt='n',
     main=title2,xlab='',yaxs='i',cex.lab=1.2)
mtext('Lag (yr)',side=1,line=1.5,cex=1.2)
for (k in 1:5){
  mtext(xTickLabs[k], side = 1, line = 0, outer = FALSE, 
        at = k, adj = NA, padj = NA, cex = NA, col = NA, 
        font = NA)
}
dev.off()

#----FIGURE 03. 1x1. TWO-Y-AXIS TIME PLOT TO HELP GUIDE USER IN CHOICE OF 
# CALIRATION PERIOD FOR MULTI-SITE-RECONSTRUCTION (MSR) MODEL. 
# 
# Before extension by tsmEndtend(), the signal-screeed SSRs generally end in
# different years. These time plots cover the tail years of the SSRs: from 
# the last year with data for all SSRs through the last year of the SSR with
# most recent data. One plot is the maximum adjusted R-squared of available
# SSRs. The other plot is the number of available SSRs. Annotated on the plot
# is also the ending year of the observed hydro variable. The user will eventually
# need to select the end year of calibration of the MSR model. This year
# cannot be later than the last year of the observed hydro series. The year could
# be as late as the last year SSR at any site, but this may not be a good
# idea if the adjusted R-squared of the most recent SSR is small. 


#--- Prepare data needed for the plots
#
# Status. Have full-length observed hydro in 1-col matrices, yrv, v. 
# Have the adjusted R-squared of signal-screened models in SSRdf2[,7].
# Have the corresponding end years of the SSRs in  SSRdf2$EndR. 
# Have time series of SSRs in matrices Y3,yrY3

ResSD <- SignalDrop1(SSRdf2$EndR,SSRdf2[,7])
x1 <- ResSD$x1;   x2 <- ResSD$x2;  yrx1 <- ResSD$yrx1; yrx2 <- ResSD$yrx1; 

# Want x axis to start year before and end year after the relevant period
yrsEnd <- unique(SSRdf2$EndR) # unique ending years of screened SSRs
xHead <- max(yrv) # head of arrow here
yrLo <- min(yrsEnd)
yrHi <- (xHead - yrLo) + xHead
xlims <- c(yrLo-0.05,yrHi+0.05) # limits for x axis
ylims <- c(min(x1)-0.2,max(x1)+0.2) # limits for y axis

# Arrow
yHead <- ylims[1]+ diff(ylims)/2
yTail <- yHead + diff(ylims)/10
xTail <- xHead+ (xlims[2]-xHead)/2 # x position of tail

#--- yyplot
png(filename=paste(outputDir,"Figure03-SSR3.png",sep=""), width = 960, height = 480)
par(mar=c(5,5,4,5),cex.main=1.4)
plot(yrx1,x1,xaxt='n',yaxt='n',type="b",pch=1,col='red3',xlim=xlims,ylim=ylims,cex=1,
     main='Drop in Signal Strength with Loss of SSRs (Chronologies) on Recent End',xlab='Year',
     ylab='Number of SSRs',cex.lab=1.2)
# vertical line at last year of hydro series
axis(1,at=seq(yrLo,yrHi,1))
axis(2,at=seq(min(x1),max(x1),1))
abline(v=max(yrv),col='Magenta',lty=2)

# arrow to the vertical line
arrows(xTail,yTail,xHead,yHead)
txtTemp<- paste('Last year of',HydName)
text(xTail,yTail,txtTemp,adj=c(0,0.5),cex=1.4)
  
par(new = T)
plot(yrx2,x2, pch=16, type='b',lty=0,axes=F, xlab=NA, ylab=NA, xlim=xlims,
     cex.lab=1.3,cex.axis=1.2)
axis(side = 4)
mtext(side = 4, line = 3, 'Maximum adjusted R-squared')
legend("topright",
       legend=c("N of SSRs","Max Adj R-squared"),
       lty=c(1,0), pch=c(1, 16), col=c("red3", "black"),cex=1.3)
dev.off()
rm(txtTemp,xTail,xHead,yTail,yHead,yrsEnd,yrLo,yrHi,xlims,ylims,x1,x2,yrx1,yrx2)


#----FIGURE 04. 1x1. SCATTERPLOT OF OBSERVED HYDRO (SEASONALIZED HYDROLOGIC SERIES) ON THE
# MEAN OF THE SIGNAL-SCREENED SSRS. THIS PLOT WILL HELP USER DECIDE WHAT TO CHOOSE AS METHOD
# FOR THE MULTI-SITE RECONSTRUCTION (MSR)
#
# Status. 
# Y3, yrY3 are matrices with the screened SSRs
# V is matrix of hydro series, with year as col 1 and data as col 2

# Prepare the two series for the scatterplot
w <-rowMeans(Y3, na.rm=TRUE) # average of the SSRs (vector)
yrw = yrY3 # 
W<- as.matrix(cbind(yrw,w)) # bind yrw and w into a time series matrix
ResPC<- PeriodCommon(W,V)# get common period of W and V
yrgo1 <- ResPC$tgo; yrsp1 <- ResPC$tsp # start and end years of common period 
w1 <- ResPC$X[,2]; yrw1 <- ResPC$X[,1]; # mean-SSR and its year, as vectors, for common period 
# with observed hydro
v1 <- ResPC$Y[,2]; yrv1 <- ResPC$Y[,1]; # observed hydro and its year, as vectors, for common period 

# Pearson correlation of observed hydro with mean of SSRs
r = cor(v1,w1)
rStr<- paste('r=',toString(round(r,digits=2)))
strMain=paste('Scatter of Mean of',as.character(n4),
              'Single-Site Reconstructions (SSRs) of',Dtype,'on Observed',Dtype,
              '\n(Fits are straight-line (red) and loess (blue))') # plot title
rm(r)

# Strings for plots
ylabTemp <- paste('Mean SSR',LabUnits)


#--- scatterplot
png(filename=paste(outputDir,"Figure04-SSR4.png",sep=""), width = 760, height = 480)
par(mar=c(5,5,6,4))
# Call function from package car for the scatterplot. This function as called will regress w1 on
# v1 and plot the lease-squared-fit straight line. Also plotted is a loess (local regression) curve
# and loess curves to represent variand of the loess estimat. The loess curves use a a span of 2/3 and
# are estimated by R function loess.R. For the error bars, the negative and positive departures from
# the loess estimate of the mean are squared and themselves fit wit a loess curve. The plotted lines are
# at the square root of those squared-departure fits. Because the two bordering smoothed line represent the
# typical positive and negative departure, the confidence interval can be heuristically considered 
# a 50% confidence interval.
scatterplot(v1,w1,boxplots=FALSE,
            regLine=list(method=lm, lty=1, lwd=2, col="red"),
            ylab=ylabTemp,
            xlab=paste('Observed',Dtype,LabUnits),
            main=strMain,cex.lab=1.2)
text(min(v1),max(w1),rStr,adj=c(0,1),cex=2)
dev.off()
rm (ylabTemp)


################################################################################
#
# COMBINE SSR'S INTO A FINAL MULTI-SITE RECONSTRUCTION (MSR
# Method depends on settings for methMSR and PCApredictors. If PCA is involved in
# reconstruction, method might also depend on settings of nkHowPCA, PCoption 
# and nPCsKeep
#
# methMSR has three possible values: (1) simple linear regression, 
# (2) stepwise multiple linear regression on SSRs (PCApredictors=false) or their 
# PCs (PCApredictors=true), and (3) PCA analog
#
# Simple linear regrssion is done by calling function RecSLR1
# The other methods ared done by calling function RecMLR1

ReconMethods <- c("Simple linear regression of y on mean of SSRs", 
                  "Multiple linear regression of y on SSRs or their PCs",
                  "Principal component analog nearest neighbor")
ReconMethod <- ReconMethods[methMSR]
NextFigNumber<-5 # because SSR has already produced figure files Figure01.png to Figure04.png

# Set calibraton period of MSR to start with latest of [yrgoc;d first available year of hydro series; first available year of mean SSR[
# Set period to end with earliest of [yrgoc; last available year of hydro series' last available year of mean SSR]
if (methMSR==1){
  # Recon by simple linear regression, using RecSLR1(). 
  SpecListSLR1 <- list("PdfDescribe"=PdfDescribe,"Text"=RecListx,"u"=w,"yru"=yrw,
                       "v"=V[,2],"yrv"=V[,1],"yrsC"=yrsCalWindow,"nNeg"=nNeg,"nPos"=nPos,
                       "NcMin"=N3,"NextFigNumber"=NextFigNumber,"outputDir"=outputDir)
  # Until RecSLR1() finished, will not code the call with the above inputs; will just have that
  # function read MyData from wd.
  #save(SpecListLR,file="MyData")
  #---- stop
  
  
  # Call to function for method RecSLR1
  save(SpecListSLR1,file="MyDataSLR1")
  Z <- RecSLR1(SpecListSLR1)
  if (Z$flag>0){
    emssg<-Z$Msg
    ResTemp<-emssgUNH(emssg,outputDir)
    stop(emssg)
  }
} else if (methMSR==2 | methMSR==3){
  # Recon by regression on sreened SSRs or their PCs, with call to RecMLR1 
  #save(SpecListMLR1,file="MyData"). Note that lags have been dealt with at the SSR step. No
  # lags are included in the MSR model. But, nPos and nNeg are used in the MSR model to set
  # m in leave-m-out cross-validation, on grounds that the SSRs did use lagging.
  SpecListMLR1 <- list("Text"=RecListx,"U"=Y3,"yrU"=yrY3,"nmsU"=nms3,"jScreened"=jScreened,"v"=V[,2],"yrv"=V[,1],
                     "yrsC"=yrsCalWindow,"nNeg"=nNeg,"nPos"=nPos,"incR2a"=incR2a,"kstop"=kstop,
                     "NcMin"=N3, "PCoption"=PCoption,"f"=f,"PCApredictors"=PCApredictors,
                     "methMSR"=methMSR,"PdfDescribe"=PdfDescribe, "nPCsKeep"=nPCsKeep,"alphaR"=alphaR,
                     "ScreenAnalogPCs"=ScreenAnalogPCs, "kHowPCA"=kHowPCA,"NextFigNumber"=NextFigNumber,
                     "outputDir"=outputDir)
  save(SpecListMLR1,file="MyDataMLR1")
  Z <- RecMLR1(SpecListMLR1)
  if (Z$flag>0){
    emssg<-Z$Msg
    ResTemp<-emssgUNH(emssg,outputDir)
    stop(emssg)
  }
}

