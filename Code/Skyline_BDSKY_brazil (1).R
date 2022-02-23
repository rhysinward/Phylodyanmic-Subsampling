# Install the package by uncommenting the next line and running it
#
#install.packages('devtools')
#devtools::install_github("laduplessis/bdskytools")
#
# If you cannot install the package, source the R-files in the folder and make sure the 
# packages boa and RColorBrewer are installed

# Load the package (if successfully installed)
library(bdskytools)


# Extract the data and HPDs from the logfile
############################################

# Set the working directory to the directory where your log files are stored 
# (on RStudio navigate to Session > Set Working Directory > Choose Directory)
# or change "hcv_bdsky.log" to the path of the BDSKY logfile on your computer
fname <- "inverse_Brazil.log"
lf    <- readLogfile(fname, burnin=0.1)

Re_sky    <- getSkylineSubset(lf, "reproductiveNumber")
Re_hpd    <- getMatrixHPD(Re_sky)
delta_hpd <- getHPD(lf$becomeUninfectiousRate)


# The non-gridded intervals
###########################
# Since the intervals in this analysis are equidistant between the origin and the present 
# and the origin was also estimated, they should probably not be plotted in this way.
#
# Use this only when the interval times in bdsky are fixed 
# (needs to be manually set in the xml at the moment)
# or when the origin is fixed (automatically fixes the shift-times)
#
# When doing this do NOT plot with type='smooth' as it gives a misleading result!!!
plotSkyline(1:10, Re_hpd, type='step', ylab="R")

# Why is the result below misleading?
plotSkyline(1:10, Re_hpd, type='smooth', ylab="R")


# Extract gridded HPDs  
######################
# The trimegrid specifies the years into the past to grid the skyline to.
# We grid the skyline to every 4th year for the past 400 years (100 grid cells)
# If the timegrid is too fine it will not be that smooth anymore, 
# but it's fast enough to quickly interpolate a grid of 400 (every year for the past 400 years)
timegrid       <- seq(0,0.167,length.out=69)
Re_gridded     <- gridSkyline(Re_sky,    lf$origin, timegrid)
Re_gridded_hpd <- getMatrixHPD(Re_gridded)

save(Re_gridded_hpd, file = "inverse.R")

# Plotting the skyline
#####################
# The calendar times of the grid (the most recent sample is from 1993)
times     <- 2021.099-timegrid
plotSkyline(times, Re_gridded_hpd, type='smooth', xlab="Time", ylab="R")
plotSkyline(times, Re_gridded_hpd, type='lines', xlab="Time", ylab="R")


# Pretty skyline plots
######################
# Use this function for a more polished almost publication-ready plot
plotSkylinePretty(times, Re_gridded_hpd, axispadding=0.0, col=pal.dark(corange), fill=pal.dark(corange, 0.5), col.axis=pal.dark(corange),
xlab="Time", ylab=expression("R"[e]), side=2, yline=4, xline=2, xgrid=TRUE, ygrid=TRUE, gridcol=pal.dark(cgray), ylims=c(0,4))

# Add a line at Re = 1 for comparison
abline(h=1, col=pal.dark(cred), lty=2)

# Plot both Re and delta skylines on one set of axes
# Delta only has a dimension of 1, so the skyline is not really insightful
# Can also use this to compare skylines of the same parameter between different models (eg. changing the priors or number of shifts)
par(mar=c(5,4,4,4)+0.1)
plotSkylinePretty(range(times), as.matrix(delta_hpd), type='step', axispadding=0.0, col=pal.dark(cblue), fill=pal.dark(cblue, 0.5), col.axis=pal.dark(cblue),
ylab=expression(delta), side=4, yline=2, ylims=c(0,1), xaxis=FALSE)
plotSkylinePretty(times, Re_gridded_hpd, type='smooth', axispadding=0.0, col=pal.dark(corange), fill=pal.dark(corange, 0.5), col.axis=pal.dark(corange),
xlab="Time", ylab=expression("R"[e]), side=2, yline=2.5, xline=2, xgrid=TRUE, ygrid=TRUE, gridcol=pal.dark(cgray), ylims=c(0,3), new=TRUE, add=TRUE)

