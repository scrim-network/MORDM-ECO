# FigureA1RandomSeedPlot.R
# Riddhi Singh
# rus197@psu.edu, Penn State 
# March, 2014

# R code to identify robust solutions w.r.t objective space for each random seed
# Then plot the robust solutions
# there are 2 possible approaches:
# There are two methods which can be chosen by the rdm_method flag
# method 1 - maximize the minimum performance across all objectives per pareto set
# method 2 - minimize the variation between different objectives (normalized) and maximize the mean across all
#            or minimize the coefficient of variation per pareto set

#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

if (!require("plotrix")) {
  install.packages("plotrix", dependencies = TRUE)
  library(plotrix)
}
#########################################################
### B) Main body
#########################################################

# Mise en place
rm(list = ls())
graphics.off()

#load the data file
files    = c("LakeProblem_5ObjStochMM11April_Type4_RSeed")
path     = "Data/"

#import all data, normalize etc.
datalist = list()
mins     = matrix()
maxs     = matrix()
nseeds   = 10

#number of objectives to plot
n.objs   = 20
nvars    = 0
endcols  = nvars+n.objs
startcols= nvars+1
indtouse = array(0, dim = nseeds)

for (rseed in 1:nseeds)
{
  filename  <- paste(path,files,rseed-1,"_Par_Pareto.txt",sep="")
  print(filename)
  my.rawdata <- read.table(filename)
  
  ##write the maxs and mins
  tempmins =  apply(my.rawdata[,startcols:endcols],2,function(x) min(x))
  tempmaxs =  apply(my.rawdata[,startcols:endcols],2,function(x) max(x))
  
  if (rseed==1)
  { mins = matrix(tempmins,1)
    maxs = matrix(tempmaxs,1)
  }
  else
  {
    mins = rbind(mins, tempmins);
    maxs = rbind(maxs, tempmaxs);
  }
  datalist[length(datalist)+1] = list(my.rawdata[,startcols:endcols]);
  rm(my.rawdata)
  rm(tempmins,tempmaxs)
}

##write the maxs and mins
minval =  apply(mins,2,function(x) min(x))
maxval =  apply(maxs,2,function(x) max(x))

indtouse  = c(377 ,351 ,115 ,246 ,386 ,352  ,58 ,259 ,191 ,218)

newmat    = array(0, dim=c(nseeds, n.objs))
for (rseed in 1:nseeds)
{
  #assign array to be plotted
  newmat[rseed,]         = unlist(datalist[[rseed]][indtouse[rseed],])
}

#########################################################
### C) Customizing and plotting the plot
#########################################################
#initialize plot
#start plot 
plot.new()
par(family="sans")
par(cex=1.5)
par(mar=c(5,5,2,2)+0.1)
par(bg=NA)

#define number of points
nYears         = 90
xvals          = seq(1,100)

plot.window(xlim=c(1,nYears), ylim=c(0,0.1))

#determining plotting stuff for the surface
maxall = apply(newmat[,1:20],2,function(x) max(x))
minall = apply(newmat[,1:20],2,function(x) min(x))
maxall = rep(maxall,each=5)
minall = rep(minall,each=5)
plotdat= apply(newmat, 1 , rep, each=5)
plotdat= plotdat[1:nYears,]
plotdat= t(plotdat)

plot(xvals[1:nYears],rep(0, length=nYears),type="l",xlab="Time [years]",ylab = "Pollution [-]",
     ylim=c(0,0.1),col="white", cex.axis=1,cex.lab=1, mgp=c(2.5,1,0))
apply(plotdat, 1, lines, col="black")
polygon(c(xvals[1:nYears],rev(xvals[1:nYears])),c(maxall[1:nYears],rev(minall[1:nYears])),col=rgb(0.5,0.5,0.5,0.5),border=NA)

lines(c(52,58),c(0.08,0.08),lty=1,lwd=3)
lines(c(52,58),c(0.07,0.07),lty=1,lwd=3)
text(75, 0.08, "Random Seeds", cex=1)
text(75, 0.07, "Upper/lower limits", cex=1)
dev.copy2pdf(file="AppendixA1.pdf",height=10, width=10)