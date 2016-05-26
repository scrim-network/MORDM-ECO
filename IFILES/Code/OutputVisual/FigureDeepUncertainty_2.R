# FigureDeepUncertainty_2.R
# Riddhi Singh
# rus197@psu.edu, Penn State 
# March, 2014

# R code for estimating uniformly distributed samples for sampling the parameters (mean and variance) of the lognormal distributions

# source code for the lognormal pdf and visualization from : https://stat.ethz.ch/pipermail/r-help/2003-April/032058.html accessed on 14th Nov, 2013
# edited by R Singh 

# Mise en place
rm(list = ls())
graphics.off()

#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

#define the lognormal functions
logmean <- function(location,scale) {
  mean_ln = exp(location+scale^2/2)
  return(mean_ln)
}

logmode <- function(location,scale) {
  mode_ln = exp(location-scale^2)
  return(mode_ln)
}

logmed  <- function(location,scale) {
  med_ln = exp(location)
  return(med_ln)
}

logvar  <- function(location,scale) {
  var_ln = (exp(scale^2)-1)*exp(2*location+scale^2)
  return(var_ln)
}

logskew  <- function(location,scale) {
  skew_ln = (exp(scale^2)+2)*sqrt(exp(scale^2)-1)
  return(skew_ln)
}

loglocation  <- function(meanval,varval) {
  loc_ln = log(meanval^2/(sqrt(varval+meanval^2)))
  return(loc_ln)
}

logscale  <- function(meanval,varval) {
  sca_ln = sqrt(log(1+(varval/meanval^2)))
  return(sca_ln)
}

genlognormal <- function(location, scale)
{
  num_pts    = 10000
  xrand      = rlnorm(num_pts, location, scale)  # for example
  rangeval   = range(xrand)
  grid       = seq(0,1,length=num_pts)
  pdfvals    = dlnorm(xrand, meanlog = mean(log(xrand)), sdlog = sd(log(xrand)))
  
  xrand_sort = sort(xrand,index=TRUE)
  
  #intialize
  x          = rep(0,length=length(xrand_sort))
  y          = rep(0,length=length(xrand_sort))
  logplot    = list(x=x, y=y)
  
  #assign
  logplot$x  = xrand[xrand_sort$ix]
  logplot$y  = pdfvals[xrand_sort$ix]
  
  return(logplot)
  
  #cleanup
  rm(xrand, rangeval, pdfrange, pdfvals, xrand_sort)
}
#specify the range of lognormal parameters - mean and variance
min_x               = 0.01
max_x               = 0.03
min_y               = -5    
max_y               = -6

range_x             = max_x - min_x
range_y             = max_y - min_y

#specify the total number of required samples :
ntot                = 9
nrows               = sqrt(ntot)
ncols               = sqrt(ntot)

#initialize the vectors and answer matrix
meanvals            = rep(0, length=nrows)
variancevals        = rep(0, length=ncols)
lognormal           = list(mean = mean, var = var)
lognormal$mean      = matrix(data=0, nrow=nrows, ncol=ncols)
lognormal$var       = matrix(data=0, nrow=nrows, ncol=ncols)
#now sample each parameter uniformly
int_x               = range_x/(nrows-1)
int_y               = range_y/(ncols-1)

#find the evenly spaced samples
for (row in 1:nrows)
{ 
  meanvals[row]     = min_x + int_x*(row-1)
}

for (col in 1:ncols)
{
  variancevals[col] = 10^(min_y + int_y*(col-1))
}

#generate the matrix
for (row in 1:nrows)
{ for (col in 1:ncols)
{ 
  lognormal$mean[row,col] = meanvals[row]
  lognormal$var[row,col]  = variancevals[col]
}
}

#location of original
centre_loc   = floor(sqrt(ntot)/2)+1

#initialize plot
#start plot 
plot.new()
par(family="sans")
par(cex=1.5)
par(mar=c(0,0.5,1,0.5)+0.1)
par(bg=NA)

mat = matrix(seq(1,2,length=2),2,1,byrow=TRUE)
layout(mat,widths=c(1,1),heights=c(1,1))
cols       = c(rgb(1,0,0),rgb(1,0.5,0),rgb(1,1,0),rgb(0.5,1,0),rgb(0,1,0), rgb(0,1,0.5), rgb(0,1,1), rgb(0,0.5,1), rgb(0,0,1))
indchange  = c(7,6,1,8,5,2,9,4,3)
cols       = cols[indchange]  #rearrange so that cols increase with mean
textlabs   = c("i","ii","iii","iv","v","vi","vii","viii","ix")

#plot the lognormal sample space
shift_x = 0.0025
plot(lognormal$mean+shift_x,log10(lognormal$var), type = 'p',pch=16, col=cols,xlab="",ylab="",cex=5,
     cex.lab=1.5,ylim=c(-6.8,-4.5),xlim=c(0.00,0.05),bty='n',axes=FALSE)
points(lognormal$mean[centre_loc,centre_loc]+shift_x,log10(lognormal$var[centre_loc,centre_loc]), type='p',pch=16,cex=3)
text(lognormal$mean+shift_x+0.0015,log10(lognormal$var)-0.1, textlabs[rev(indchange)],cex=1.5)
#plot the axes
axis(1,at=meanvals+shift_x, labels = meanvals,cex.axis=1.5,pos=-6.1-0.15)
axis(2,at=log10(variancevals), labels = log10(variancevals),cex.axis=1.5,pos=0.008+shift_x)

#plot label
text(0.004,-4.6,'(a)',cex=1.5)
#plot x and y label
text(0.004+shift_x,-5.5,'log10(Variance) [-]',srt=90,cex=1.5)
text(0.02+shift_x,-6.5-0.15,'Mean [-]',cex=1.5)
#plot title
text(0.02+shift_x,-4.7,'Deep uncertainty scenarios',cex=1.5)

#now plot these lognormals
for (col in 1:ncols)
{ for (row in 1:nrows)
{
  ind        = row + nrows*(col-1)
  location   = loglocation(lognormal$mean[row,col],lognormal$var[row,col])
  scale      = logscale(lognormal$mean[row,col],lognormal$var[row,col])
  plotdata   = genlognormal(location, scale)
  
  if ((row==1) && (col==1)) { 
    plot(plotdata$x, plotdata$y,xlab="Random pollution inflow [-]", ylab="Probability density [-]",xlim = c(-0.005,0.05),ylim=c(-150,500),
         cex.lab=1.5,col=cols[ind],type="l", mgp=c(2,1,0),lwd=4,bty='n',axes=FALSE,lty=1)
  } else {
    lines(plotdata$x,plotdata$y,lwd=4,col=cols[ind],lty=1)
  } 
}
}

#replot the original 
location   = loglocation(lognormal$mean[centre_loc,centre_loc],lognormal$var[centre_loc,centre_loc])
scale      = logscale(lognormal$mean[centre_loc,centre_loc],lognormal$var[centre_loc,centre_loc])
plotdata   = genlognormal(location, scale)
lnorg      = length(plotdata$x)
sel        = seq(1,lnorg,length=1000)
lines(plotdata$x[sel],plotdata$y[sel],lwd=4,col="black",lty=2)

#plot the axes
axis(1,at=seq(0.00,0.05,length=6), labels = seq(0.00,0.05,length=6),cex.axis=1.5,pos=0)
# axis(2,at=seq(0,500,length=2), labels = c('',''), cex.axis=1.5,pos=0)

#plot label
text(-0.001,500,'(b)',cex=1.5)
#plot x and y label
text(0.02,-120,'Mean [-]',cex=1.5)
#plot title
text(0.02,500,'Associated lognormals',cex=1.5)

#add legend
points(0.032,420, type='p',pch=16,cex=2,col="red")
points(0.0335,420, type='p',pch=16,cex=2,col="green")
points(0.035,420, type='p',pch=16,cex=2,col="blue")
points(0.0335,370, type='p',pch=16,cex=2,col="black")
text(0.044,420, "Uncertainty scenarios",cex=1.2) 
text(0.0430,370, "Known uncertainty",cex=1.2) 

dev.copy2pdf(file="FigureLognormals.pdf",width=10,height=10)