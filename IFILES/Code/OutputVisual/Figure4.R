# Figure4.R
# Riddhi Singh
# rus197@psu.edu, Penn State 
# March, 2014

# R code for parallel coordinate plots

# Mise en place
rm(list = ls())
graphics.off()

#load the data file
files       = c("LakeProblem_5ObjStochMM11AprilSelectedRSeed_Type4_Reev_Type4","LakeProblem_5ObjStochMM11AprilSelectedRSeed_Type4_Reev_TypeAll")
path        = "Data/"

#import all data, normalize etc.
datalist    = list();
mins        = matrix();
maxs        = matrix();
lnrun       = length(files);

#number of objectives to plot
n.objs      = 5;
nvars       = 0;
endcols     = nvars+n.objs;
startcols   = nvars+1;
trans_power = 1
for (filenum in 1:lnrun)
{
  # read in the data
  filename      = paste(path,files[filenum],".txt",sep="")
  print(filename)
  my.rawdata     = read.table(filename)
  #apply -1 to correct for the sign on the objectives
  my.rawdata[,1] = (-1)*my.rawdata[,1]
  my.rawdata[,3] = my.rawdata[,3]^(trans_power)
  ##write the maxs and mins
  tempmins       =  apply(my.rawdata[,startcols:endcols],2,function(x) min(x))
  tempmaxs       =  apply(my.rawdata[,startcols:endcols],2,function(x) max(x))
  
  if (filenum==1)
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
minval     =  apply(mins,2,function(x) min(x))
maxval     =  apply(maxs,2,function(x) max(x))
print(minval)
print(maxval)

#now start plotting
#assign labels
lbs         = c("Phosphorus","Expected utility","Current DM","Future DM","Reliability");
lbs         = lbs[1:n.objs]

#define thresholds for robustness 
rangeval    = maxval - minval 
#threshold percentages for different objectives
norm_thres  = c(0.90,0.95,0.70,0.90,0.98)
#start plot 
plot.new()
par(family="sans")
par(cex=1.5)
par(mar=c(1,2,1,2)+0.1)
par(bg=NA)

layout(matrix(seq(1,3,length=3),3,1,byrow=TRUE))

minvalr = round(minval, digits=2)
maxvalr = round(maxval, digits=2)

labelx= c(seq(1,n.objs,length=n.objs))


#Plot 1 and 2 - ParallelCoordinate Plots
for (filenum in 1:lnrun)
{
  trial = datalist[[filenum]]
  for (i in 1:n.objs)
  {  
    colval = (trial[,i]-minval[i])/(maxval[i]-minval[i])
    
    if (i==1)
      newmat = matrix(colval,length(colval))
    else
      newmat = cbind(newmat, colval)
    
    rm(colval)
  }
  
  colnames(newmat) = lbs
  sz               = dim(datalist[[filenum]]);
  colnum           = sz[1];
  #load robustness file
  if (filenum==1) {
    robustness = unlist(read.table("Type4RobustnessSelectedSeed.txt"))
  } else {
    robustness = unlist(read.table("TypeAllRobustnessSelectedSeed.txt"))
  }
#   print(range(robustness))
  alphaval  = array(0, dim = c(length(robustness),1))
  lwid      = array(0, dim = c(length(robustness),1))
  colv      = array(0, dim = c(length(robustness),3))
  count     = 0
  for (ind in 1:length(robustness)) {
      if (robustness[ind]>0.70 && newmat[ind,5]>0.90) {
          alphaval[ind]  = 1
          lwid[ind]      = 1
          colv[ind, ]    = rgb(0.5,0.5,0.5)
          count          = count+1
      } else {
          alphaval[ind]  = 0.15
          lwid[ind]      = 1
          colv[ind, ]    = rgb(0.8,0.8,0.8)
      }
  }
  

  for (linenum in 1:colnum)
  {
    colv[linenum] = adjustcolor(colv[linenum], alphaval[linenum]) 
  }
  newmatnew  = newmat;
  
  if (filenum==1)
  { plot(newmatnew[1,], type="l",lwd = 0, col = colv[filenum,],cex.lab=1,ylim=c(-0.3,1.5),xlim=c(0.7,n.objs+1.5),axes=FALSE,xlab="",ylab="")
  } else if (filenum==2) {
    plot(newmatnew[1,], type="l",lwd = 0, col = colv[filenum,],cex.lab=1,ylim=c(-0.3,1.5),xlim=c(0.7,n.objs+1.5),axes=FALSE,xlab="",ylab="")
  }
  for (ind in 1:length(robustness)) {
    if (ind!=30 && ind!=92 && ind!=342) {
  lines(newmatnew[ind, ], col=colv[ind], lwd= lwid[ind])
    }
 }
print("robustness is ")
 print(robustness[30])

  #plotting the solutions that worked for the deep uncertainty
  lines(newmatnew[30,],   col=rgb(0,1,0, alphaval[30]),       lwd=4,  lty=1)
  if (filenum==1) {
    lines(newmatnew[92,], col=rgb(1,0,0, alphaval[92]),       lwd=4,  lty=1)   #red = 92, green= 30, blue =342
  } else {
    lines(newmatnew[92,], col=rgb(1,0,0, alphaval[92]+0.05),  lwd=4,  lty=1)   #red = 92, green= 30, blue =342  
  }
  lines(newmatnew[342,],  col=rgb(0,0,1, alphaval[342]), lwd=4, lty=1)
  
  for (i in 1:n.objs)
  {
    segments(i,0,i,1)
  }
  rm(newmatnew, alphaval, lwid, colv)
  #plot single objective Bentham
  if (filenum==1)
  {   files = c("LakeProblem_1ObjStochMM11April_Type4_Reev_Type4.txt")  
  } else if (filenum==2) {
    files = c("LakeProblem_1ObjStochMM11April_Type4_Reev_TypeAll.txt") }
  
  filename     = paste(path,files[1],sep="")
  stochBen     = read.table(filename)
  stochBen     = stochBen[, startcols:endcols]
  stochBen[,1] = -1*stochBen[,1]
  stochBen[,3] = stochBen[,3]^(trans_power)
  trial = stochBen
  
  for (i in 1:n.objs)
  {
    colval = (trial[,i]-minval[i])/(maxval[i]-minval[i])
    if (i==1)
      newmat = matrix(colval,length(colval))
    else
      newmat = cbind(newmat, colval)
    rm(colval)
  }
  if (filenum==1) {
    robustness  = 0.5960676
  } else {
    robustness  = 0.4126665
  }
  if (robustness>0.70 && newmat[1,5]>0.90) {
    alphaval  = 1
    lwid      = 1
    colv    = rgb(0.5,0.5,0.5)
    count          = count+1
  } else {
    alphaval  = 0.05
    lwid      = 1
    colv    = rgb(0.8,0.8,0.8)
  }

  newmatstoch = newmat
  lines(newmatstoch[1,], col=colv, lwd=4,lty=2)
  rm(alphaval, colv, lwid, count)
  #   par(new=T)
  #add the satisficing box for the standard risk scenario (single pdf)
  if (filenum==1)
  { wid = 0.19
    polygon(c(0.80,0.80,1.1,1.9,2.1,2.9,3.1,3.9,4.1, 4.9, 5+wid,5+wid), c(1+wid,norm_thres[1],norm_thres[1],norm_thres[2]-0.025,norm_thres[2]-0.025,norm_thres[3],norm_thres[3],
                                                                          norm_thres[4],norm_thres[4],norm_thres[5],norm_thres[5], norm_thres[5]+wid), border='black', lwd =2,lty=2 ) }
  #   finally plot the labels
  for (i in 1:n.objs)
  { if (i==1) { text(i,-0.1, 'Max',cex=1.2) 
  } else {
    text(i,-0.1, 'Min',cex=1.2) }
  }
  for (i in 1:n.objs)
  { if (i==1) { text(i,1.1, 'Min',cex=1.2) 
  } else {
    text(i,1.1, 'Max',cex=1.2) }
  }
  
  #apply titles
  if (filenum==1) {
    text(3,1.3, 'Known PDF',cex=1.2)
  } else if (filenum==2) {
    text(3,1.4, 'Multiple PDFs (average)',cex=1.2)
  }
  
  #apply plot labels
  if (filenum==1) {
    text(0.7,1.3, '(a)',cex=1.2)
  } else if (filenum==2) {
    text(0.7,1.4, '(b)',cex=1.2)
  }
  
  #write the labels
  for (i in 1:n.objs)
  {
    text(i,-0.27, paste("",lbs[i],sep=""),family="sans",cex=1.2)      
  }
  #apply preference arrow and associated text
  arrows(0.7,0,0.7,1)
  text(0.6,0.4,"Preference",srt=90,cex=1.2)
  text(6.1,0.4,"%Robust > 70%",cex=1.2)
  text(6.1,0.6,"Reliability > 90%",cex=1.2)
  lines(c(5.2,5.6),c(0.5,0.5),col=rgb(0.4,0.4,0.4),lwd=1.2)
  rm(trial,newmat,colnum,sz)
}

#clean everything again
rm(list=setdiff(ls(),"trans_power"))     #delete everything except the power for transformation 

#plot 3 - policies
#load the data file
files    =c("LakeProblem_5ObjStochMM11AprilSelectedRSeed_Type4_Par_Pareto","LakeProblem_1ObjStochMM11April_Type4_Par_Pareto")
path        = "Data/"

#import all data, normalize etc.
datalist = list();
mins     = matrix();
maxs     = matrix();
lnrun    = length(files);

#number of objectives to plot
nvars    = 20;
endcols  = nvars;
startcols= 1;

for (filenum in 1:lnrun)
{
  #   filename =files[filenum];
  filename  <- paste(path,files[filenum],".txt",sep="")
  print(filename)
  my.rawdata <- read.table(filename)
  
  ##write the maxs and mins
  tempmins =  apply(my.rawdata[,startcols:endcols],2,function(x) min(x))
  tempmaxs =  apply(my.rawdata[,startcols:endcols],2,function(x) max(x))
  
  if (filenum==1)
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
# 
# plot.new()
par(xpd=T)
par(family = "sans")
par(font=1)
par(cex=1.0)
par(bg=NA)

#define border of the figure, giving some room for the legend on the right
par(mar=c(2,3,1,3)+0.1)
#define distance between x and y axis titles and plot
par(mgp=c(0.0,0.4,0))

#indices of chosen solutions to be highlighted
indexnum1      = 30
indexnum2      = 92
indexnum3      = 342
# indexrdm       = 392
sortdata       = datalist[[1]]       
newmat1        = sortdata[indexnum1,]
newmat2        = sortdata[indexnum2,]
newmat3        = sortdata[indexnum3,]
# newmatnew      = sortdata[indexrdm,]

newmat1        = rep(newmat1[1,1:20], each=5)
newmat2        = rep(newmat2[1,1:20], each=5)
newmat3        = rep(newmat3[1,1:20], each=5)

# newmatrdm      = matrix(0,100,length(indexrdm))
# for (i in 1:length(indexrdm))
#   newmatrdm[,i]        = as.numeric(rep(newmatnew[i,], each=5))

newmat1        = as.numeric(newmat1)
newmat2        = as.numeric(newmat2)
newmat3        = as.numeric(newmat3)
# newmatrdm      = as.numeric(newmatrdm)

terminal_year  = 90
xvals          = seq(1,terminal_year)

newsortdata    = apply(sortdata, 1 , rep, each=5)
newsortdata    = newsortdata[1:terminal_year,]
#plotting the policies
plot(xvals,newsortdata[,1],type="l",lwd=2,xlab="Time [years]",ylab = "",ylim=c(-0.020,0.1),xlim=c(0,95),
     col="green", cex.axis=0.8,cex.lab=0.8,axes=FALSE,)

#overlay all other solutions
apply(newsortdata, 2 , lines, lwd=2, col=rgb(0.8,0.8,0.8))
#overlay chosen solutions
lines(seq(1,terminal_year),newmat1[1:terminal_year],lwd=4,col=rgb(0,1,0))
lines(seq(1,terminal_year),newmat2[1:terminal_year],lwd=4,col=rgb(1,0,0))
lines(seq(1,terminal_year),newmat3[1:terminal_year],lwd=4,col=rgb(0,0,1))

# for (i in 1:length(indexrdm))
#   lines(seq(1,terminal_year),newmatrdm[1:terminal_year,i],lwd=4,col=rgb(0,1,0))


#overlay the Bentham single objective
newmat4       = datalist[[2]]
newmat4       = rep(newmat4[1,1:20],each=5)
lines(seq(1,terminal_year),newmat4[1:terminal_year],lwd=4,col=rgb(0.4,0.4,0.4),lty=2)

#plot title
text(45,0.1, 'Strategies',cex=0.8)

#apply plot labels
text(-8,0.11, '(c)',cex=0.8)

#put y label
text(-8, 0.045,"Pollution [-]", srt=90,cex=0.8)

axis(1, at=c(0,50,90),lab=c('1','50','90'),pos=0,cex.axis=0.8)
axis(2, at=c(0,0.05,0.1),pos=0,cex.axis=0.8)
# dev.off()
filename      = paste("Figure7Revisions2015Proof",trans_power,".pdf",sep="")

dev.copy2pdf(file=filename,height=10,width=8)
