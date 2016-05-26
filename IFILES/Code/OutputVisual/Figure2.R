# Figure2.R
# Riddhi Singh
# rus197@psu.edu, Penn State 
# June, 2014

# R code to plot solutions shaded based on their robustness for the well known
# uncertainty case


#########################################################
### A) Installing and loading required packages
#########################################################
if (!require("plotrix")) {
  install.packages("plotrix", dependencies = TRUE)
  library(plotrix)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

#########################################################
### B) Main body
#########################################################

# Mise en place
rm(list = ls())
graphics.off()

#load the data file
files       = c("LakeProblem_5ObjStochMM11AprilSelectedRSeed_Type4_Reev_Type4")
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
  #   filename =files[filenum];
  filename      = paste(path,files[filenum],".txt",sep="")
  my.rawdata     = read.table(filename, colClasses="numeric")
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

#assign the minimum reliability to that of the single objective
minval[5]  = 0.98033800000000004
minmax     = rbind(minval,maxval)

print(round(minval,2))
print(round(maxval,2))
#write minvals and maxvals for baseline to output file
write.table(minmax,sep = " ",eol="\n ",file="MinMaxValsBaseline.txt",row.names=FALSE, col.names=FALSE)

#assign labels
lbs         = c("Phosphorus","Expected utility","Current SH","Future SH","Reliability"
                ,"Economic activity","Inertia");
lbs         = lbs[1:n.objs]

#define thresholds for robustness 
rangeval      = maxval - minval 
#threshold percentages for different objectives
norm_thres     = read.table("Exploratory_RobustnessThresholds.txt",colClasses="numeric")

norm_thres[1]  = -norm_thres[1]
norm_thres[1]  = (norm_thres[1]-minval[1])/(maxval[1]-minval[1])
norm_thres[1]  = max(unlist(c(0, norm_thres[1])))
#setting utility thresholds to 0%
norm_thres[2]  = 0
norm_thres[3]  = 0
norm_thres[4]  = 0
norm_thres[5]  = (norm_thres[5]-minval[5])/(maxval[5]-minval[5])
norm_thres[6]  = (norm_thres[6]-minval[6])/(maxval[6]-minval[6])
norm_thres[7]  = (norm_thres[7]-minval[7])/(maxval[7]-minval[7])

#now start plotting
#start plot 
plot.new()
par(family="sans")
par(cex=1.5)
par(mar=c(0,0,0,0)+0.1)
par(bg=NA)

layout(matrix(seq(1,1,length=1),1,1,byrow=TRUE))

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
  sz               = dim(datalist[[filenum]])
  colnum           = sz[1]
  
  #method to find the compromise: maximize the minimum, consider only objectives over 
  #which the algorithm was optimized
  #find the minimum objective of each set
  findcompromise   = array(0, dim=c(colnum,5))
  findcompromise   = newmat[,1:5]
  min_obj          = apply(findcompromise, 1, min)
  #find the maximum of the minimums
  ind              = sort(min_obj,decreasing=TRUE, index.return=TRUE)
  indcompromise1   = ind$ix[1]     #assign the maximum of the minimums -  this is the compromise solution
  
  #load robustness file
  if (filenum==1) {
    robustness = unlist(read.table("Type4RobustnessSelectedSeedAll.txt",colClasses="numeric"))
  } else {
    robustness = unlist(read.table("TypeAllRobustnessSelectedSeedAll.txt",colClasses="numeric"))
  }
  
  alphaval    = array(0.25, dim = c(length(robustness),1))
  lwid        = array(1, dim = c(length(robustness),1))
  count       = 0
  num         = 100
  colall      = array(0, dim = c(num,3))
  colbase     = array(0, dim = c(num,3))
  
  colall      = (rainbow(num, s=1, v=1, start=0, end = 1/1.4,alpha=1)) # rainbow colors
  sortvals    = sort(robustness, index.return=TRUE)
  
  indtoassign         = array(0, dim=c(length(robustness),1))
  indtoassign         = round(robustness*100)
  zeros               = which(indtoassign==0)
  ones                = which(robustness==1)
  indtoassign[zeros]  = 1 
  
  for (linenum in 1:colnum)
  {
    colbase[linenum] = adjustcolor(colall[indtoassign[linenum]], alphaval[linenum]) 
    if (linenum==indcompromise1 || linenum==342 || linenum==92) 
    {colbase[linenum] = adjustcolor(colall[indtoassign[linenum]], 1)}
  }
  newmatnew  = newmat
  
#########################################################################################################
##this part finds the compromise solution after removing the solutions that are not robust#####
#   #now estimate the compromise solution after removing those which are not robust
#   surviving  = array(0, dim=c(length(ones),n.objs))
#   survIndex  = array(0, dim=c(length(ones),1))
#   count      = 0
#   for (linenum in 1:colnum)
#   {
#     if (robustness[linenum]==1) {
#       count             = count+1
#       survIndex[count]  = linenum
#       surviving[count,] = newmat[linenum,1:n.objs]
#     }     
#   }
#   
#   write.table(survIndex,sep = " ",eol=" ",file="WellKnownPDFRobust.txt",row.names=FALSE, col.names=FALSE)
#   rm(min_obj, ind)
#   #find the minimum objective of each set
#   min_obj          = apply(surviving[,1:5], 1, min)
#   #find the maximum of the minimums
#   ind              = sort(min_obj,decreasing=TRUE, index.return=TRUE)
#   indcompromise2   = ind$ix[1]     #assign the maximum of the minimums -  this is the compromise solution
#########################################################################################################  
  
  plot(newmatnew[1,], type="l",lwd = 0, col = colbase[1],cex.lab=1,ylim=c(-0.72,1.3),xlim=c(0.7,n.objs+0.2),axes=FALSE,xlab="",ylab="")
  
  for (ind in 1:length(robustness)) { lines(newmatnew[ind, ], col=colbase[ind], lwd= lwid[ind]) }
  
  print(round(newmatnew[indcompromise1,],3))
  print(round(newmatnew[342,],3))
  print(round(newmatnew[92,],3))
  
  #highlight the compromise solution before robustness assessment
  lines(newmatnew[indcompromise1, ], col=colbase[indcompromise1], lwd= 4)
  points(newmatnew[indcompromise1, ], col=rgb(0,0,0), pch=21, bg=colbase[indcompromise1], cex=3, lwd=4)
  
  #highlight the compromise solution after robustness assessment
  #   lines(surviving[indcompromise2, ], col="yellow", lwd= 1.5)
  #   points(newmatnew[indcompromise1, ], col=colbase[ind], pch=1, cex=1, lwd=1 )
  
  #highlight the min phosphorus solution
  lines(newmatnew[342,],  col=colbase[342], lwd=4)
  points(newmatnew[342,], col=rgb(0,0,0), pch=23, bg = colbase[342], cex=3, lwd=4)
  
  #highlight the max utility solution
  lines(newmatnew[92,],   col=colbase[92],    lwd=4)
  points(newmatnew[92,],  col=rgb(0,0,0), pch=22, bg= colbase[92], cex=3, lwd=4)
  
  for (i in 1:n.objs) { segments(i,0,i,1) }
  
  rm(newmatnew, alphaval, lwid)

  #plot single objective Bentham
  if (filenum==1)
  {   files = c("LakeProblem_1ObjStochMM11April_Type4_Reev_Type4.txt")  
  } else if (filenum==2) {
    files = c("LakeProblem_1ObjStochMM11April_Type4_Reev_TypeAll.txt") }
  
  filename     = paste(path,files[1],sep="")
  stochBen     = read.table(filename,colClasses="numeric")
  stochBen     = stochBen[, startcols:endcols]
  stochBen[,1] = -1*stochBen[,1]
  
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
  print(round(newmat[1,],3))
  if (filenum==1) {
    robustness  = 0
  } else {
    robustness  = 0
  }
  #   robustness  = (robustness-0.66)/(1-0.66)
  indtocol = round(robustness*100)
  if (indtocol == 0) {indtocol = 1}
  newmatstoch = newmat
  
  lines(newmatstoch[1,], col=colall[indtocol], lwd=4)
  points(newmatstoch[1,], col=rgb(0,0,0), pch = 24, bg = colall[indtocol], lwd=4, cex=3)
  rm(count)
  
  #add the satisficing box for the standard risk scenario (single pdf)
#   if (filenum==1)
#   { wid = 0.1
#     polygon(c(0.80,0.80,1.1,1.9,2.1,2.9,3.1,3.9,4.1, 4.9,5.1,5.9, 6+wid,6+wid), 
#             c(0.95+wid,0,0,norm_thres[2],norm_thres[2],norm_thres[3],
#               norm_thres[3],norm_thres[4],norm_thres[4],norm_thres[5],norm_thres[5], 
#               norm_thres[6],norm_thres[6],0.95+wid), 
#               border='black', lwd =2,lty=2 ) }
  
  #instead of a window, add tolerable thresholds 
#   points(1, norm_thres[1],pch=4, cex=3, lwd=3)
#   points(5, norm_thres[5],pch=4, cex=3, lwd=3)
  
  #finally plot the labels
  for (i in 1:n.objs)
  { if (i==1) { text(i,-0.1, bquote("Max"[b]),cex=1.2) 
  } else {
    text(i,-0.1, bquote("Min"[b]),cex=1.2) }
  }
  for (i in 1:n.objs)
  { if (i==1) { text(i,1.105, bquote("Min"[b]),cex=1.2) 
  } else {
    text(i,1.105, bquote("Max"[b]),cex=1.2) }
  }
  
  #apply titles
  if (filenum==1) {
    text(n.objs/2+0.5,1.25, 'Well-characterized uncertainty',cex=1.2)
  } else if (filenum==2) {
    text(3,1.4, 'Multiple PDFs (average)',cex=1.2)
  }
  
  #write the labels
  for (i in 1:n.objs)
  {
    text(i,-0.18, paste("",lbs[i],sep=""),family="sans",cex=1.2)      
  }
  
  #apply preference arrow and associated text
  arrows(0.7,0,0.7,1)
  text(0.6,0.4,"Preference",srt=90,cex=1.2)

  #apply legends
  ylow      = -0.72
  yhigh     = -0.62
  lwidth    = 0.2
  xlocs     = c(2,3.3)     #specify the location of the legend
 
  lines(c(xlocs[1]-lwidth/2, xlocs[1]+lwidth/2), c(yhigh,yhigh), col=rgb(0,0,0), lwd= 2)
  points(c(xlocs[1]-lwidth/2,xlocs[1]+lwidth/2),c(yhigh,yhigh), col=rgb(0,0,0), pch=23, cex=2, lwd=2)
  text(xlocs[1]+lwidth, yhigh, "Low phosphorus", cex=1,pos=4)

  lines(c(xlocs[1]-lwidth/2, xlocs[1]+lwidth/2), c(ylow, ylow), col=rgb(0,0,0), lwd= 2)
  points(c(xlocs[1]-lwidth/2,xlocs[1]+lwidth/2),c(ylow, ylow), col=rgb(0,0,0), pch=21, cex=2, lwd=2)
  text(xlocs[1]+lwidth, ylow, "Compromise", cex=1, pos=4)

  lines(c(xlocs[2]-lwidth/2, xlocs[2]+lwidth/2), c(yhigh,yhigh), col=rgb(0,0,0), lwd= 2)
  points(c(xlocs[2]-lwidth/2,xlocs[2]+lwidth/2),c(yhigh,yhigh), col=rgb(0,0,0), pch=22, cex=2, lwd=2)
  text(xlocs[2]+lwidth, yhigh, "Utility", cex=1, pos=4)

  lines(c(xlocs[2]-lwidth/2, xlocs[2]+lwidth/2), c(ylow,ylow), col=rgb(0,0,0), lwd= 2)
  points(c(xlocs[2]-lwidth/2,xlocs[2]+lwidth/2), c(ylow,ylow), col=rgb(0,0,0), pch=24, cex=2, lwd=2)
  text(xlocs[2]+lwidth, ylow, "MEU (single objective)", cex=1, pos=4)
  
  text(xlocs[1], yhigh+0.08, "Legend", cex=1)

  polygon(c(xlocs[1]-0.2,xlocs[2]+1.05,xlocs[2]+1.05,xlocs[1]-0.2),c(ylow-0.06, ylow-0.06, yhigh+0.12, yhigh+0.12), 
          border='gray', lwd =1,lty=1 ) 
  
  #colorbar
  leghei      = 0.05
  legwidth    = 3.0
  legy        = -0.34
  legx        = 1.5
  par("fg"="black")
  color.legend(legx,legy,legx+legwidth,legy+leghei/1.2,legend=c("",""),
             rect.col=colall,gradient="x",cex=1)
  text(legx,legy-0.028,round(min(0.0*100)),srt=0,cex=1)
  text(legx+legwidth,legy-0.028,round(max(1*100)),srt=0,cex=1)
  text(legx+legwidth/2, legy+leghei/2,"%Robust",srt=0,cex=1,pos=3)
  text(legx+legwidth/2, legy-1.3*leghei,"%Robust = %SOWs where phosphorus < 0.21 & reliability > 0.99 & economic activity > 0.5*Optimal",
     srt=0,cex=1, pos=1)

  rm(trial,newmat,colnum,sz)
}

#clean everything again
rm(list=setdiff(ls(),"trans_power"))     #delete everything except the power for transformation 

filename      = paste("Figure2b",".pdf",sep="")
dev.copy2pdf(file=filename,height=9,width=14)