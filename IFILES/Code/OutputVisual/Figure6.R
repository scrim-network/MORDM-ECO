# Figure6.R
# Riddhi Singh
# rus197@psu.edu, Penn State 
# Febrauary, 2014
# 
# R code for plotting the solutions that are robust across a range of probabilities
# 
# # Mise en place
rm(list = ls())
graphics.off()

#load the data file

files  =c("LakeProblem_5ObjStochMM11AprilSelectedRSeed_Type4_lakestate_Type_green",
          "LakeProblem_5ObjStochMM11AprilSelectedRSeed_Type4_lakestate_Type_red",
          "LakeProblem_5ObjStochMM11AprilSelectedRSeed_Type4_lakestate_Type_blue",
          "LakeProblem_1ObjStochMM11April_Type4_lakestate_Type")           
path   = "Data/"
# import all data, normalize etc.
datalist = list()
mins     = matrix()
maxs     = matrix()
lnrun    = length(files)

#number of objectives to plot
endcols  = 100
startcols= 1

#order
order    = c(3, 2, 1, 6, 5, 4, 9, 8, 7)
for (filenum in 1:lnrun)
{
  for (dist in 1:9)
  {
    filename  <- paste(path,files[filenum],order[dist]-1, ".txt",sep="")
    print(filename)
    policydata <- read.table(filename)
  
    ##write the maxs and mins
    tempmins =  apply(policydata[,startcols:endcols],2,function(x) min(x))
    tempmaxs =  apply(policydata[,startcols:endcols],2,function(x) max(x))
  
    if (filenum==1)
    { mins = matrix(tempmins,1)
      maxs = matrix(tempmaxs,1)
    }
    else
    {
      mins = rbind(mins, tempmins);
      maxs = rbind(maxs, tempmaxs);
    }
    datalist[length(datalist)+1] = list(policydata[,startcols:endcols]);
    rm(tempmins,tempmaxs)
  }
}

########################################################
## A) Installing and loading required packages
########################################################

if (!require("plotrix")) {
  install.packages("plotrix", dependencies = TRUE)
  library(plotrix)
}


##write the maxs and mins
minval =  apply(mins,2,function(x) min(x))
maxval =  apply(maxs,2,function(x) max(x))
minval =  0
maxval = max(maxval)

#lognormal stuff
meanx    <- c(-4.610145,-4.620737,-4.652825,-3.913271,-3.915960,-3.924369,-3.507113,-3.508312,-3.512083)
meansd   <- c(0.09975135,0.17644567,0.3087235,0.04996879,0.08873899,0.1571388,0.03332408,0.05922401,0.1051182)
#the standard distribution
stdrand  <- rlnorm(3000, meanx[5],meansd[5])
stdrange <- range(stdrand)
stdpdfr  <- dlnorm(stdrange[1]:stdrange[2], meanlog = mean(log(stdrand)), sdlog = sd(log(stdrand)))
stdpdfv  <- dlnorm(stdrand, meanlog = mean(log(stdrand)), sdlog = sd(log(stdrand)))

pdf("Figure6.pdf", height=12, width=22)

#define border of the figure, giving some room for the legend on the right
par(family="sans")
par(cex=1.5)
par(mar=c(0,2.5,0,2.5))
par(mgp=c(1,2,0))
#define layout
widm    = 0.5
widl    = 0.1
mat     = matrix(seq(1,72,length=72),8,9,byrow=TRUE)
mat[8,] = 64
colr    = c(0,1,0.85)
layout(mat, widths=c(0.2,1,0.4,0.15,1,0.4,0.15,1,0.4,0.4),heights=c(0.1,widm,widl,widm,widl,widm,0.20,0.26))

#define plot labels
labs = c('i','ii','iii','iv','v','vi','vii','viii','ix','x','xi','xii')

terminal_year = 90
samples_plot  = 10000
xvals         = seq(1,terminal_year)

#plot the top row of blank figures
plot(1,type="n",axes=F,xlab="",ylab="")
plot(1,type="n",axes=F,xlab="",ylab="")
plot(1,type="n",axes=F,xlab="",ylab="")
plot(1,type="n",axes=F,xlab="",ylab="")
plot(1,type="n",axes=F,xlab="",ylab="")
plot(1,type="n",axes=F,xlab="",ylab="")
plot(1,type="n",axes=F,xlab="",ylab="")
plot(1,type="n",axes=F,xlab="",ylab="")
plot(1,type="n",axes=F,xlab="",ylab="")
#now plot individual figures
for (i in 1:9)
{
  
  if ((i==1) || (i==4) || (i==7)) {
    #plot y labels
    basey = 4
    basex = 8
    plot(seq(0,basex),rep(basey,each=basex+1),type="l",ylim=c(-5,13),xlim=c(0,15),lwd=6,col = "white",axes=FALSE, xlab="",ylab="",lty=2)
    text(basex,basey,"Phosphorus [-]",cex=3,srt=90)
  } 

  num           = i
  
  #postprocess the input lake states to include only the first 90 years out of 100 years from each solution type
  yvalsg        = datalist[[num]]
  yvalsg        = yvalsg[,1:terminal_year]
  
  yvalsr        = datalist[[9+num]]
  yvalsr        = yvalsr[,1:terminal_year]
  
  yvalsb        = datalist[[18+num]]
  yvalsb        = yvalsb[,1:terminal_year]
  
  yvals4        = datalist[[27+num]]
  yvals4        = yvals4[,1:terminal_year]
  
  useyval1 = yvalsg[1:samples_plot,1:terminal_year]
  useyval2 = yvalsr[1:samples_plot,1:terminal_year]
  useyval3 = yvalsb[1:samples_plot,1:terminal_year]
  useyval4 = yvals4[1:samples_plot,1:terminal_year]
  
  alphaval = 0.3
  if (i<7)
  {  plot(xvals,useyval1[1,],type="l",xlab="",ylab="",lwd= 4,ylim=c(minval,0.62),col=rgb(colr[1],colr[2],colr[3],0),
         cex.axis=4,cex.lab=4,xaxt="n",yaxt="n")
  }  else { 
    plot(xvals,useyval1[1,],type="l",lwd= 4,ylab = "",ylim=c(minval,2.1),col=rgb(colr[1],colr[2],colr[3],0),
         xlab="",cex.axis=4,cex.lab=4,xaxt="n",yaxt="n")
  }
  
  if (i<7)
    {
    apply(useyval1, 1, lines, col=rgb(0,0,1,alphaval))
    apply(useyval2, 1, lines, col=rgb(colr[1],colr[2],colr[3],alphaval))
    apply(useyval3, 1, lines, col=rgb(1,0,0,alphaval))
    apply(useyval4, 1, lines, col=rgb(1,0,0,alphaval), lty=2)
    lines(seq(0,terminal_year),rep(0.5,terminal_year+1), col="black",lty=2,lwd=2)
    axis.break(2,0.55,style="slash",brw=0.05)
  } else {
    apply(useyval1, 1, lines, col=rgb(0,0,1,alphaval))
    apply(useyval2, 1, lines, col=rgb(colr[1],colr[2],colr[3],alphaval))
    apply(useyval3, 1, lines, col=rgb(1,0,0,alphaval))
    apply(useyval4, 1, lines, col=rgb(1,0,0,alphaval), lty=2)
    lines(seq(0,terminal_year),rep(0.5,terminal_year+1), col="black",lty=2,lwd=2)
  }
  
  if ((i==7) || (i==8) || (i==9))
  axis(1,at=c(1,55,terminal_year), labels = c(1,50,terminal_year),cex.axis=3)
  
  if (i<7) {
    axis(2,at=c(0,0.4,0.60), labels = c(0,0.4,2),cex.axis=3)
  } else {
    axis(2,at=c(0,1,2), labels = c(0,1,2),cex.axis=3)
  }

  if (i==1)
    text(45,0.55,"Irreversible threshold",cex=3)
  
  #lognormal stuff
  xrand    <- rlnorm(3000, meanx[i], meansd[i])  # for example
  rangeval <- range(xrand)
  pdfrange <- dlnorm(rangeval[1]:rangeval[2], meanlog = mean(log(xrand)), sdlog = sd(log(xrand)))
  pdfvals  <- dlnorm(xrand, meanlog = mean(log(xrand)), sdlog = sd(log(xrand)))
  
  index1 = (sort(xrand, index.return=TRUE))
  index2 = (sort(stdrand, index.return=TRUE))
  index1 = unlist(index1[2])
  index2 = unlist(index2[2])
  
  plot(pdfvals[index1], xrand[index1],type="l",yaxt="n",xaxt="n",xlab='',ylab='', ylim=c(0,0.04), xlim=c(0,400),lwd=4)
  lines(stdpdfv[index2], stdrand[index2],col=rgb(0.7,0.7,0.7),lwd=4)
  axis(2,at=c(0,0.04), labels = c(0,0.04),cex.axis=3)
  text(325,0.037,labs[i],cex=4)

  if ((i==1) || (i==2) || (i==4) || (i==5) || (i==7) || (i==8)) {
    plot(1,type="n",axes=F,xlab="",ylab="")
  } else if ((i==3) || (i==6)) {
    plot(1,type="n",axes=F,xlab="",ylab="")
    plot(1,type="n",axes=F,xlab="",ylab="")
    plot(1,type="n",axes=F,xlab="",ylab="")
    plot(1,type="n",axes=F,xlab="",ylab="")
    plot(1,type="n",axes=F,xlab="",ylab="")
    plot(1,type="n",axes=F,xlab="",ylab="")
    plot(1,type="n",axes=F,xlab="",ylab="")
    plot(1,type="n",axes=F,xlab="",ylab="")
    plot(1,type="n",axes=F,xlab="",ylab="")
  } else {}
  
  #plot x labels
  if (i==9)
  { 
    basey = 4.5
    basex = 8
    plot(1,type="n",axes=F,xlab="",ylab="")
    plot(seq(basex-1,basex),rep(basey,each=2),type="l",ylim=c(4,6),xlim=c(5,10),lwd=6,col = "white",axes=FALSE, xlab="",ylab="",lty=2)
    text(basex,basey,"Time [years]",cex=3)
    plot(1,type="n",axes=F,xlab="",ylab="")
    plot(1,type="n",axes=F,xlab="",ylab="")
    plot(seq(basex-1,basex),rep(basey,each=2),type="l",ylim=c(4,6),xlim=c(5,10),lwd=6,col = "white",axes=FALSE, xlab="",ylab="",lty=2)
    text(basex,basey,"Time [years]",cex=3)
    plot(1,type="n",axes=F,xlab="",ylab="")
    plot(1,type="n",axes=F,xlab="",ylab="")
    plot(seq(basex-1,basex),rep(basey,each=2),type="l",ylim=c(4,6),xlim=c(5,10),lwd=6,col = "white",axes=FALSE, xlab="",ylab="",lty=2)
    text(basex,basey,"Time [years]",cex=3)
  }
}

#apply legends
plot(1,type="n",axes=F,xlab="",ylab="")
plot(0.5,0.5,col=rgb(1,1,1,1),xlim=c(0,1),
     mgp=c(0,0,0),ylim=c(0,0.9),axes=F,xlab="",ylab="")
ylow2     = 0.05
ylow      = 0.2
yhigh     = 0.6
lwidth    = 0.04
legwidth  = 0.05
xlocs     = seq(0.25,0.7,length=3)    #specify the location of the legend
xlocs     = c(0.18, xlocs)

lines(c(xlocs[2]-lwidth/2, xlocs[2]+lwidth/2), c(yhigh,yhigh), col=rgb(1,0,0), lwd= 4)
text(xlocs[2]+lwidth+0.05, yhigh, "Low phosphorus", cex=3)
  
lines(c(xlocs[2]-lwidth/2, xlocs[2]+lwidth/2), c(ylow, ylow), col=rgb(0,0,1), lwd= 4)
text(xlocs[2]+lwidth+0.035, ylow, "Compromise", cex=3)
  
lines(c(xlocs[3]-lwidth/2, xlocs[3]+lwidth/2), c(yhigh,yhigh), col=rgb(colr[1],colr[2],colr[3],1), lwd= 4)
text(xlocs[3]+lwidth+0.01, yhigh, "Utility", cex=3)
  
lines(c(xlocs[3]-lwidth/2, xlocs[3]+lwidth/2), c(ylow,ylow), col=rgb(1,0,0), lwd= 4, lty=2)
text(xlocs[3]+lwidth+0.05, ylow, "Single objective", cex=3)
  
lines(c(xlocs[4]-lwidth/2,xlocs[4]+lwidth/2), c(yhigh,yhigh), col=rgb(0.7,0.7,0.7), lwd=4)
text(xlocs[4]+lwidth+0.07,yhigh, "Baseline distribution",cex=3)
  
lines(c(xlocs[4]-lwidth/2,xlocs[4]+lwidth/2), c(ylow,ylow),  col="black", lwd=4)
text(xlocs[4]+lwidth+0.065,ylow, "Testing distribution",cex=3)

#colorbar
text(xlocs[1]-0.015, yhigh, "Legend", cex=3)
  
polygon(c(xlocs[1]-0.05,xlocs[3]+0.42,xlocs[3]+0.42,xlocs[1]-0.05),c(ylow2-0.01, ylow2-0.01, yhigh+0.2, yhigh+0.2), 
          border='gray', lwd =1,lty=1 ) 

dev.off()