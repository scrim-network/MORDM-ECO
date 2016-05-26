# FigureA2.R
# Riddhi Singh
# rus197@psu.edu, Penn State 
# March, 2014

# R code for plotting runtime dynamics

# Mise en place
rm(list = ls())
graphics.off()

#Installing and loading required packages
if (!require("animation")) {
  install.packages("animation", dependencies = TRUE)
  library(animation)
}

if (!require("scatterplot3d")) {
  install.packages("scatterplot3d", dependencies = TRUE)
  library(scatterplot3d)
}
if (!require("aqfig")) {
  install.packages("aqfig", dependencies = TRUE)
  library(aqfig)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("plotrix")) {
  install.packages("plotrix", dependencies = TRUE)
  library(plotrix)
}
#load the data file
files       = c("runtime5Obj1ConstStochMM_Type4RSeed_0_0")
path        = "Data"

#import all data, normalize etc.
datalist    = list();
mins        = matrix();
maxs        = matrix();
lnrun       = length(files);

#number of objectives to plot
n.objs      = 5;
nvars       = 20;
endcols     = nvars+n.objs;
startcols   = nvars+1;

for (filenum in 1:lnrun)
{
  #   filename =files[filenum];
  filename      = paste(path,files[filenum],".txt",sep="")
  print(filename)
  my.rawdata     = file(filename, open="r")
  textline       = readLines(my.rawdata)
  ln             = length(textline)
  
  #first 4 lines are headers
  counter        = 4;
  nruns          = 0;
  
  counter        = counter+1;
  tempvar        = strsplit(textline[counter]," ") 
  
  #declare vectors 
  nfe            = vector();
  elt            = vector();
  sbx            = vector();
  de             = vector();
  pcx            = vector();
  spx            = vector();
  undx           = vector();
  um             = vector();
  
  #declare objective list
  objlist        = list()
  while (counter<ln)
  {
    if (tempvar=="#")
    { 
      counter      = counter+1;
      tempvar      = strsplit(textline[counter]," ")
      tempvar      = as.numeric(unlist(tempvar))
      nruns        = nruns+1;
      nfe[nruns]   = tempvar[1];
      elt[nruns]   = tempvar[2];
      sbx[nruns]   = tempvar[3];
      de[nruns]    = tempvar[4];
      pcx[nruns]   = tempvar[5]
      spx[nruns]   = tempvar[6];
      undx[nruns]  = tempvar[7];
      um[nruns]    = tempvar[8];
      objlist[[nruns]] = list();
      nsets        = 0;
      counter      = counter+2;
      tempvar      = textline[counter]
      if (tempvar=="//Empty")      #to handle empty sets
      {
        counter    = counter+1;
        tempvar    = textline[counter]
      } else {
        tempvar      = strsplit(textline[counter]," ")
        tempvar      = as.numeric(unlist(tempvar))
        while ((length(tempvar)==25) && (counter<ln))
        {
          nsets       = nsets+1
          objlist[[nruns]][nsets] = list(tempvar[startcols:endcols]);
          counter     = counter+1;
          tempvar     = textline[counter]
          tempvar     = strsplit(tempvar," ")
          tempvar     = unlist(tempvar)
         if (length(tempvar)==25)
           {tempvar      = as.numeric(tempvar)}
        }
      }
    }
  }
  rm(my.rawdata)
}

# now plotting the runtime dynamics
# number of nfe's
lnx = length(objlist)
# 
png(file="PngFolder/Runtime%07d.png", width=400, height=400)
for (nrun in 1:lnx)
        { 
           lnt     = length(objlist[[nrun]])           
           tempar  = array(data=NA, dim=c(lnt,n.objs))
           sortedar= array(data=NA, dim=c(lnt,n.objs))
           
           #to handle empty sets, sets  with only one front, and the multiple fronts
           if (lnt==0) {
             temp         = 0;
             my_palette   =rgb(1,1,1)
             sd3d         = scatterplot3d(0, 0, 0,            
                                          xlab='Expected utiltiy [-]', ylab='', zlab='Phosphorus in the lake [-]',
                                          xlim = c(0,0.45),ylim=c(-0.01,0.01),zlim=rev(c(0,0.35)),color=my_palette, angle=60,mgp=c(0.5,0.0,0.0),tck=0.02,
                                          mar=c(4,4,1,1),cex.lab=1.0)             
             sd3d$points3d(0.50,0.01,0.00,pch=16, col= "red",cex=2)
             rm(sd3d)
           } else if (lnt==1) {
                         temp             = objlist[[nrun]][1]
                         temp             = unlist(temp)
                         temp[2]          = -temp[2]
                         temp[3]          = -temp[3]
                         temp[4]          = -temp[4]
                         temp[5]          = -temp[5]
                        
                          #define colors
                          if (temp[4]<0){
                          my_palette      = rgb(0,0,1)
                          }  else {
                          my_palette      = rgb(1,0,0)
                          }
                            
                         sd3d =   scatterplot3d(temp[2], temp[4], temp[1],            
                                     xlab='Expected utiltiy [-]', ylab='', zlab='Phosphorus in the lake [-]',
                                     xlim = c(0,0.45),ylim=c(-0.01,0.01), zlim=rev(c(0,0.35)),color=my_palette, angle=60,mgp=c(0.5,0.0,0.0),tck=0.02,
                                     mar=c(4,4,1,1),cex.lab=1.0)
                         sd3d$points3d(0.50,0.01,0.00,pch=16, col= "red",cex=2)
                         rm(sd3d)
           } else if (lnt>1) {
                  for (nset in 1:lnt)
                    {
                      temp             = objlist[[nrun]][nset]
                      temp             = unlist(temp)
                      temp[2]          = -temp[2]
                      temp[3]          = -temp[3]
                      temp[4]          = -temp[4]
                      temp[5]          = -temp[5]
                      tempar[nset,1]  = temp[1]
                      tempar[nset,2]  = temp[2]
                      tempar[nset,3]  = temp[3]
                      tempar[nset,4]  = temp[4]
                      tempar[nset,5]  = temp[5]
                    }
                   ind = sort(tempar[,3],decreasing=FALSE, index.return=TRUE)
                  sortedar[,1] = tempar[ind$ix,1]
                  sortedar[,2] = tempar[ind$ix,2]
                  sortedar[,3] = tempar[ind$ix,3]
                  sortedar[,4] = tempar[ind$ix,4]
                  sortedar[,5] = tempar[ind$ix,5]
          
                  #define colors
                   my_palette   = colorRampPalette(c("blue", "green", "red"))(n =lnt)
                       
                   xlim        = c(0, max(tempar[,2]))
                   sd3d        = scatterplot3d(sortedar[,2], sortedar[,4], sortedar[,1],            
                                              xlab='Expected utiltiy [-]', ylab='', zlab='Phosphorus in the lake [-]',
                                              xlim = c(0,0.45),ylim=c(-0.01,0.01),zlim=rev(c(0,0.35)),color=my_palette, angle=60,mgp=c(0.5,0.0,0.0),tck=0.02,
                                              mar=c(4,4,1,1),cex.lab=1.0)
                  sd3d$points3d(0.50,0.01,0.00,pch=16, col= "red",cex=2)
                  rm(sd3d)
           }
                  text(4.0,1,sprintf("NFE is %i",nfe[nrun]))
                  text(6.7,1.0,"Future gen [-]",srt=60,cex=1.0)
                  #colorbar
                  par("fg"="black")
                  color.legend(-0.5,7.7,0,9.0,legend=c("",""),
                               rect.col=my_palette,gradient="y",cex=1)
                  text(0.2,7.7,"min",srt=0,cex=0.7)
                  text(0.2,9.0,"max",srt=0,cex=0.7)
        
                  rm(temp, tempar, my_palette)
    }