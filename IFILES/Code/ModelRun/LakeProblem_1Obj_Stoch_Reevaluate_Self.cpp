/* LakeProblem_1Obj_Stoch_Reevaluate_Self.cpp
   
   Riddhi Singh, Feb, 2014
   The Pennsylvania State University
   rus197@psu.edu

   Reevlauating the solutions obtained by optimizing the single objective Bentham utiltiy with stochastic uncertainty.
   This reevaluation is against the assumed lognormal distributions, represented by 10000 states of the world (SOWs). 
   Once the objective functions are evaluated for each SOW, the average is taken to estimate the performance across all SOWs.
   This is following the 'Laplace's principle of insufficient reason'.
 
   Inputs 
   Parameters related to lake
      b : proportion of phosphorous retained in the lake each year    
      q : steepness of the sigmoid curve, large values give a steeper curve

   Parameters related to utility function
      delta   : discount rate
      alpha   : utility from pollution
      beta    : eutrophic cost

   State variables 
     lake_stateX : Phosphorous concentration at previous time step

   Decision Vector 
      vars : anthropogenic pollution flow at previous time step 
      100 years divided into 20 time steps, decision being made every 5 years

  Input Files
      1. SOWs_Type4.txt : a file containing all the SOWs against which the solutions are to be evaluated.
      2. LakeProblem_1ObjStochMM11April_Type4_Par_Pareto.txt: 
      File containing the decision vectors from the optmization. There are 10 vectors, from each random seed.
      
  Outputs
  
   Objectives estimated for different assumptions of uncertiainty:
   1. Phosphorous in the lake
   2. Bentham Utility
   3. Utility in first generation 
   4. Utility from 50-100 generations
   5. Reliability
   6. Inertia
   7. Bentham Utility with inertia

   Output Files
   1. LakeProblem_1ObjStochMM11April_Type4_Reev_Type4.txt -
   Contains the reevaluated objectives for each case. 
   Note that since we had 10 random seeds per optimization, but only the first seed was reevaluated. 
   Only the first random seed is used for figures in the paper. 
   The rest are for exploratory purposes only or constitute supplementary material.
*/

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <sstream>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <algorithm>                            //for sorting
using namespace std;

#define PI 3.14159265358979323846
#define nDays 100
#define q 2
#define b 0.42
#define alpha 0.4
#define beta 0.08
#define delta 0.98
#define samples 10000

//define reliability related parameters
#define pcrit 0.5
#define reliab_thres 0.90

//define inertia objective parameters
#define k1 0.004
#define k2 0.04
#define inertia_thres 0.01 

//initialize the stochastic term for the pollution flow
#define pol_sigma 0.00177827941
#define pol_mean  0.02

//initialize the stochastic term for the discount factor
#define dis_sigma 0.05
#define dis_mean 0.0

//decision space control and precision control parameters
int interval = 5;
int precis   = 3;

double passInitX= 0;

int nobjs;
int nvars;
double nat_flowmat [10000][nDays];

void lake_problem(double* vars, double* objs) 
{
  //run the stochastic case now for objective function 3
  double * ofs1 = new double [samples];    // Phosphrous in the lake
  double * ofs2 = new double [samples];    // Bentham utility - discounted sum of utilities
  double * ofs3 = new double [samples];    // utility of first generation
  double * ofs4 = new double [samples];    // utility of 50-100 generations
  double * ofs5 = new double [samples];    // reliability, Prob(Ph>Pcrit)
  double * ofs6 = new double [samples];    // Average inertia
  double * ofs7 = new double [samples];    // Bentham utilty with inertia

  for (int sample=0; sample<samples;sample++)
    {
      ofs1[sample]    = 0.0;
      ofs2[sample]    = 0.0;
      ofs3[sample]    = 0.0;
      ofs4[sample]    = 0.0;
      ofs5[sample]    = 0.0;
      ofs6[sample]    = 0.0;
      ofs7[sample]    = 0.0;
    }
   
  int linetouse [samples];
  srand (time(NULL));
  for (int sample=0; sample<samples;sample++) {
    //no random number - we are testing all SOWs here
    linetouse[sample] =sample;
  }
  for (int sample=0; sample<samples;sample++)
    {   
      double *nat_flow = new double [nDays];
      int index = linetouse[sample];
      // get the random natural flow from the States of the world file
      for (int i=0;i<nDays;i++){
	nat_flow[i] = nat_flowmat[index][i];
      }
      
      //initialize the random walk for discounted utility 
      boost::mt19937 rng(sample); 
      boost::normal_distribution<> nd2(dis_mean, dis_sigma);
      boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor2(rng, nd2);
            
      double * dis_val = new double [nDays];
      
      for (int i=1;i<nDays;i++)
	dis_val[i]  = var_nor2();
      
      int i;
      double *lake_stateX = new double [nDays];
      double * utility    = new double [nDays];
      double * poll       = new double [nDays];
      double * phos       = new double [nDays];
      double * pol_flow   = new double [nDays];
      double * dis_fac    = new double [nDays];
      double * usevars    = new double [nDays];
      double * inertia    = new double [nDays];
      double * utilityNoIn= new double [nDays];

      for (i=0; i < nDays; i++)
	{ lake_stateX[i] = 0.0;
	  utility[i]     = 0.0;
	  poll[i]        = 0.0;
	  phos[i]        = 0.0;
	  usevars[i]     = 0.0;
	  inertia[i]     = 0.0;
	}
                        
      //assigning 'nYears/interval' decision vectors to 'nYears' years
      for (i=0; i<nDays; i++)
	{
	  int pos = floor(i/interval);
	  usevars[i] = vars[pos];
	}   
      //now add natural pollution
      for (i=0; i < nDays; i++)
	{
	  dis_fac[i]     = delta + dis_val[i];
	  pol_flow[i]    = usevars[i]  + nat_flow[i];
	}
      
      //round off to the required precision
      for (i=0; i<nDays; i++)
	{
	  pol_flow[i] = round(pol_flow[i]*pow(10,(double)precis))/(pow(10,(double)precis));
	  usevars[i]  = round(usevars[i]*pow(10,(double)precis))/(pow(10,(double)precis));
	}
      
      //estimate inertia
      double tempin;
      for (i=0; i<nDays-1; i++)
	{
	  tempin = abs(usevars[i+1]-usevars[i]);
	  if (tempin<inertia_thres)
	    inertia[i] = tempin*k1;
	  else
	    inertia[i] = tempin*k2;
	}	  
      //set the last day to 0
      inertia[nDays-1] = 0.0;

      for (i=0; i<nDays; i++)
	{
	  if (i==0)
	    {
	      lake_stateX[i] = passInitX*(1-b)+pow(passInitX,q)/(1+pow(passInitX,q))+pol_flow[i];
	    }
	  else
	    {
	      lake_stateX[i] = lake_stateX[i-1]*(1-b)+(pow(lake_stateX[i-1], q))/(1+pow(lake_stateX[i-1],q))+pol_flow[i];
	    }

	  utility[i]     = alpha*usevars[i] - beta*pow(lake_stateX[i],2)-inertia[i];
	  utilityNoIn[i] = alpha*usevars[i] - beta*pow(lake_stateX[i],2);
	  phos[i]        = beta*pow(lake_stateX[i],2);
	  poll[i]        = alpha*usevars[i];
	}

      for (i=0; i<nDays; i++)
	{
	  ofs1[sample] = ofs1[sample] + lake_stateX[i];
	  ofs2[sample] = ofs2[sample] + utilityNoIn[i]*pow(delta,(i));
	  if (i==0)
	    ofs3[sample] = utilityNoIn[i];
	  if (i>49)
	    ofs4[sample] = ofs4[sample]+utilityNoIn[i];
	  //estimate the reliability matrix
	  if (lake_stateX[i]<pcrit)
	    ofs5[sample] = ofs5[sample]+1;

	  ofs6[sample] = ofs6[sample] + inertia[i];
	  ofs7[sample] = ofs7[sample] + utility[i]*pow(delta,(i));
	}
      
      ofs1[sample] = ofs1[sample]/nDays;
      ofs4[sample] = ofs4[sample]/(nDays-50);
      ofs6[sample] = ofs6[sample]/nDays;

      delete [] lake_stateX;
      delete [] utility;
      delete [] poll;
      delete [] phos;
      delete [] pol_flow;
      delete [] dis_fac;
      delete [] dis_val;
      delete [] nat_flow;
      delete [] usevars;
      delete [] inertia;
    }
  double dumofs1 = 0.0;
  double dumofs2 = 0.0;
  double dumofs3 = 0.0;
  double dumofs4 = 0.0;
  double dumofs5 = 0.0;
  double dumofs6 = 0.0;
  double dumofs7 = 0.0;

  for (int sample=0;sample<samples;sample++)
    {
      dumofs1  = dumofs1 + ofs1[sample];   //Phosphorous in the lake
      dumofs2  = dumofs2 + ofs2[sample];   //Benthams's utility
      dumofs3  = dumofs3 + ofs3[sample];   //Utiltiy of first generation
      dumofs4  = dumofs4 + ofs4[sample];   //Utility of 50-100 generations
      dumofs5  = dumofs5 + ofs5[sample];   //reliability 
      dumofs6  = dumofs6 + ofs6[sample];   //inertia 
      dumofs7  = dumofs7 + ofs7[sample];   //utility estimates with inertia
   }

  objs[0] = dumofs1/samples;
  objs[1] = dumofs2/samples;
  objs[2] = dumofs3/samples;
  objs[3] = dumofs4/samples;

  double reliability  = dumofs5/(nDays*samples);
  if (reliability>1)
    exit(EXIT_FAILURE);

  objs[4]  = reliability;
  
  objs[5]  = dumofs6/samples;
  objs[6]  = dumofs7/samples;

  objs[0]  = objs[0];      //want to minmize phosphorous in the lake
  objs[1]  = objs[1];      //want to maximize Bentham's utility
  objs[2]  = objs[2];      //want to maximize utility of first generation
  objs[3]  = objs[3];      //want to maximize utility of last 50 generations
  objs[4]  = objs[4];      //want to maximize reliability
  objs[5]  = objs[5];      //want to minimize inertia
  objs[6]  = objs[6];      //want to maximize Bentham's utiltiy with inertia
    
  //cout<<"This step acieved"<<" "<<"cutoff of "<<cutoff<<" first term is "<< dummy<<"\n";
  // cout<<" and the objective functions were OF1 is "<<objs[0]<<" "<<" OF2 is "<<objs[1]<<" OF3 is "<<objs[2]<<" OF4 is  "<<objs[3]<<" OF5 is "<<objs[4]<<"\n";
  delete [] ofs1;
  delete [] ofs2;
  delete [] ofs3;
  delete [] ofs4;
  delete [] ofs5;
  delete [] ofs6;
  delete [] ofs7;
}

int main(int argc, char* argv[]) 
{

 nvars = nDays/interval;
 nobjs = 7;
 double *vars = new double [nvars];
 double *objs = new double [nobjs];
 
 string path           = "/Users/rus197/Desktop/PostDoc2013/LakeProblem/Lake_Code/Carpenter1999/BORG/Borg-1.6Parallel/Output/April11Redo/";
 string SOWheader      = "SOWs_Type";
 string myheader       = "LakeProblem_1ObjStochMM11April_Type4";   //for single objective stochastic
 //string myheader       = "LakeProblem_1ObjDetMM11April_Type4";        // for single objective deterministic
 string myext          = "_Par_Pareto";
 string outext         = "_Reev";
 string ext            = ".txt";
 
 int nfiles            = 9;
 string sowfilename;
 string myfilename;
 string outfilename;
 
 for (int solfnum=4; solfnum<5; solfnum++)
   {
     
     //now open the output file and the file whose solutions are to be evaluated for these states of the world 
     string tempind;
     ostringstream convert;
     convert << solfnum;
     tempind = convert.str();
     
     sowfilename = path + SOWheader + tempind + ext;
     cout<<" Random numbers file is " << sowfilename<<endl;

     const char *sowfname = sowfilename.c_str();
     FILE *SOWfile;
     SOWfile   = fopen(sowfname,"r");
     
     int linenumSOW =0;
     int maxSizeSOW =5000;
     
     if (SOWfile==NULL) perror("Error opening file");
     else 
       {
	 char bufferSOW [maxSizeSOW];
	 while ( fgets(bufferSOW, maxSizeSOW, SOWfile)!=NULL) 
	   { linenumSOW++;
	     if (bufferSOW[0]!='#')
	       {
		 char *pEndSOW;
		 char *testbufferSOW = new char [maxSizeSOW];
		 for (int i=0; i <maxSizeSOW; i++)
		   testbufferSOW[i] = bufferSOW[i];
		 for (int cols =0;cols<nDays;cols++)              // use nDays for reading in the uncertain daily data
		   {
		     nat_flowmat[linenumSOW-1][cols] = strtod(testbufferSOW, &pEndSOW);
		     testbufferSOW  = pEndSOW;	
		   }				
	       }
	   }
       }
     fclose(SOWfile);
  
     int filer;
     int linenum = 0;
     int maxSize = 30000;
     
     myfilename  = path + myheader + myext + ext;
     const char *myfname = myfilename.c_str();
     FILE * myfile;
     myfile = fopen(myfname,"r");
 
  
     outfilename  = path + myheader + "Type" + tempind + outext + ext;
     const char *outfname = outfilename.c_str();
     FILE * outfile;
     outfile = fopen(outfname,"w");

     cout<<"Input file name is "<<myfilename << " and output file is " <<outfilename<<endl;
     double trial;
 
     if (myfile==NULL) perror("Error opening file");
     else 
       {
	 char buffer [maxSize];
	 while ( fgets(buffer, maxSize, myfile)!=NULL) 
	   { linenum++;
	     if ( (buffer[0]=='#') && linenum<13)
	       {
	     fprintf(outfile, "%s", buffer );
	       }
	     else if (buffer[0]!='#')
	       {
		 char *pEnd;
		 char *testbuffer = new char [maxSize];
		 for (int i=0; i <maxSize; i++)
		   testbuffer[i] = buffer[i];
		 for (int cols =0;cols<nvars;cols++)
		   {
		     vars[cols] = strtod(testbuffer, &pEnd);
		     testbuffer  = pEnd;
		     
		     //    fprintf(outfile,"%.17g ",vars[cols]);
		   }
		 //evaluate the lake with this vector
		 lake_problem(vars, objs);
		 //spit the output to the output file
		 for (int nobj =0; nobj< nobjs;nobj++)
		   {
		     fprintf(outfile,"%.17g ", objs[nobj]);
		   }
		 fprintf(outfile,"\n");
	       }
	   }
       }     
     fclose(myfile);
     fclose(outfile);
   }
}          