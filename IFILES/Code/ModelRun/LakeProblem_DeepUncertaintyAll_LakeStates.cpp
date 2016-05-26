/* LakeProblem_DeepUncertaintyAll_LakeStates.cpp
   
   Riddhi Singh, Feb, 2014
   The Pennsylvania State University
   rus197@psu.edu

   Estimating the time series of phosphorus in the lake using 3 solutions from the Pareto approximate set.
   The Pareto set is obtained after 5 objective optimization under stochastic uncertainty. 
   The solutions are labelled 'robust (green)','low phosphorus (blue)', and 'max. util (red)'. 
   These names are based on the objective(s) they mainly maximize (focus on).

  
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

  Solution Number
  The solution number from the Pareto set file. 
  Since there are three solutions, three numbers are used one by one.

  Input Files
      1. SOWs_Type0.txt to SOWs_Type8.txt : Files containing the SOWs against which the solutions are to be evaluated.
      2. LakeProblem_5ObjStochMM11AprilSelectedRSeed_Type4_Par_Pareto.txt: 
      File containing the decision vectors from the optmization done before. 
      Selected decision vectors from this file are used to generate output files for each solutions

  Outputs
  
  Output Files
   Objectives estimated for different assumptions of uncertiainty:
   1. LakeProblem_5Obj_1Const_StochMM11Apr_Type4_lakestate_Type_green0.txt to LakeProblem_5Obj_1Const_StochMM11Apr_Type4_lakestate_Type_green8.txt
   2. LakeProblem_5Obj_1Const_StochMM11Apr_Type4_lakestate_Type_red0.txt to LakeProblem_5Obj_1Const_StochMM11Apr_Type4_lakestate_Type_red8.txt
   3. LakeProblem_5Obj_1Const_StochMM11Apr_Type4_lakestate_Type_blue0.txt to LakeProblem_5Obj_1Const_StochMM11Apr_Type4_lakestate_Type_blue8.txt

   Time series of phosphorous in the lake for each uncertainty distributin (there are 9) 
   Each distribution is saved in a single file, with 10000 samples per distribution, there are 9 files per chosen solution.
*/

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <sstream>
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

//decision space control and precision control parameters
int interval = 5;
int precis   = 3;

//define reliability related parameters
#define pcrit 0.5
#define reliab_thres 0.90

double passInitX= 0;

int nobjs;
int nvars;
double nat_flowmat [10000][nDays];
double alllake_objs [10000][7] = {0};

void lake_problem(double* vars) 
{
  //run the stochastic case now for objective function 
  
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
      
      int i;
      double *lake_stateX = new double [nDays];
      double * utility    = new double [nDays];
      double * poll       = new double [nDays];
      double * phos       = new double [nDays];
      double * pol_flow   = new double [nDays];
      double * usevars    = new double [nDays];
      double * inertia    = new double [nDays];
      
      for (i=0; i < nDays; i++)
	{ lake_stateX[i] = 0.0;
	  utility[i]     = 0.0;
	  poll[i]        = 0.0;
	  phos[i]        = 0.0;
	  usevars[i]     = 0.0;
	  inertia[i]     = 0.0;
	}

      //estimate inertia
      for (i=0; i<nDays/interval-1; i++)
	{
	  inertia[i] = sqrt(pow((vars[i+1]-vars[i])*100/vars[i],2));
	}	  
      //set the last day to 0
      inertia[nDays/interval] = 0.0;
      // find the maximum inertia 
      sort(inertia, inertia+nDays/interval);                           // sorting inertia in ascending order 
      alllake_objs[sample][6]  = inertia[nDays/interval-1];                       // assign the maximum inertia to the test function

      //assigning 'nYears/interval' decision vectors to 'nYears' years
      for (i=0; i<nDays; i++)
	{
	  int pos = floor(i/interval);
	  usevars[i] = vars[pos];
	}   
      //now add natural pollution
      for (i=0; i < nDays; i++)
	{
	  pol_flow[i]    = usevars[i]  + nat_flow[i];
	}
      
      //round off to the required precision
      for (i=0; i<nDays; i++)
	{
	  pol_flow[i] = round(pol_flow[i]*pow(10,(double)precis))/(pow(10,(double)precis));
	  usevars[i]  = round(usevars[i]*pow(10,(double)precis))/(pow(10,(double)precis));
	}      
      
      //estimate average input pollution
      for (i=0; i<nDays; i++)
	{
	        alllake_objs[sample][5] =       alllake_objs[sample][5] + usevars[i];
	}
            alllake_objs[sample][5] =       alllake_objs[sample][5]/nDays;

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

	  utility[i]     = alpha*usevars[i] - beta*pow(lake_stateX[i],2);
	  phos[i]        = beta*pow(lake_stateX[i],2);
	  poll[i]        = alpha*usevars[i];
	}

      for (i=0; i<nDays; i++)
	{
	  alllake_objs[sample][0] =       alllake_objs[sample][0] + lake_stateX[i];
          alllake_objs[sample][1] =       alllake_objs[sample][1] + utility[i]*pow(delta,(i));
	  if (i==0)
	    alllake_objs[sample][2] = utility[i];
	  if (i>49)
            alllake_objs[sample][3] =     alllake_objs[sample][3]+utility[i];
	  //estimate the reliability matrix
	  if (lake_stateX[i]<pcrit)
	    alllake_objs[sample][4] =   alllake_objs[sample][4]+1;
	}
      
      alllake_objs[sample][0] =       alllake_objs[sample][0]/nDays;
      alllake_objs[sample][3] =       alllake_objs[sample][3]/(nDays-50);
      alllake_objs[sample][4] =       alllake_objs[sample][4]/nDays;

      delete [] lake_stateX;
      delete [] utility;
      delete [] poll;
      delete [] phos;
      delete [] pol_flow;
      delete [] nat_flow;
      delete [] usevars;
      delete [] inertia;
    }
}

int main(int argc, char* argv[]) 
{
  nvars = nDays/interval;
  nobjs = 7;
  double *vars = new double [nvars];
  double *objs = new double [nobjs];

  string path           = "/Data/";
  string solheader      = "LakeProblem_5ObjStochMM11AprilSelectedRSeed_Type";
  string SOWheader      = "SOWs_Type";
  string solext         = "_Par_Pareto";
  string ext            = ".txt";
  string outheader      = "LakeProblem_5ObjStochMM11AprilSelectedRSeed_Type";
  string outext         = "_allSOWOBJS_Type";
  
  string SOWfilename;
  string outfilename;
  string solfilename;
  string tempind;

  int nfiles  = 9;
 for (int solfnum=0; solfnum<nfiles; solfnum++)
   {
     
     //now open the output file and the file whose solutions are to be evaluated for these states of the world 
     ostringstream convert;
     convert << solfnum;
     tempind = convert.str();
     
     SOWfilename = path + SOWheader + tempind + ext;
     cout<<" Random numbers file is " << SOWfilename<<endl;
     
     const char *sowfname = SOWfilename.c_str();
     FILE *SOWfile;
     SOWfile   = fopen(sowfname,"r");
     
     int linenumSOW =0;
     int maxSizeSOW =45000;
     
     if (SOWfile==NULL) { cout<< sowfname<<endl; perror("Error opening file");}
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
     //random states of the world file closed
     int maxSize = 30000;
     
     int solnum = 342;    //other solutions are: 30 (green), 92 (red), and 342 (blue)
     string sol = "_blue";
     //now open the output file and the file whose solutions are to be evaluated for these states of the world 
     outfilename = path + outheader + "4" + outext + sol + tempind + ext;
     solfilename = path + solheader + "4" + solext + ext;
     
     const char *outfname = outfilename.c_str();
     FILE *outfile;
     outfile   = fopen(outfname,"w");
     
     const char *solfname = solfilename.c_str();
     FILE *solfile;
     solfile   = fopen(solfname,"r");
     
     if (solfile==NULL) { cout<< solfname <<endl; perror("Error opening file");}
     else 
       {
	 int linenum = 0;
	 
	 char buffer [maxSize];
	 while ( fgets(buffer, maxSize, solfile)!=NULL) 
	   { linenum++;
	     if (linenum==solnum)
	       {
		 char *pEnd;
		 char *testbuffer = new char [maxSize];
		 for (int i=0; i <maxSize; i++)
		   testbuffer[i] = buffer[i];
		 for (int cols =0;cols<nvars;cols++)      // +1 for identifier
		   {
		     vars[cols] = strtod(testbuffer, &pEnd);
		     testbuffer  = pEnd;
		   }
		 //evaluate the lake with this vector
		 lake_problem(vars);
		 //cout<<"reached here"<<endl;
		 //cout<<"sample is "<<alllake_states[0][0]<<endl;
		 for (int samplenum =0; samplenum< samples;samplenum++)
		   {
		     for (int objnum=0; objnum<7 ; objnum++)
		       fprintf(outfile,"%.17g ", alllake_objs[samplenum][objnum]);
		     fprintf(outfile,"\n");
		   }
		 fclose(outfile);	  	      
	       }   		
	   }
       }
     fclose(solfile);	          
   }
}
