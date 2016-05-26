/* LakeProblem_5Obj_1Const_StochRT.cpp
   
   Riddhi Singh, Jan, 2014
   The Pennsylvania State University
   rus197@psu.edu

   Running BORG parallel for the stochastic lake model extended from Carpenter et al., 1999 and 
   save the runtime dynamics for visualizing the evolution of the Pareto approximate surface
   
   Stochasticity is introduced by:
   1. Natural variability around anthropogenic pollution flow

   Inptus 
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
      vars : anthropogenic pollution flow at previous time step (this was aval in R and MATLB versions)

  Outputs
     Utility and discounted at a given time step
     utility : utiltiy at every time step
     npv_util: discounted utility - this is also the objective function
     
     Updated lake_stateX 

 Objectives
 1. Phosphorous in the lake
 2. Bentham Utility
 3. Utility in first generation
 4. Utility from 50-100 generations
 5. Reliability

Additional features:
1. Reduced decision space - 20 decision vectors repeating themselves for 5 years each
*/

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <sstream>
#include <time.h>
#include <mpi.h>
#include "borgmm.h"
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
#define samples 100

//define reliability related parameters
#define pcrit 0.5
#define reliab_thres 0.90

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
int nconsts;
double nat_flowmat [10000][nDays];

void lake_problem(double* vars, double* objs, double* consts) 
{
  //run the stochastic case now for objective function 3
  double * ofs1 = new double [samples];    // Phosphrous in the lake
  double * ofs2 = new double [samples];    // utility discounted
  double * ofs3 = new double [samples];    // utility of first generation
  double * ofs4 = new double [samples];    // utility of 50-100 generations
  double * ofs5 = new double [samples];    // reliability, Prob(Ph>Pcrit)
  
  for (int sample=0; sample<samples;sample++)
    {
      ofs1[sample]    = 0.0;
      ofs2[sample]    = 0.0;
      ofs3[sample]    = 0.0;
      ofs4[sample]    = 0.0;
      ofs5[sample]    = 0.0;
    }

  int linetouse [samples];
  srand (time(NULL));
  for (int sample=0; sample<samples;sample++) {
    //pick a random number based on time
    linetouse[sample] = rand() % 10001;
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

      for (i=0; i < nDays; i++)
	{ lake_stateX[i] = 0.0;
	  utility[i]     = 0.0;
	  poll[i]        = 0.0;
	  phos[i]        = 0.0;
	  usevars[i]     = 0.0;
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
	  
	  utility[i] = alpha*usevars[i] - beta*pow(lake_stateX[i],2);
	  phos[i]    = beta*pow(lake_stateX[i],2);
	  poll[i]    = alpha*usevars[i];	    
	}

      for (i=0; i<nDays; i++)
	{
	  ofs1[sample] = ofs1[sample] + lake_stateX[i];
	  ofs2[sample] = ofs2[sample] + utility[i]*pow(delta,(i));
	  if (i==0)
	    ofs3[sample] = utility[i];
	  if (i>49)
	    ofs4[sample] = ofs4[sample]+utility[i];
	  //estimate the reliability matrix
	  if (lake_stateX[i]<pcrit)
	    ofs5[sample] = ofs5[sample]+1;
	}
      
      ofs1[sample] = ofs1[sample]/nDays;
      ofs4[sample] = ofs4[sample]/(nDays-50);

      delete [] lake_stateX;
      delete [] utility;
      delete [] poll;
      delete [] phos;
      delete [] pol_flow;
      delete [] dis_fac;
      delete [] dis_val;
      delete [] nat_flow;
      delete [] usevars;
    }
  double dumofs1 = 0.0;
  double dumofs2 = 0.0;
  double dumofs3 = 0.0;
  double dumofs4 = 0.0;
  double dumofs5 = 0.0;

  for (int sample=0;sample<samples;sample++)
    {
      dumofs1  = dumofs1 + ofs1[sample];   //Phosphorous in the lake
      dumofs2  = dumofs2 + ofs2[sample];   //Benthams's utility
      dumofs3  = dumofs3 + ofs3[sample];   //Utility of first generation
      dumofs4  = dumofs4 + ofs4[sample];   //Utility of last 50 generations
      dumofs5  = dumofs5 + ofs5[sample];   //reliability estimator
   }

  objs[0] = dumofs1/samples;
  objs[1] = dumofs2/samples;
  objs[2] = dumofs3/samples;
  objs[3] = dumofs4/samples;

  double reliability  = dumofs5/(nDays*samples);
  if (reliability>1)
    exit(EXIT_FAILURE);

  objs[4]  = reliability;
  
  if (reliability>reliab_thres)
    consts[0]= 0.0;
  else
    consts[0] = reliab_thres-reliability;
  
  objs[0]    = objs[0];      //want to minmize phosphorous in the lake
  objs[1]    = -objs[1];     //want to maximize Benthams's utility
  objs[2]    = -objs[2];      //want to maximize utility of first generation
  objs[3]    = -objs[3];     //want to maximize utility of last 50 generations
  objs[4]    = -objs[4];     //want to maximize reliability

  delete [] ofs1;
  delete [] ofs2;
  delete [] ofs3;
  delete [] ofs4;
  delete [] ofs5;
}

int main(int argc, char* argv[]) 
{
  nvars   = nDays/interval;
  nobjs   = 5;
  nconsts = 1;

  for (int i=0;i<samples;i++)
    for (int j=0;j<nDays;j++)
      nat_flowmat[i][j] = 0.0;
  
  FILE * myfile;
  myfile = fopen("SOWs_Type4.txt","r");
	
  int linenum =0;
  int maxSize =5000;
  
  if (myfile==NULL) perror("Error opening file");
  else 
    {
      char buffer [maxSize];
      while ( fgets(buffer, maxSize, myfile)!=NULL) 
	{ linenum++;
	  if (buffer[0]!='#')
	    {
	      char *pEnd;
	      char *testbuffer = new char [maxSize];
	      for (int i=0; i <maxSize; i++)
		testbuffer[i] = buffer[i];
	      for (int cols =0;cols<nDays;cols++)                // use nDays not nvars, since now they are different
		{
		  nat_flowmat[linenum-1][cols] = strtod(testbuffer, &pEnd);
		  testbuffer  = pEnd;	
		}				
	    }
	}
    }
  fclose(myfile);
  
  int i, j;
  int rank;
  char runtime[256];
  char timing[256];  

  // All multi-master runs need to call startup, specify the number of
  // islands (one master per island), and set the runtime limits.
  BORG_Algorithm_ms_startup(&argc, &argv);
  BORG_Algorithm_ms_islands(8);
  BORG_Algorithm_ms_max_time(4);
  BORG_Algorithm_output_frequency(100);
  BORG_Algorithm_ms_max_evaluations(1000000);

  // Enable global Latin hypercube initialization to ensure each island
  // gets a well sampled distribution of solutions.
  BORG_Algorithm_ms_initialization(INITIALIZATION_LATIN_GLOBAL);
  
  // Define the problem.  Problems are defined the same way as the
  // serial example (see dtlz2_serial.c).
  BORG_Problem problem = BORG_Problem_create(nvars, nobjs, nconsts, lake_problem);
  
  for (j=0; j<nvars; j++) {
    BORG_Problem_set_bounds(problem, j, 0.0, 0.1);
  }
  
  /*  for (j=0; j<nobjs; j++) {
    BORG_Problem_set_epsilon(problem, j, 0.001);
    }*/
  //setting epsilons individually  
  BORG_Problem_set_epsilon(problem, 0, 0.01);   //Phosphorous in the lake
  BORG_Problem_set_epsilon(problem, 1, 0.01);   //Benthams's utility
  BORG_Problem_set_epsilon(problem, 2, 0.001);   //Utility of first generation
  BORG_Problem_set_epsilon(problem, 3, 0.001);   //Utility of last 50 generations
  BORG_Problem_set_epsilon(problem, 4, 0.01);   //Reliability estimator

  // Get the rank of this process.  The rank is used to ensure each
  // parallel process uses a different random seed.
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  

  // When running experiments, we want to run the algorithm multiple
  // times and average the results.
  for (i=0; i<1; i++) {
    // Save runtime dynamics to a file.  Each master node will
    // write its local runtime dynamics to this file.  The %%d
    // gets replaced by the index of the master.
  
    sprintf(runtime, "runtime5Obj1ConstStochMM_Type411AprilRD_%d_%%d.txt", i);
    BORG_Algorithm_output_runtime(runtime);
    
    // Seed the random number generator.
    BORG_Random_seed(37*i*(rank+1));
    
    // Run the multi-master Borg MOEA on the problem.
    BORG_Archive result = BORG_Algorithm_ms_run(problem);
		
    // Only the controller process will return a non-NULL result.
    // The controller aggregates all of the Pareto optimal
    // solutions generated by each master.  Then print the Pareto
    // optimal solutions to the screen.
    if (result != NULL) {
      BORG_Archive_print(result, stdout);
      BORG_Archive_destroy(result);
    }
  }
  
  // Shutdown the parallel processes and exit.
  BORG_Algorithm_ms_shutdown();
  BORG_Problem_destroy(problem);
	
  return EXIT_SUCCESS; 
}
