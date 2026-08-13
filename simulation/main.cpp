//convo with dave, find s shape funciton at middle of range
//  main.cpp
//  New IIASA Model
//
//  Created by Anna Christina Vinton on 6/26/18.
//  Copyright © 2018 Anna Christina Vinton. All rights reserved.
//
//need to record and do stats first...after t=250, min(n),mean(n),mean(u),var(u). How do the mean traits look at each location on the x axis for only 1 rep (illulstrative), show distribution at end..
#include "Population.hpp"
#include "Individual.hpp"
#include "SpatialContainer.hpp"
#include <iostream>
#include <time.h>
#include <fstream>
#include <string>
#include <array>
#include <list>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include <cstdlib>
#include <stdlib.h> //abs
#include <math.h>
#include <vector>
#include <numeric>
#include <string>
#include <functional>
#include <random>
#include<chrono> //record run time
#include<math.h>
#include <algorithm>
#include <cstdlib>
#include <iomanip>
///time to add patchiness term
using namespace std;


//ofstream fout;
///LANDSCAPE VALUE

//double sigma = 1.0;
//double r, s = 2.0 * sigma * sigma;
 // double sum = 0.0;
double S=1.2; //set from env SLOPE at start of main
//SKIP RUN_15
double C=0; //0,.2,.4,.6,.8,1 curvature
int CEN=0; //0 if s=/=0 otherise 50


double environmentalvalue(double x, double y,double pertvalue,double repvalue,double patch){
    return CEN + S*x + (1/2)*C*x*x + pertvalue + patch;
   // return 50 + pertvalue;
}


static double envd(const char* k,double d){const char*v=getenv(k);return v?atof(v):d;}
static int envi(const char* k,int d){const char*v=getenv(k);return v?atoi(v):d;}

int main() {
    S   = envd("SLOPE", 1.2);
    MUT_SD = envd("MUTSD", 0.0);
    DISP_SD = envd("DISP", 1.5);
    int TPERT = envi("TPERT", 450);
    int PERT   = envi("PERT", 12);
    int NREPS  = envi("NREPS", 1);
    int TMAX   = envi("TMAX", 150);
    int NFOUND = envi("NFOUND", 700);
    int SEED   = envi("SEED", 54321);
    int MAXPERT=17;
 //perreps from 0 to 8, 0,4,and 8.
    for(int pertrep=PERT;pertrep<PERT+1;pertrep+=1){
    int MAXREPS=NREPS;
       
    for(int rep=0; rep<MAXREPS; rep++){ //
 //get the start time
 // auto start = std::chrono::steady_clock::now();

    int iteration;

    ///MAIN LOOP-INITIALIZE POPULATION, CHOOSE NEXT EVENT
  
        int MAXTIME = TMAX;

    //initialize the population
    Population thePopulation;
    Individual theIndividual;
    ////
        int iterationbegin;
        iterationbegin=NFOUND;
    for(int i=0; i<iterationbegin; i++) {
   // double  x = rand() % 100 + 1;
    static mt19937 rnd(SEED+7);
    uniform_real_distribution<> d(0, 98);
   double  x = d(rnd);
   double y = d(rnd);
   double patch = 0;
//    double patch = theIndividual.patchValue(x,y);    
 //  double u = d(rnd);
   
///this is integer but pick new random function that produces double!!
        //thePopulation.AddIndividual(x,y,u,i);
        thePopulation.AddIndividual(x,y,CEN+S*x + (1/2)*C*x*x,i,patch); //start perfectly adapted
    }
    
    srand48(SEED+rep*1000);
    iteration=iterationbegin;
    while(thePopulation.currenttime<MAXTIME){
       
            if(thePopulation.currenttime<TPERT){
            thePopulation.pertvalue=0;
            thePopulation.repvalue=rep;
           thePopulation.pertname=pertrep;
        }else{
            thePopulation.pertvalue=pertrep;
            thePopulation.pertname=pertrep;
        };
        thePopulation.iteration = iteration;
        thePopulation.NextEvent(iteration);
        thePopulation.repvalue=rep;
       
        
  /*      if(thePopulation.numindividuals==0 || thePopulation.time_counter==500||thePopulation.time_counter==550||thePopulation.time_counter==600||thePopulation.time_counter==650||thePopulation.time_counter==700||thePopulation.time_counter==750||thePopulation.time_counter==800||thePopulation.time_counter==850||thePopulation.time_counter==900||thePopulation.time_counter==950||thePopulation.time_counter==1000||thePopulation.time_counter==1050){
            //  AT SOME POINT NEED TO TAKE MEAN/SKEWNESS OF TRAITS..
        //    fout.open("newtry.rtf",ios::app);
            ofstream fout;
            fout.open ("test.txt");
fout<<rep<<","<<pertrep<<","<<thePopulation.time_counter<<","<<thePopulation.numindividuals<<endl;
            fout.close();
        }
  */
       /*
        if(thePopulation.numindividuals==0 || thePopulation.currenttime>=MAXTIME){
       cout<<"persistence time:"<<thePopulation.currenttime<<","<<thePopulation.time_counter<<","<<"rep:"<<rep<<","<<thePopulation.numindividuals<<endl;
        }
   */
     
        
         if(thePopulation.numindividuals==0) break;
         if(thePopulation.currenttime>MAXTIME) break;
        
        iteration++;
    }
  //  auto end = std::chrono::steady_clock::now();
  //  double elapsed_time= double (std::chrono::duration_cast <std::chrono::seconds> (end-start).count());
  // cout<<"elapsed time (sec):"<<elapsed_time<<endl;
    }
   
} //added bracket for reps
} //bracket for pertrep
