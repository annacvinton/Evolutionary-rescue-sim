//
//  Individual.cpp
//  New IIASA Model
//
//  Created by Anna Christina Vinton on 6/26/18.
//  Copyright © 2018 Anna Christina Vinton. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include <random>
#include <chrono>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <numeric>
#include <string>
#include <functional>
#include <algorithm>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <cstdlib>
#include "Individual.hpp"
#include "SpatialContainer.hpp"
//#include "Population.hpp"

using namespace std;

class SpatialContainer;


double environmentalvalue(double x, double y,double pertvalue,double repvalue,double patch);

double Individual::renew_environmentalvalue(){
    envvalue=environmentalvalue(x, y, *pPertvalue,*pRepvalue, patch);
    return envvalue;
}

///PATCH FUNCTION


double Individual::patchValue(double x,double y){
    int firstx=round(x);
    //
    int newy=round(y);
   // if(newy>0){
       newx=firstx+newy*100;
   //   }
    //  else{
    //      newx=firstx;
    //  }
    //
   // string file_name = "/Users/annavinton/Dropbox/ac0var25.txt";
    const char* envp = getenv("LANDSCAPE");
    string file_name = envp ? string(envp) : string("landscape.txt");
    ///
    static vector<string> text;                 // cached; file never changes mid-run
    if (text.empty()) {
        ifstream input_stream(file_name);
        if (!input_stream) cerr << "Can't open input file!";
        string line;
        while (getline(input_stream, line)) text.push_back(line);
    }
if(newx==10000){newx=9999;};
    string Text = text[newx];
    double patch;
    if ( ! (istringstream(Text) >> patch) ) patch = 0;
  //  cout<<patch<<","<<x<<","<<newx<<endl;
    return patch;
    
}

///BIRTHRATE
double Individual::birth(SpatialContainer &myPopulation){
   
    birthvalue=0.54; //change from 0.5 to 0.52;
    return birthvalue;
}
///DEATHRATE

double Individual::death(SpatialContainer &myPopulation){
    vector<Individual*> theList;
    theList = myPopulation.GetNeighborhoodIndividuals(*this);
    
    double compxyu;
    int num1;
    
    vector<Individual*>::iterator it;
    double sum=0;
    int num=0;
    for ( it = theList.begin(); it != theList.end(); it++)
    {
        double difx;
        double dify;
        double difu;
        double sds = 1.3;
      //  int sdc = 5;
 
        // Access the object through iterator
        Individual &currentIndividual = **it;
        difx = this->x - currentIndividual.x;
        dify = this->y - currentIndividual.y;
        difu = this->u - currentIndividual.u;
                
        compxyu= 1/((exp(sqrt((difx*difx)+(dify*dify))))/(sds*sds));
        
        num1=1;
        //num = number of individuals in neighborhood
        num+=num1;
        sum+=compxyu;
       // cout<<num<<","<<sum<<endl;
    }
  
    double k0=10;
    double sdk=3; 
    //  renew_environmentalvalue(); //change 5 
   deathvalue = sum / (k0 * exp(-((u - envvalue)*(u - envvalue)) / (sdk * sdk)));
  
 //   cout<<num<<","<<sum<<","<<"death"<<deathvalue<<","<<u - envvalue<<endl;
    
    return deathvalue;
}
/////

Individual::Individual(){
    dirty=true;
    x=0;
    y=0;
    u=0;
    id=0;
    patch=0;
  //  renew_environmentalvalue(); //change 4
}
//change 2 below
double MUT_SD = 0.0;   // set from env MUTSD in main
double DISP_SD = 1.5;  // set from env DISP in main
default_random_engine de(12345); //seed
    normal_distribution<double> nd(0.0, 1.5);//dispersal sd //mean followed by stdiv
 normal_distribution<double> nd0(0.0, 1.5);
 normal_distribution<double> nd1(0.0, 1.5);//mutation sd

Individual Individual::operator+(Individual other_individual) const
{
    // (locals that shadowed the globals removed -- they reseeded the RNG on every birth)
    double randomnumx= nd(de)*(DISP_SD/1.5);
    double randomnumy=nd0(de)*(DISP_SD/1.5);
    double randomnumu = (MUT_SD>0.0) ? nd1(de)*(MUT_SD/1.5) : 0.0;   // MUT_SD=0 -> no mutation
    
    Individual offspring;
    double offspringnumx;
    double offspringnumy;
    ////check that the spatial container includes 100 and 0, below updated 
    offspringnumx = (x+other_individual.x)/2 + randomnumx;
    if(offspringnumx >=0 && offspringnumx <=99 ){
        offspring.x=offspringnumx;
    } else if(offspringnumx>99){
        offspring.x=99;
    }
    else if(offspringnumx<0){
        offspring.x=0;
    }
   offspringnumy= (y+other_individual.y)/2 + randomnumy;
    if(offspringnumy>=0 && offspringnumy <99 ){
        offspring.y=offspringnumy;
    } else if(offspringnumy>=99){
        offspring.y=offspringnumy-99;
    }
    else if(offspringnumy<0){
        offspring.y=offspringnumy+99;
    }
   
   offspring.u= (u+other_individual.u)/2 + randomnumu;
    
    offspring.id = *pIteration;
    offspring.pIteration=pIteration;
    
    offspring.pRepvalue=pRepvalue;
    
    offspring.pPertvalue=pPertvalue;
   // offspring.renew_environmentalvalue();
    
    offspring.pPertname=pPertname;
    offspring.patch=offspring.patchValue(offspring.x,offspring.y);   // FIX: set patch FIRST
    offspring.renew_environmentalvalue();                            // then compute envvalue
  //  offspring.renew_environmentalvalue(); //change 3   
    return offspring;
}

void Individual::speak(){
    cout <<x<<","<<y<<","<<u<<","<<id<<patch<<endl;
}

