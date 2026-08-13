//
//  Population.cpp
//  New IIASA Model
//
//  Created by Anna Christina Vinton on 6/26/18.
//  Copyright © 2018 Anna Christina Vinton. All rights reserved.
//
#include <cstdlib>
#include <cmath>
#include <vector>
#include "Individual.hpp"
#include "Population.hpp"
#include <cmath>
#include <iostream>
#include <time.h>
#include <fstream>
#include <array>
#include <random>
#include <math.h>
#include <tgmath.h>  //round

using namespace std;

ifstream fin;
ofstream fout;


extern double S;   // gradient slope, defined in main.cpp

// central moments -> mean, sample sd, skewness, excess kurtosis
static void moments(const vector<double>&v,double&mn,double&sd,double&sk,double&ku){
    size_t n=v.size();
    if(n==0){mn=sd=sk=ku=NAN;return;}
    double s=0; for(double x:v) s+=x; mn=s/n;
    double m2=0,m3=0,m4=0;
    for(double x:v){double d=x-mn,d2=d*d; m2+=d2; m3+=d2*d; m4+=d2*d2;}
    m2/=n; m3/=n; m4/=n;
    sd = (n>1)? sqrt(m2*n/(n-1.0)) : NAN;
    sk = (m2>0)? m3/pow(m2,1.5) : NAN;
    ku = (m2>0)? m4/(m2*m2)-3.0 : NAN;
}

// nearest-neighbour distance, expanding ring search over the spatial grid
static double nnDist(SpatialContainer &sc, Individual *me){
    int bi=(int)floor(me->x/MAXINT), bj=(int)floor(me->y/MAXINT);
    if(bi>=sc.max_i) bi=sc.max_i-1; if(bj>=sc.max_j) bj=sc.max_j-1;
    double best=1e18; int found=0;
    for(int r=1;r<=sc.max_i && (!found || r<=2);r++){
        for(int i=bi-r;i<=bi+r;i++) for(int j=bj-r;j<=bj+r;j++){
            if(i<0||j<0||i>=sc.max_i||j>=sc.max_j) continue;
            for(map<int,Individual>::iterator it=sc.myworld[i][j].begin();
                it!=sc.myworld[i][j].end(); ++it){
                Individual *o=&(it->second);
                if(o==me) continue;
                double dx=o->x-me->x, dy=o->y-me->y, d=sqrt(dx*dx+dy*dy);
                if(d<best){best=d; found=1;}
            }
        }
    }
    return found? best : NAN;
}


// An individual q sums p into its competition iff their boxes are adjacent, so a
// birth or death at (px,py) invalidates exactly the 3x3 block of boxes around it.
static void markDirty(SpatialContainer &sc,double px,double py){
    int bi=(int)floor(px/MAXINT), bj=(int)floor(py/MAXINT);
    if(bi>=sc.max_i) bi=sc.max_i-1; if(bj>=sc.max_j) bj=sc.max_j-1;
    for(int i=bi-1;i<=bi+1;i++) for(int j=bj-1;j<=bj+1;j++){
        if(i<0||j<0||i>=sc.max_i||j>=sc.max_j) continue;
        for(map<int,Individual>::iterator it=sc.myworld[i][j].begin();
            it!=sc.myworld[i][j].end(); ++it) it->second.dirty=true;
    }
}

Population::Population(){
    numindividuals = 0;
    iteration = 0;
}


void Population::AddIndividual(double x,double y,double u,double id,double patch){
    
    Individual myIndividual;
    
    myIndividual.x = x;
    myIndividual.y = y;
    myIndividual.u = u;
    myIndividual.id = id;
    myIndividual.pIteration = &iteration;
    myIndividual.pertvalue=pertvalue;
    myIndividual.repvalue=repvalue;
    myIndividual.pcurrenttime=&currenttime;
    myIndividual.pPertvalue =&pertvalue;
    myIndividual.pRepvalue=&repvalue;
    myIndividual.pPertname=&pertname;
    myIndividual.patch= myIndividual.patchValue(myIndividual.x,myIndividual.y);  // FIX: patch FIRST
    myIndividual.renew_environmentalvalue();                                     // then envvalue
  
    
    myPopulation.insertIndividual(myIndividual);
    numindividuals++;
    
    }

void Population::NextEvent(int iteration){
    
     int i;
    
    //CALCULATE TOTAL EVENT RATE
    
    totalevent=0; //initialize event rate
   
    vector<Individual*> theList;
    theList = myPopulation.GetAllIndividuals();
    
    //Create an iterator of std::list
    vector<Individual*>::iterator it;
    
    int counter=0;
    double totalu=0;
    // Make iterate point to begining and incerement it one by one till it reaches the end of list.
    for (it = theList.begin(); it != theList.end(); it++)
    {
        counter++;
        // Access the object through iterator
        // Update public birtvalue and deathvalue of object
       
        (*it)->birth(myPopulation);
        if((*it)->dirty){ (*it)->death(myPopulation); (*it)->dirty=false; }
        totalevent += (*it)->birthvalue + (*it)->deathvalue;
        totalu+=(*it)->u;
    }
  //  cout<<numindividuals<<","<<counter<<","<<totalevent<<endl;
    
    double r2=drand48();
    double  timestep=(-log(r2)/totalevent);
        currenttime += timestep;
//cout<<currenttime<<","<<numindividuals<<","<<timestep<<endl;
    
    //SELECT EVENT TO OCCUR
    //DRAW RANDOM NUMBER BETWEEN 0 AND TOTAL EVENTS FROM UNIFORM DIST
    double myrandomnumber;
    myrandomnumber=drand48()*totalevent;
    
    double tempevent;
    tempevent = 0;
    i=0;
    
    // Make iterate point to begining and incerement it one by one till it reaches the end of list.

    it = theList.begin();
    do{
        //TOTAL EVENT IS NOW TEMP EVENT
        
        tempevent += (*it)->birthvalue+(*it)->deathvalue;
        
        //IS RANDOM NUMBER IS LESS THEN TEMP EVENT, MOVE TO NEXT INDIVIDUAL
        if(tempevent<myrandomnumber){it++;}
    }
    while(tempevent < myrandomnumber);
    
        //CALCULATE WHETHER BIRTH OR DEATH HAPPENS
    
    if(tempevent - (*it)->deathvalue > myrandomnumber){
      
   //     cout<<"birth"<<","<<*it->pIteration<<endl;
       
        
        {   Individual kid = **it + **it;
            myPopulation.insertIndividual(kid);
            markDirty(myPopulation, kid.x, kid.y);   }
        /*we defined + here in individual.cpp, now it is adding the same individual*/
       
        numindividuals ++;
       
        } else {
            {   double dx=(*it)->x, dy=(*it)->y;
                myPopulation.eraseIndividual(**it);
                markDirty(myPopulation, dx, dy);   }
            numindividuals--;
        //    cout<<"death"<<","<<*it->pIteration<<endl;
        }
    

   
  if(time_counter<currenttime){
        time_counter += 10; //instead of 1
     // cout<<repvalue<<","<<time_counter<<","<<numindividuals<<","<<pertvalue<<","<<endl;
      
   
          
          ///this is most recent print code
      
          if(time_counter>0){
      //  cout<<numindividuals<<","<<time_counter<<","<<pertvalue<<endl;
              vector<Individual*> cur = myPopulation.GetAllIndividuals();
              vector<double> uu, mal, nnd, xx;
              uu.reserve(cur.size()); mal.reserve(cur.size());
              nnd.reserve(cur.size()); xx.reserve(cur.size());
              for (size_t q=0;q<cur.size();q++){
                  Individual *ind = cur[q];
                  uu.push_back(ind->u);
                  xx.push_back(ind->x);
                  mal.push_back(fabs(S*ind->x + ind->patch - ind->u));
                  double d = nnDist(myPopulation, ind);
                  if(d==d) nnd.push_back(d);
              }
              double um,us,uk1,uk2, mm,ms,mk1,mk2, nm,ns,nk1,nk2, xm,xs,xk1,xk2;
              moments(uu,um,us,uk1,uk2);
              moments(mal,mm,ms,mk1,mk2);
              moments(nnd,nm,ns,nk1,nk2);
              moments(xx,xm,xs,xk1,xk2);
              {
                  ofstream fout;
                  const char* op=getenv("OUTFILE"); fout.open(op?op:"out.txt",ios::app);
                  fout<<pertvalue<<","<<pertname<<","<<repvalue<<","<<time_counter<<","
                      <<cur.size()<<","<<um<<","<<us<<","<<uk1<<","<<uk2<<","
                      <<mm<<","<<ms<<","<<mk1<<","<<mk2<<","
                      <<nm<<","<<ns<<","<<nk1<<","<<nk2<<","
                      <<xm<<","<<xs<<","<<xk1<<","<<xk2<<endl;
                  fout.close();
              }
              if(getenv("FULLOUT")){
                  ofstream gout;
                  const char* fp=getenv("FULLOUT"); gout.open(fp,ios::app);
                  for(size_t q=0;q<cur.size();q++)
                      gout<<pertvalue<<","<<pertname<<","<<repvalue<<","<<time_counter<<","
                          <<cur.size()<<","<<cur[q]->x<<","<<cur[q]->y<<","<<cur[q]->u<<","
                          <<cur[q]->id<<","<<cur[q]->patch<<endl;
                  gout.close();
              }
              if(false){
                  ofstream fout;
               //   cout <<pertvalue<<","<<repvalue<<","<<time_counter<<","<<numindividuals<<","<<it->x<<"," <<it->y<<","<<it->u<<","<<it->id<<","<<abs(it->x - it->u)<< endl;
                  
                  fout.close();
              }
       
          }
      ///////////end most recent print code  ///

  }
      /*for(int i=0; i+=50;){
           if(i==time_counter){
               
   cout <<repvalue<<","<<pertvalue<<","<<time_counter<<","<<numindividuals<<endl;
   fout.open("testing3.txt",ios::app);
    // Make iterate point to begining and incerement it one by one till it reaches the end of list.
  //  for (it = theList.begin(); it != theList.end(); it++)/    {
    // Access the object through iterator
        
   //   fout <<currenttime<<","<<it->x<<"," <<it->y<<","<<it->u<<","<<it->id<<","<<abs(it->x - it->u)<<","<<numindividuals <<","<<pertvalue<<","<<repvalue<< endl;
    
    fout <<repvalue<<","<<pertvalue<<","<<time_counter<<","<<numindividuals<<endl;
    
           //     fout <<currenttime<<","<<numindividuals<<","<<repvalue<< endl;
   //         fout <<currenttime<<","<<it->x<<"," <<it->u<< endl;
               }
    fout.close();
  } */

}


