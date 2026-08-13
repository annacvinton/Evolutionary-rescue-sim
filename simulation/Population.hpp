//
//  Population.hpp
//  New IIASA Model
//
//  Created by Anna Christina Vinton on 6/26/18.
//  Copyright © 2018 Anna Christina Vinton. All rights reserved.
//

#ifndef Population_hpp
#define Population_hpp
//#define MAXINDIVIDUALS 10000


#include <stdio.h>
#include <array>
#include <fstream>
#include "Individual.hpp"
#include "SpatialContainer.hpp"

using namespace std;

class Population {
public:
   //variables
    int iteration;
   
    
    double currenttime=0;
    double time_counter=0;
    double timestep=0;
    double pertvalue=0;
    double repvalue=0;
    double patch=0;
    double pertname=0;
   
    
    
    int numindividuals;
    
    double totalevent;
 ofstream fout; //define file object

    
    SpatialContainer myPopulation;
    
    Population();
    
      void NextEvent(int iteration);
    
   void AddIndividual(double x,double y,double u, double id,double patch);
    
};
#endif /* Population_hpp */
