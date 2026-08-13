//
//  SpatialContainer.hpp
//  New IIASA Model
//
//  Created by Anna Christina Vinton on 6/28/18.
//  Copyright © 2018 Anna Christina Vinton. All rights reserved.
//

#ifndef SpatialContainer_hpp
#define SpatialContainer_hpp

#include <stdlib.h>
#include <array>
#include <stdio.h>
#include <map>
#include <list>
#include <vector>
#include "Individual.hpp"
#define SIZEX 100
#define SIZEY 100
#define MAXINT 10
#define POTENTIALMAXX 100
#define POTENTIALMAXY 100



using namespace std;

class Individual;

class SpatialContainer {
public:
    int max_i;
    int max_j;
    map<int,Individual> myworld[POTENTIALMAXX][POTENTIALMAXY] ;
   
    SpatialContainer();
   
    //Methods
    void insertIndividual(Individual passIndividual);
    void eraseIndividual(Individual passIndividual);
    
    vector<Individual*> GetAllIndividuals();
  
  vector<Individual*> GetNeighborhoodIndividuals(Individual &passIndividual);
    

   
};


#endif /* SpatialContainer_hpp */

