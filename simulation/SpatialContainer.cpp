//
//  SpatialContainer.cpp
//  New IIASA Model
//
//  Created by Anna Christina Vinton on 6/28/18.
//  Copyright © 2018 Anna Christina Vinton. All rights reserved.
//
#include <cmath>
#include <ctgmath>
#include <math.h>
#include <map>
#include "SpatialContainer.hpp"
#include <list>
using namespace std;

class Individual;

/* insert, delete, get all, and get nearby individuals
 */
SpatialContainer::SpatialContainer(){
    max_i= ceil(SIZEX/MAXINT);
    max_j= ceil(SIZEX/MAXINT);
    
}

void SpatialContainer::insertIndividual(Individual passIndividual){
    
    int  boxi=floor(passIndividual.x/MAXINT);
    int  boxj=floor(passIndividual.y/MAXINT);
    
    if (boxi==max_i) boxi = max_i-1;
    if (boxj==max_j) boxj = max_j-1;
   
    myworld[boxi][boxj].insert(pair<int,Individual>(passIndividual.id,passIndividual));
   
}


void SpatialContainer::eraseIndividual (Individual passIndividual){
    //passIndividual;
    int  boxi=floor(passIndividual.x/MAXINT);
    int  boxj=floor(passIndividual.y/MAXINT);
    
    if (boxi==max_i) boxi = max_i-1;
    if (boxj==max_j) boxj = max_j-1;
    
    map<int, Individual>::iterator it;
    
    it = myworld[boxi][boxj].find(passIndividual.id);
    
    myworld[boxi][boxj].erase(it);
    
}

// This is unnecessarily slow - would be better to implement an iterator approach
// similar to the standard containers (see how it is done for map below).
vector<Individual*> SpatialContainer::GetAllIndividuals(){
  vector<Individual*>  myList;
    myList.reserve(4096);
    
    for (int i=0;i<max_i;i++){
        
        for (int j=0;j<max_j;j++){
            // Create a map iterator and point to beginning of map
            map<int, Individual>::iterator it = myworld[i][j].begin();
            
            // Iterate over the map using Iterator till end.
            while (it != myworld[i][j].end())
            {
                // Accessing KEY from element pointed by it.
               // string word = it->first;
                
                // Accessing VALUE from element pointed by it.
               myList.push_back(&(it->second));
                
                // Increment the Iterator to point to next entry
                it++;
            }
    }
}
    return myList;
}



///GET NEIGHBORHOOD OF INDIVIDUALS INSTEAD OF GET ALL INDIVIDUALS

vector<Individual*> SpatialContainer::GetNeighborhoodIndividuals(Individual &passIndividual){
    vector<Individual*>  myList;
    myList.reserve(256);
    int  boxi=floor(passIndividual.x/MAXINT);
    int  boxj=floor(passIndividual.y/MAXINT);
    //if on border, getting an incorrect neighborhood
    //GET NEIGHBORHOOD, LOOP THROUGH NEIGHBORHOOD
    if (boxi==max_i) boxi = max_i-1;
    if (boxj==max_j) boxj = max_j-1;
    

    int lower_value_i = boxi-1;
    int upper_value_i = boxi+1;
    int lower_value_j = boxj-1;
    int upper_value_j = boxj+1;
    
    // Need to correct for edge effects
    if (lower_value_i < 0) lower_value_i=0;
    if (upper_value_i == max_i) upper_value_i = max_i-1;
    if (lower_value_j < 0) lower_value_j = 0;
    if (upper_value_j == max_j) upper_value_j = max_j-1;
    
    
    for (int i=lower_value_i;i<upper_value_i+1;i++){
        
        for (int j=lower_value_j;j<upper_value_j+1;j++){
            // Create a map iterator and point to beginning of map
            map<int, Individual>::iterator it = myworld[i][j].begin();
            
            // Iterate over the map using Iterator till end.
            while (it != myworld[i][j].end())
            {
                // Accessing KEY from element pointed by it.
                // string word = it->first;
                
                // Accessing VALUE from element pointed by it.
                myList.push_back(&(it->second));
                
                // Increment the Iterator to point to next entry
                it++;
            }
        }
    }
    return myList;
}

///WEIGHTED LOCAL POP SIZE- ITERATE OVER INDIVIDUALS IN NEIGHBORHOOD AND SUM COMPETITIVE INFLUENCE

