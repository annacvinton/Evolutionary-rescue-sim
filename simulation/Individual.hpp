//
//  Individual.hpp
//  New IIASA Model
//
//  Created by Anna Christina Vinton on 6/26/18.
//  Copyright © 2018 Anna Christina Vinton. All rights reserved.
//


#ifndef Individual_hpp
#define Individual_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>

#include "SpatialContainer.hpp"

class SpatialContainer;

extern double MUT_SD;
extern double DISP_SD;

class Individual {

public:
    ///PATCH VALUE
    int size;
    double patchValue(double x,double y);
   // double patchValue(double x);
    
    // variables
    double x;
    double y;
    double u;
    double pertvalue;
    double repvalue;
    double pertname;

    std::vector<double> vecOfRandomNums(double);
    
    int newx;
    
    //std::string Text;
    
    double envvalue;
    int id;
    
    double patch;
    double psd;
  
    int *pIteration;
    
    double *pPertvalue;
    double *pRepvalue;
   
    double *pPertname;
    
    
    double *pcurrenttime;
    
  
    
    // constructor
    Individual();
    
    //methods
    
    Individual operator+(Individual other_individual) const;
    
    double renew_environmentalvalue();
    
   
    int random();
    
    void speak();
    
  
    
    ///BIRTHRATE
    double birth(SpatialContainer &myPopulation); //& makes it not copy
    
    ///DEATHRATE
    double death(SpatialContainer &myPopulation);
    
    
    
    double birthvalue;
    double deathvalue;
    bool dirty;   // death rate needs recomputing

};

#endif /* Individual_hpp */
