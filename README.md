# Evolutionary Rescue in Complex Landscapes

## SUMMARY

???


## RESEARCH QUESTIONS
1. Resistance to perturbation
-> "What makes populations resistant to perturbations?"

2. Recovery from perturbation
	a) "What drives speed of recovery of population sizes?"
	b) "What drives ability to recover (i.e. ratio of pre-/post-perturbation population size) of populations?"

3. Resilience of populations
 -> "What combination of biological and spatial combinations allow for survival of perturbations and to what level of perturbation magnitude?"

4. Lightning in a bottle / Non-linear effects on survival
-> "Are there combinations of biological or spatial parameters which allow survival of especially large perturbations?"

# PRIMARY CONTACT & COLLABORATORS
## PRIMARY CONTACT
Anna C Vinton, PhD
NSF Postdoctoral Fellow in Biology
Department of Zoology, Oxford University
She/her 
Website: AnnaVinton.com
Twitter: @TheoreticalAnna

## COLLABORATORS, order is alphabetic and does not represent final publication order
- David Vasseur, david.vasseur@yale.edu
- Erik Kusch, erik.kusch@bio.au.dk
- Ulf Dieckmann, dieckmann@iiasa.ac.at
- Åke Brännström, , ake.brannstrom@umu.se

# DATA
## REPOSITORY
Currently residing on a MEGA.nz drive shared between Anna Vinton and Erik Kusch

## DATA SETS
Filenames describe simulation runs - made up of (and separated by "_"):
- AC ... the number following this indicates autocorrelation in the landscape (0-4)
- DI ... the number following this indicates dispersal kernel (1.5, 2)
- MU ... the number following this indicates whether mutations happened (1) or not (0)
- SL ... the number following this indicates slope in the landscape (.8, 1, 1.2)
- VA ... the number following this indicates variance in the landscape (25, 50, 100)

## VARIABLES
- pert.value .... perturbation magnitude at time t [numeric]
- pert.name ..... perturbation maximum for simulation run (remains constant over time) [numeric]
- rep ........... replicate for treatment [numeric]
- t ............. time step at which samples were taken [numeric]
- n ............. population size at time t [numeric]
- x ............. x-coordinate of individual at time t [numeric]
- y ............. y-coordinate of individual at time t [numeric]
- u ............. trait value of individual (remains constant over time and affects fitness) [numeric]
- id ............ identifier for individuals (remains constant over time) [numeric]
- patch ......... patchiness of environment (remains constant over time, but not in space) [numeric]

# PREVIOUS WORK
Anna Vinton's PhD thesis (link to that available?)

# FUNDING
???
