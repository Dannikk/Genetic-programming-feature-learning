/*
	This file is part of CGP-Library
	Copyright (c) Andrew James Turner 2014 (andrew.turner@york.ac.uk)

	CGP-Library is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as published
	by the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	CGP-Library is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with CGP-Library.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "cgp.h"
#include <time.h>



int learn_features(const int dimension, const int ext_iter) {
    struct parameters* params = NULL;
    struct dataSet* trainingData = NULL;
    struct chromosome* chromo = NULL;
    struct chromosome** models = NULL;

    int numInputs = 82;
    int numNodes = 100;
    int numOutputs = 19;
    int nodeArity = 2;

    int numThreads = 4;

    int numGens = 20000;
    double targetFitness = 10;
    int updateFrequency = 500;

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);
    addNodeFunction(params, "add,sub,mul,div,sin,pow,exp,1");
    setNumThreads(params, numThreads);

    setTargetFitness(params, targetFitness);

    setUpdateFrequency(params, updateFrequency);



    return 0;
}


//struct chromosome* runGPFL(struct parameters* params, struct dataSet* data, int numGens) {
//
//  int i;
//  int gen;
//
//  /* bestChromo found using runCGP */
//  struct chromosome* bestChromo;
//
//  /* arrays of the parents and children */
//  struct chromosome** parentChromos;
//  struct chromosome** childrenChromos;
//
//  /* storage for chromosomes used by selection scheme */
//  struct chromosome** candidateChromos;
//  int numCandidateChromos;
//
//   /* error checking */
//   if (numGens < 0) {
//	 printf("Error: %d generations is invalid. The number of generations must be >= 0.\n Terminating CGP-Library.\n", numGens);
//	 exit(0);
//   }
//
//   if (data != NULL && params->numInputs != data->numInputs) {
//	 printf("Error: The number of inputs specified in the dataSet (%d) does not match the number of inputs specified in the parameters (%d).\n", data->numInputs, params->numInputs);
//	 printf("Terminating CGP-Library.\n");
//	 exit(0);
//   }
//
//   if (data != NULL && params->numOutputs != data->numOutputs) {
//	 printf("Error: The number of outputs specified in the dataSet (%d) does not match the number of outputs specified in the parameters (%d).\n", data->numOutputs, params->numOutputs);
//	 printf("Terminating CGP-Library.\n");
//	 exit(0);
//   }
//
//   /* initialise parent chromosomes */
//   parentChromos = (struct chromosome**)malloc(params->mu * sizeof(struct chromosome*));
//
//   for (i = 0; i < params->mu; i++) {
// 	parentChromos[i] = initialiseChromosome(params);
//   }
//
//   /* initialise children chromosomes */
//   childrenChromos = (struct chromosome**)malloc(params->lambda * sizeof(struct chromosome*));
//
//   for (i = 0; i < params->lambda; i++) {
// 	childrenChromos[i] = initialiseChromosome(params);
//   }
//
//   /* intilise best chromosome */
//   bestChromo = initialiseChromosome(params);
//
//   /* determine the size of the Candidate Chromos based on the evolutionary Strategy */
//   if (params->evolutionaryStrategy == '+') {
// 	numCandidateChromos = params->mu + params->lambda;
//   }
//   else if (params->evolutionaryStrategy == ',') {
// 	numCandidateChromos = params->lambda;
//   }
//   else {
// 	printf("Error: the evolutionary strategy '%c' is not known.\nTerminating CGP-Library.\n", params->evolutionaryStrategy);
// 	exit(0);
//   }
//
//   /* initialise the candidateChromos */
//   candidateChromos = (struct chromosome**)malloc(numCandidateChromos * sizeof(struct chromosome*));
//
//   for (i = 0; i < numCandidateChromos; i++) {
// 	candidateChromos[i] = initialiseChromosome(params);
//   }
//
//   /* set fitness of the parents */
//   for (i = 0; i < params->mu; i++) {
// 	setChromosomeFitness(params, parentChromos[i], data);
//   }
//
//   /* show the user whats going on */
//   if (params->updateFrequency != 0) {
// 	printf("\n-- Starting CGP --\n\n");
// 	printf("Gen\tfitness\n");
//   }
//
//   /* for each generation */
//   for (gen = 0; gen < numGens; gen++) {
//
// 	/* set fitness of the children of the population */
// #pragma omp parallel for default(none), shared(params, childrenChromos,data), schedule(dynamic), num_threads(params->numThreads)
// 	for (i = 0; i < params->lambda; i++) {
// 	  setChromosomeFitness(params, childrenChromos[i], data);
// 	}
//
// 	/* get best chromosome */
// 	getBestChromosome(parentChromos, childrenChromos, params->mu, params->lambda, bestChromo);
//
// 	/* check termination conditions */
// 	if (getChromosomeFitness(bestChromo) <= params->targetFitness) {
//
// 	  if (params->updateFrequency != 0) {
// 		printf("%d\t%f - Solution Found\n", gen, bestChromo->fitness);
// 	  }
//
// 	  break;
// 	}
//
// 	/* display progress to the user at the update frequency specified */
// 	if (params->updateFrequency != 0 && (gen % params->updateFrequency == 0 || gen >= numGens - 1)) {
// 	  printf("%d\t%f\n", gen, bestChromo->fitness);
// 	}
//
// 	/*
// 		Set the chromosomes which will be used by the selection scheme
// 		dependant upon the evolutionary strategy. i.e. '+' all are used
// 		by the selection scheme, ',' only the children are.
// 	*/
// 	if (params->evolutionaryStrategy == '+') {
//
// 	  /*
// 		  Note: the children are placed before the parents to
// 		  ensure 'new blood' is always selected over old if the
// 		  fitness are equal.
// 	  */
//
// 	  for (i = 0; i < numCandidateChromos; i++) {
//
// 		if (i < params->lambda) {
// 		  copyChromosome(candidateChromos[i], childrenChromos[i]);
// 		}
// 		else {
// 		  copyChromosome(candidateChromos[i], parentChromos[i - params->lambda]);
// 		}
// 	  }
// 	}
// 	else if (params->evolutionaryStrategy == ',') {
//
// 	  for (i = 0; i < numCandidateChromos; i++) {
// 		copyChromosome(candidateChromos[i], childrenChromos[i]);
// 	  }
// 	}
//
// 	/* select the parents from the candidateChromos */
// 	params->selectionScheme(params, parentChromos, candidateChromos, params->mu, numCandidateChromos);
//
// 	/* create the children from the parents */
// 	params->reproductionScheme(params, parentChromos, childrenChromos, params->mu, params->lambda);
//   }
//
//   /* deal with formatting for displaying progress */
//   if (params->updateFrequency != 0) {
// 	printf("\n");
//   }
//
//   /* copy the best best chromosome */
//   bestChromo->generation = gen;
//   /*copyChromosome(chromo, bestChromo);*/
//
//   /* free parent chromosomes */
//   for (i = 0; i < params->mu; i++) {
// 	freeChromosome(parentChromos[i]);
//   }
//   free(parentChromos);
//
//   /* free children chromosomes */
//   for (i = 0; i < params->lambda; i++) {
// 	freeChromosome(childrenChromos[i]);
//   }
//   free(childrenChromos);
//
//   /* free the used chromosomes and population */
//   for (i = 0; i < numCandidateChromos; i++) {
// 	freeChromosome(candidateChromos[i]);
//   }
//   free(candidateChromos);
//
//   return bestChromo;
//}


int main(void) {

    time_t timeStart, timeEnd;
    double totalTime;
    struct parameters* params = NULL;
    struct dataSet* trainingData = NULL;
    struct chromosome* chromo = NULL;

    /*int numInputs = 1;
    int numNodes = 15;
    int numOutputs = 1;
    int nodeArity = 2;*/

    int numInputs = 82;
    int numNodes = 100;
    int numOutputs = 19;
    int nodeArity = 2;

    int numThreads = 4;

    int numGens = 2000;
    /*double targetFitness = 0.1;*/
    double targetFitness = 100;
    int updateFrequency = 500;

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

    addNodeFunction(params, "add,sub,mul,div,sin,pow,exp,1");
    setNumThreads(params, numThreads);

    setTargetFitness(params, targetFitness);

    setUpdateFrequency(params, updateFrequency);

    printParameters(params);

    //printf("|__|_%d", params->numInputs);

    timeStart = time(NULL);

    // Note: you may need to check this path such that it is relative to your executable
    /*trainingData = initialiseDataSetFromFile("../../dataSets/symbolic.data");*/
    trainingData = initialiseDataSetFromFile("../../dataSets/ProbenBenchmarks/soybean/soybean1.txt");

    chromo = runCGP(params, trainingData, numGens);

    timeEnd = time(NULL);
    totalTime = difftime(timeEnd, timeStart);

    printChromosome(chromo, 0);

    printf("%d thread time: %.f seconds\n", numThreads, totalTime);

    freeDataSet(trainingData);
    freeChromosome(chromo);
    freeParameters(params);

    return 0;
}
