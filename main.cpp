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

#include <cstdio>
#include "cgp.h"
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include <Windows.h>
#include <tchar.h>
#include <cstring>

using namespace cv;
using namespace std;


double partialModelError(struct parameters *params, struct chromosome *chromo, struct dataSet *data){

    int i, p;
    double totalError = 0;
    double imageError = 0;

    double meanImage = 0;
    double meanPred = 0;

    double variancePred = 0;
    double correlation = 0;

    double a, b;

    double *imagePred;

    const int numImages = getNumImages(params);
    const int numRes = getImageResolution(params);

    imagePred = (double*)malloc(numRes*sizeof(double));

    if(getNumChromosomeInputs(chromo) !=getNumDataSetInputs(data)){
        printf("Error: the number of chromosome inputs must match the number of inputs specified in the dataSet.\n");
        printf("Terminating.\n");
        // return
        exit(0);
    }

    if(getNumChromosomeOutputs(chromo) != getNumDataSetOutputs(data)){
        printf("Error: the number of chromosome outputs must match the number of outputs specified in the dataSet.\n");
        printf("Terminating.\n");
        exit(0);
    }

    for (i=0; i<numImages; i++) {
        // calculation of predicted
        for (p=0; p<numRes; p++)
            imagePred[p] = 0;
        for (p=0; p<numRes; p++) {
            executeChromosome(chromo, getDataSetSampleInputs(data, i*numImages + p));
            imagePred[p] = getChromosomeOutput(chromo,0);
        }

        // MSE calculation for one image
        imageError = 0;
        for (p=0; p<numRes; p++) {
            imageError += pow(getDataSetSampleOutput(data,i*numImages + p,0) - imagePred[p], 2);
        }

        // ls coefs calculation

        // mean calculation
        meanImage = 0;
        meanPred = 0;
        for (p=0; p<numRes; p++){
            meanImage += getDataSetSampleOutput(data,i*numImages + p,0);
            meanPred += imagePred[p];
        }

        // calculation of variance of predicted image and correlation
        variancePred = 0;
        correlation = 0;
        for (p=0; p<numRes; p++) {
            variancePred += pow(imagePred[p] - meanPred, 2);
            correlation +=  (imagePred[p] - meanPred)*
                            (getDataSetSampleOutput(data,i*numImages + p,0) - meanImage);
        }

        // save A and B coefficients (linear scaling)
        b = correlation / variancePred;
        a = meanImage - b*meanPred;
        setA(chromo, i, a);
        setB(chromo, i, b);

        totalError += imageError;
    }

    free(imagePred);


    return totalError / numImages;
}


struct dataSet* loadDataSetFromImages(char* sourseDir, int numImages, int width, int height, bool logging){
    vector<string> imageDirs = vector<string>();

    string source = R"(C:\Users\nikit\CLionProjects\opencv_test\images\)";
    char *src2find = strdup(R"(C:\Users\nikit\CLionProjects\opencv_test\images\*)");

    int resolution = width*height;

    double** outputs;
    double** inputs;

    Mat image;

    WIN32_FIND_DATA FindFileData;
    HANDLE hf;

    hf=FindFirstFile(src2find, &FindFileData);

    if (hf!=INVALID_HANDLE_VALUE){
        do{
            if (strcmp(FindFileData.cFileName, ".") != 0 &&
                strcmp(FindFileData.cFileName, "..") != 0)  {
                imageDirs.emplace_back(source + FindFileData.cFileName);
            }
        }
        while (FindNextFile(hf,&FindFileData)!=0);
        FindClose(hf);
    }

    if (logging){
        cout << "Images were found:" << endl;
        for(const string& dir: imageDirs)
            cout << "\t" << dir << endl;
    }

//    int width = 2048;
//    int height = 2048;

    if (numImages > int(imageDirs.size())) {
        cout << "Warning!\n"
                "\tExpected number of images is %d, but real number of images is %d" <<
                numImages << int(imageDirs.size()) << "\n";
        numImages = int(imageDirs.size());
    }

    outputs = new double*[numImages * resolution];
    inputs = new double*[numImages * resolution];

    for (int i = 0; i < numImages; i++) {
        if (logging){
            cout << "Reading: " << imageDirs[i] << " ..." << endl;
        }
        image = imread(imageDirs[i], IMREAD_GRAYSCALE);

        uchar *arr = image.isContinuous() ? image.data : image.clone().data;

        for (int k = 0; k < resolution; k++) {
            inputs[i * resolution + k] = new double[2]{double(k / width),
                                                       double(k % width)};
            outputs[i * resolution + k] = new double[1]{double(arr[k]) / 255.0};
        }
    }

    int d = 1362585;
    int c = 0;
    double *testArr = new double[4194404];
    for (int i=0; i < 1000; i++){
        for (int j=0; j < 1000; j++){
//            cout << *outputs[2048*i + 2048*665 + 665 + j] << ", ";
//            cout << 2048*i + 2048*665 + 665 + j << "_";
            testArr[c] = *outputs[2048*(i) + 2048*665 + 665 + j];
//            cout << testArr[c] << "_";
            c++;
        }
//        cout << endl;
    }

    double *tarr = new double[49]{1, 0, 0, 0, 0, 0, 0,
                                  0, 1, 0, 0, 0, 0, 0,
                                  0, 0, 1, 0, 0, 0, 0,
                                  0, 0, 0, 1, 0, 0, 0,
                                  0, 0, 0, 0, 1, 0, 0,
                                  0, 0, 0, 0, 0, 1, 0,
                                  0, 0, 0, 0, 0, 0, 1};

    cout << "recording finished" << endl;
    cv::Mat greyImg = cv::Mat(1000, 1000, CV_64F, testArr);
    greyImg *= 255.0;
//    std::memcpy(greyImg.data, &tarr, 7 * 7 * sizeof(uint8_t));
    imwrite("test.png", greyImg);
    greyImg *= 1 / 255.0;
    std::string greyArrWindow = "Grey Array Image";
    cv::namedWindow(greyArrWindow, cv::WINDOW_NORMAL);
    cv::imshow(greyArrWindow, greyImg);
    waitKey(0);


    return initialiseDataSetFromArrays(2, 1, numImages, inputs[0], outputs[0]);
}

int learn_features(const int dimension, const int ext_iter) {
    char *src2find = strdup(R"(C:\Users\nikit\CLionProjects\opencv_test\images\)");
    const int numImages = 1;
    const int width = 2048;
    const int height = 2048;

    struct parameters* params = nullptr;
    struct dataSet* trainingData = NULL;
    vector<struct chromosome*> chromos;
    struct chromosome** models = NULL;

    int numInputs = 2;
    int numNodes = 100;
    int numOutputs = 1;
    int nodeArity = 2;

    int numThreads = 1;
    int numGens = 2000;
    double targetFitness = 0.1;
    int updateFrequency = 500;

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);
    addNodeFunction(params, "add,sub,mul,div,sin,pow,exp,1");
    setNumThreads(params, numThreads);
    setTargetFitness(params, targetFitness);
    setUpdateFrequency(params, updateFrequency);

//    important detail for GPFL
    setNumImages(params, numImages);
    setImageResolution(params, width*height);

    cout << "Hello" << endl;

    trainingData = loadDataSetFromImages(src2find, numImages, width, height, true);

    for (int m = 0; m < ext_iter; m++) {
        chromos.push_back(runCGP(params, trainingData, numGens));



    }


    freeDataSet(trainingData);
//    freeChromosome(chromo);
    freeParameters(params);
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

    cout << "Hello" << endl;
    /*time_t timeStart, timeEnd;
    double totalTime;
    struct parameters* params = NULL;
    struct dataSet* trainingData = NULL;
    struct chromosome* chromo = NULL;*/

    /*int numInputs = 1;
    int numNodes = 15;
    int numOutputs = 1;
    int nodeArity = 2;*/

    /*int numInputs = 82;
    int numNodes = 100;
    int numOutputs = 19;
    int nodeArity = 2;

    int numThreads = 4;

    int numGens = 2000;
    double targetFitness = 0.1;
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
    trainingData = initialiseDataSetFromFile("../../dataSets/symbolic.data");
    trainingData = initialiseDataSetFromFile("../../dataSets/ProbenBenchmarks/soybean/soybean1.txt");

    chromo = runCGP(params, trainingData, numGens);

    timeEnd = time(NULL);
    totalTime = difftime(timeEnd, timeStart);

    printChromosome(chromo, 0);

    printf("%d thread time: %.f seconds\n", numThreads, totalTime);

    freeDataSet(trainingData);
    freeChromosome(chromo);
    freeParameters(params);*/

    learn_features(2, 2);

    return 0;
}
