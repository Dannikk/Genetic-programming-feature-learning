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
#include <time.h>
#include <cstdlib>
#include <cmath>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include <Windows.h>
#include <tchar.h>
#include <cstring>

using namespace cv;
using namespace std;


double partialModelError(struct parameters *params, struct chromosome *chromo, struct dataSet *data){

    int i, p, tmp;
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

        meanImage = 0;
        meanPred = 0;
        imageError = 0;
        for (p=0; p<numRes; p++) {
/*            tmp = i*numImages + p;
            cout << tmp << endl;
            cout << getDataSetSampleInputs(data, i*numImages + p)[0] << "; " << getDataSetSampleInputs(data, i*numImages + p)[1] << endl;*/
            executeChromosome(chromo, getDataSetSampleInputs(data, i*numRes + p));
            imagePred[p] = getChromosomeOutput(chromo,0);

            // mean calculation
            meanImage += getDataSetSampleOutput(data,i*numRes + p,0) / numRes;
            meanPred += imagePred[p] / numRes;
        }
/*        cout << "______ for " << i << " image" << endl;
        for (int t=0; t < 10; t++){
            for (int tt = 0; tt < 2; tt++){
                cout << getDataSetSampleInputs(data, i*numRes + t)[tt] << "; ";
            }
            cout << endl;
        }*/


        // calculation of variance of predicted image and correlation
        variancePred = 0;
        correlation = 0;
        for (p=0; p<numRes; p++) {
            variancePred += pow(imagePred[p] - meanPred, 2) / numRes;
            correlation +=  (imagePred[p] - meanPred)*
                            (getDataSetSampleOutput(data,i*numRes + p,0) - meanImage) / numRes;
        }

        if ((variancePred - 0) < pow(10, -15)){
            totalError = pow(10, 15);
            break;
        }
        // save A and B coefficients (linear scaling)
        b = correlation / variancePred;
        a = meanImage - b*meanPred;
        setA(chromo, i, a);
        setB(chromo, i, b);

        // MSE calculation for one image
        for (p=0; p < numRes; p++){
            imageError += pow(getDataSetSampleOutput(data,i*numRes + p,0) - (a + imagePred[p] * b), 2);
        }
        imageError /= numRes;

        totalError += imageError;
    }

    free(imagePred);


    return totalError / numImages;
}


struct dataSet* loadDataSetFromImages(char* sourseDir, int numImages, int width, int height, bool logging){
    vector<string> imageDirs = vector<string>();

    string source = R"(C:\Users\nikit\CLionProjects\GPFL\images2gpfl\)";
    char *src2find = strdup(R"(C:\Users\nikit\CLionProjects\GPFL\images2gpfl\*)");

    int numRes = width * height;

    double* outputs;
    double* inputs;
    struct dataSet* dataSet2return;

    Mat image;

    WIN32_FIND_DATA FindFileData;
    HANDLE hf;

    hf=FindFirstFile(src2find, &FindFileData);

    if (hf!=INVALID_HANDLE_VALUE){
        do{
            cout << "\t\t" << "FindFileData.cFileName: " << FindFileData.cFileName << endl;
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

    if (numImages > int(imageDirs.size())) {
        cout << "Warning!\n"
                "\tExpected number of images is " << numImages << ", but real number of images is " <<
                int(imageDirs.size()) << "\n";
        numImages = int(imageDirs.size());
    }

    outputs = new double[numImages * numRes];
    inputs = new double[2 * numImages * numRes];

    for (int i = 0; i < numImages; i++) {
        if (logging){
            cout << "Reading: " << imageDirs[i] << " ..." << endl;
        }
        image = imread(imageDirs[i], IMREAD_GRAYSCALE);

        uchar *arr = image.isContinuous() ? image.data : image.clone().data;

        for (int k = 0; k < numRes; k++) {
            inputs[2*i * numRes + 2 * k] = double(k / width);
            inputs[2*i * numRes + 2 * k + 1] = double(k % width);

            /*if (logging) {
                if (k >= 0 && k < 10) {
                    cout << 2 * i * numRes + 2 * k << ": " << inputs[2 * i * numRes + 2 * k] << "; "
                         << 2 * i * numRes + 2 * k + 1 << ": " << inputs[2 * i * numRes + 2 * k + 1] << endl;
                }
            }*/
//            cout << i * numRes + k << ": " << inputs[i * numRes + k] << "; " << i * numRes + k + 1 << ": " << inputs[i * numRes + k + 1] << endl;
            outputs[i * numRes + k] = double(arr[k]) / 255.0;
        }
        cout << endl;

    }

    /*int c = 0;
    double *testArr = new double[4194304];
    for (int i=0; i < 2000; i++){
        for (int j=0; j < 2000; j++){
//            cout << *outputs[2048*i + 2048*665 + 665 + j] << ", ";
//            cout << 2048*i + 2048*665 + 665 + j << "_";
            testArr[c] = outputs[2048*(i) + j];
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


    cv::Mat greyImg = cv::Mat(2000, 2000, CV_64F, testArr);
    greyImg *= 255.0;
//    std::memcpy(greyImg.data, &tarr, 7 * 7 * sizeof(uint8_t));
    imwrite("test.png", greyImg);
    greyImg *= 1 / 255.0;
    std::string greyArrWindow = "Grey Array Image";
    cv::namedWindow(greyArrWindow, cv::WINDOW_NORMAL);
    cv::imshow(greyArrWindow, greyImg);
    waitKey(0);*/


    dataSet2return = initialiseDataSetFromArrays(2, 1, numImages * numRes, inputs, outputs);
    cout << "Recording finished!" << endl;

    delete[] inputs;
    delete[] outputs;
    /*delete[] testArr;
    delete[] tarr;*/

    /*for (int i=0; i < numImages; i++){
        cout << "______ for " << i << " image" << endl;
        for (int t=0; t < 10; t++){
            for (int tt = 0; tt < 2; tt++){
                cout << i*numRes + t + tt << ": " << getDataSetSampleInputs(dataSet2return, i*numRes + t)[tt] << "; ";
            }
            cout << endl;
        }
    }*/

    return dataSet2return;
}

void save_pictures(struct dataSet* data, vector<struct chromosome*> chromos, const struct parameters* params) {
    int width = getWidth(params);
    int height = getHeight(params);

    double* imageArray = new double[width*height];
    double res, Sum=0;
    string name;
    int numRes = getImageResolution(params);
    int numImages = getNumImages(params);

    for (int i = 0; i < numImages; i++){
        Sum = 0;
        for (int p=0; p < numRes; p++) {
            double pixel = 0;
            for (auto chrm: chromos) {
                executeChromosome(chrm, getDataSetSampleInputs(data, i * numRes + p));
                res = getChromosomeOutput(chrm, 0);
                pixel += getA(chrm, i) + getB(chrm, i) * res;
            }
            imageArray[p] = pixel;
            Sum += pixel;
        }
        cout << "A: " << getA(chromos[0], i) << " B: " << getB(chromos[0], i) << endl;
        cout << "Mean of (predicted) " << i << " image: " << Sum / numRes << endl;
        cv::Mat greyImg = cv::Mat(height, width, CV_64F, imageArray);
        greyImg *= 255.0;
        name = to_string(i) + "_.png";
        imwrite(name, greyImg);

        std::string greyArrWindow = "Grey Array Image";
        cv::namedWindow(greyArrWindow, cv::WINDOW_NORMAL);
        cv::imshow(greyArrWindow, greyImg / 255);
        waitKey(0);
    }

    delete[] imageArray;
}

struct dataSet* updateDataSet(dataSet* oldDataSet, chromosome* chromo, parameters* params, bool logging = false){
    int numImages = getNumImages(params);
    int numRes = getImageResolution(params);
    int width = getWidth(params);
    struct dataSet* dataSet2return;

    double* outputs = new double[numImages * numRes];
    double* inputs = new double[2 * numImages * numRes];

    for (int i = 0; i < numImages; i++) {
        double a = getA(chromo, i);
        double b = getB(chromo, i);

        for (int p = 0; p < numRes; p++) {
            inputs[2*i * numRes + 2 * p] = double(p / width);
            inputs[2*i * numRes + 2 * p + 1] = double(p % width);

            /*if (logging) {
                if (k >= 0 && k < 10) {
                    cout << 2 * i * numRes + 2 * k << ": " << inputs[2 * i * numRes + 2 * k] << "; "
                         << 2 * i * numRes + 2 * k + 1 << ": " << inputs[2 * i * numRes + 2 * k + 1] << endl;
                }
            }*/
//            cout << i * numRes + k << ": " << inputs[i * numRes + k] << "; " << i * numRes + k + 1 << ": " << inputs[i * numRes + k + 1] << endl;
            executeChromosome(chromo, getDataSetSampleInputs(oldDataSet, i*numRes + p));
            outputs[i * numRes + p] = getDataSetSampleOutput(oldDataSet, i*numRes + p,0)
                    - (a + b*getChromosomeOutput(chromo,0));
        }
    }
    freeDataSet(oldDataSet);
    dataSet2return = initialiseDataSetFromArrays(2, 1, numImages * numRes, inputs, outputs);

    delete[] inputs;
    delete[] outputs;

    return dataSet2return;
}

int learn_features(int max_int_iter, int ext_iter, bool logging = false) {
    char *src2find = strdup(R"(C:\Users\nikit\CLionProjects\GPFL\images2gpfl\)");
    const int numImages = 3;
    const int width = 512;
    const int height = 512;

    struct parameters* params = nullptr;
    struct dataSet* trainingData = NULL;
    vector<struct chromosome*> chromos;
    struct chromosome** models = NULL;
    struct chromosome* new_chromo;

    int numInputs = 2;
    int numNodes = 100;
    int numOutputs = 1;
    int nodeArity = 2;

    int numThreads = 4;
    int numGens = max_int_iter;
    double targetFitness = 0.01;
    int updateFrequency = 5;

    time_t timeStart, timeEnd;
    double runningTime;

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);
    addNodeFunction(params, "add,sub,mul,div,sin,pow,exp,1,sig,tanh");

    //    important detail for GPFL
    setCustomFitnessFunction(params, partialModelError, "partialModelError");

    setNumThreads(params, numThreads);
    setTargetFitness(params, targetFitness);
    setUpdateFrequency(params, updateFrequency);

//    important detail for GPFL
    setNumImages(params, numImages);
    setImageResolution(params, width*height);
    setWidth(params, width);
    setHeight(params, height);

    if (logging)
        cout << "- Dataset loading from images starts:" << endl;

    trainingData = loadDataSetFromImages(src2find, numImages, width, height, true);

    if (logging)
        cout << "- Dataset loading from images completed!" << endl;

    for (int m = 0; m < ext_iter; m++) {
        timeStart = time(NULL);
        new_chromo = runCGP(params, trainingData, numGens);
        timeEnd = time(NULL);
        runningTime = difftime(timeEnd, timeStart);
        if (logging)
            cout << "Time: " << runningTime << endl;

        chromos.push_back(new_chromo);
        if (logging)
            cout << "Chromosome was added" << endl;
        save_pictures(trainingData, chromos, params);
        trainingData = updateDataSet(trainingData, new_chromo, params);
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


    learn_features(20, 7, true);

    return 0;
}
