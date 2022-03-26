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

//#include <stdio>
//#include <Windows.h>
//#include <stdlib>
//#include <char.h>
//#include <string>
//#include <math>
//#include <ctime>
#include "cgp.h"
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include "Image4CGP.h"

using namespace cv;
using namespace std;


double partialModelError(struct parameters *params, struct chromosome *chromo, struct dataSet *data){

    int i, p;
    double totalError = 0;
    double imageError;

    double meanImage;
    double meanPred;

    double variancePred = 0;
    double correlation;

    double a, b;

    double *imagePred;

    const int numImages = getNumImages(params);
    const int numRes = getImageResolution(params);

    const double MAX_VALUE = pow(10, 10);
    const double EPSILON = pow(10, -15);

    double toReturn;

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

    meanPred = 0;
    //TODO: compare this
    for (p=0; p<numRes; p++) {
        executeChromosome(chromo, getDataSetSampleInputs(data, p));
        imagePred[p] = getChromosomeOutput(chromo,0);

        // mean calculation
        meanPred += imagePred[p] / numRes;
    }
    // and embedded mean function
    // ...

    // and also predict variance and embedded variance function
    for (p=0; p<numRes; p++){
        double diff = imagePred[p] - meanPred;
        variancePred += diff * diff / numRes;
    }

    for (i=0; i<numImages; i++) {
        meanImage = 0;
        imageError = 0;

        //TODO: compare this
        for (p=0; p<numRes; p++) {
            meanImage += getDataSetSampleOutput(data, i*numRes + p, 0) / numRes;
        }
        // and embedded mean function
        // ...

        // calculation of correlation
        correlation = 0;
        for (p=0; p<numRes; p++) {
            correlation +=  (imagePred[p] - meanPred) * (getDataSetSampleOutput(data,i*numRes + p,0) - meanImage) / numRes;
        }

        if (!isfinite(variancePred)){
//            cout << "Oops!: " << variancePred << endl;
//            cout << "meanPred: " << meanPred << endl;
            totalError = MAX_VALUE;
            setA(chromo, i, EPSILON);
            setB(chromo, i, EPSILON);
            continue;
        }

        if (variancePred < EPSILON) {
            totalError = MAX_VALUE;
            setA(chromo, i, EPSILON);
            setB(chromo, i, EPSILON);
            continue;
        }
        // save A and B coefficients (linear scaling)
        b = correlation / variancePred;
        a = meanImage - b*meanPred;
        setA(chromo, i, a);
        setB(chromo, i, b);

        // MSE calculation for one image
        for (p=0; p < numRes; p++){
            double part_error = getDataSetSampleOutput(data,i*numRes + p,0) - (a + imagePred[p] * b);
            imageError += part_error * part_error;
        }
        imageError /= numRes;

        totalError += imageError;
    }

    free(imagePred);

    toReturn = totalError / numImages;

    if (!isfinite(toReturn)) {
        cout << "NAN: \n\ttotalError: " << totalError << "; numImages: " << numImages << endl;
        cout << "\tvariancePred: " << variancePred << endl;
    }


    return toReturn;
}

void save_pictures(const struct dataSet* data, const vector<struct chromosome*>& chromos, const struct parameters* params,
        bool show_new_pictures = false) {

    int width = getWidth(params);
    int height = getHeight(params);

    auto* imageArray = new double[width*height];
    double res;
    string name;
    int numRes = getImageResolution(params);
    int numImages = getNumImages(params);

    for (int i = 0; i < numImages; i++){
        //Sum = 0;
        for (int p=0; p < numRes; p++) {
            double pixel = 0;
            for (auto chrm: chromos) {
                executeChromosome(chrm, getDataSetSampleInputs(data, i * numRes + p));
                res = getChromosomeOutput(chrm, 0);
                pixel += getA(chrm, i) + getB(chrm, i) * res;
            }
            imageArray[p] = pixel;
            //Sum += pixel;
        }
//        cout << "A: " << getA(chromos[0], i) << " B: " << getB(chromos[0], i) << endl;
//        cout << "Mean of (predicted) " << i << " image: " << Sum / numRes << endl;
        cv::Mat greyImg = cv::Mat(height, width, CV_64F, imageArray);
        greyImg *= 255.0;
        name = to_string(i) + "_.png";
        imwrite(name, greyImg);

        if (show_new_pictures) {
            std::string greyArrWindow = String("Image") + to_string(i);
            cv::namedWindow(greyArrWindow, cv::WINDOW_NORMAL);
            cv::imshow(greyArrWindow, greyImg / 255);
            waitKey(0);
        }
    }

    delete[] imageArray;
}


int learn_features(int max_int_iter, int ext_iter, bool logging = false) {
    //char *src2find = strdup(R"(C:\Users\nikit\CLionProjects\GPFL\images2gpfl\)");
    string path2srcimages = string("../images2gpfl");
    const int numImages = 3;
    const int width = 512;
    const int height = 512;

    struct parameters* params;
    struct dataSet* trainingData;
    vector<struct chromosome*> chromos;
    struct chromosome* new_chromo;

    int numInputs = 2;
    int numNodes = 100;
    int numOutputs = 1;
    int nodeArity = 2;

    int numThreads = 4;
    int numGens = max_int_iter;
    double targetFitness = 0.01;
    int updateFrequency = 20;

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

    cout << "Number of threads: " << numThreads << endl;

    if (logging)
        cout << "- Dataset loading from images starts:" << endl;

    trainingData = loadDataSetFromImages(path2srcimages, numImages, width, height, true);

    if (logging)
        cout << "- Dataset loading from images completed!" << endl;

    for (int m = 0; m < ext_iter; m++) {
        cout << "\nIteration " << m << endl;
        timeStart = time(nullptr);
        new_chromo = runCGP(params, trainingData, numGens);
        timeEnd = time(nullptr);
        runningTime = difftime(timeEnd, timeStart);
        if (logging)
            cout << "Time: " << runningTime << endl;

        chromos.push_back(new_chromo);
        //save_pictures(trainingData, chromos, params);
        trainingData = updateDataSet(trainingData, new_chromo, params);
    }
    save_pictures(trainingData, chromos, params, true);


    freeDataSet(trainingData);
    // TODO: create a chromosome clearing method
    freeParameters(params);
    return 0;
}


int main() {

    cout << "Hello" << endl;

    learn_features(41, 150, true);

    return 0;
}
