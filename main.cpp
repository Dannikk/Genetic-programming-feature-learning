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
#include <ctime>
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
    double meanPred = 0;

    double meanImage2;
    double meanPred2 = 0;
    double meanImagePred; // the value Image*Prediction
    double varianceImage;

    double variancePred = 0;
    double covariance;

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

    //TODO: compare this
    for (p=0; p<numRes; p++) {
        executeChromosome(chromo, getDataSetSampleInputs(data, p));
        //imagePred[p] = getChromosomeOutput(chromo,0);
        double pixel = getChromosomeOutput(chromo,0);
        imagePred[p] = pixel;

        // mean calculation
        //meanPred += imagePred[p] / numRes;
        meanPred += pixel;
        meanPred2 += pixel * pixel;
    }
    meanPred /= numRes;
    meanPred2 /= numRes;
    variancePred = meanPred2 - meanPred * meanPred;
    // and embedded mean function
    // ...

    /*// and also predict variance and embedded variance function
    for (p=0; p<numRes; p++){
        double diff = imagePred[p] - meanPred;
        variancePred += diff * diff / numRes;
    }*/

    for (i=0; i<numImages; i++) {
        meanImage = 0;
        imageError = 0;
        meanImage2 = 0;
        meanImagePred = 0;

        //TODO: compare this
        for (p=0; p<numRes; p++) {
            double pixel = getDataSetSampleOutput(data, i*numRes + p, 0);
            meanImage += pixel;
            meanImage2 += pixel * pixel;
            meanImagePred += pixel * imagePred[p];
        }
        meanImage /= numRes;
        meanImage2 /= numRes;
        varianceImage /= meanImage2 - meanImage * meanImage;
        covariance = (meanImagePred - numRes * meanImage * meanPred) / (numRes - 1); // /varianceImage/variancePred;
        // and embedded mean function
        // ...

        /*// calculation of covariance
        covariance = 0;
        for (p=0; p<numRes; p++) {
            covariance +=  (imagePred[p] - meanPred) * (getDataSetSampleOutput(data,i*numRes + p,0) - meanImage) / numRes;
        }*/

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
        b = covariance / variancePred;
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
        vector<string> file_names, bool show_new_pictures = false) {

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
        name = "../Result/" + file_names[i] + ".png";
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


int learn_features(int max_int_iter, int ext_iter, const int width, const int height,
                   const int numImages, const string& path2srcimages, bool logging = false) {
    //char *src2find = strdup(R"(C:\Users\nikit\CLionProjects\GPFL\images2gpfl\)");

    struct parameters* params;
    struct dataSet* trainingData;
    vector<struct chromosome*> chromos;
    struct chromosome* new_chromo;
    vector<string> file_names;

    int numInputs = 2;
    int numNodes = 300;
    int numOutputs = 1;
    int nodeArity = 2;

    int numThreads = 12;
    int numGens = max_int_iter;
    double targetFitness = 0.000001;
    int updateFrequency = 200;

//    time_t timeStart, timeRunning, timeEnd;
    unsigned int timeStart, timeRunning, timeEnd;
    double runningTime, updatingTime;
    double meanRunningTime = 0, meanUpdatingTime = 0;

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);
    //addNodeFunction(params, "add,sub,mul,div,sin,pow,exp,1,sig,tanh");
    //addNodeFunction(params, "add,sub,mul,div,sin,pow,exp,1,0,sq,sqrt,step,sig,tanh");
    addNodeFunction(params, "add,sub,mul,div,sig,1,min,max,not,xor,and,or,wire,step");

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

    trainingData = loadDataSetFromImages(path2srcimages, file_names, numImages, width, height, true);
    if (trainingData == nullptr){
        return 1;
    }

    if (logging)
        cout << "- Dataset loading from images completed!" << endl;

    for (int m = 0; m < ext_iter; m++) {
        cout << "\nIteration " << m << endl;
//        timeStart = time(nullptr);
        timeStart = clock();
        new_chromo = runCGP(params, trainingData, numGens);
//        timeRunning = time(nullptr);
        timeRunning = clock();
        //runningTime = difftime(timeRunning, timeStart);
        runningTime = double(timeRunning - timeStart) / 1000;
        meanRunningTime += runningTime;

        chromos.push_back(new_chromo);
        //save_pictures(trainingData, chromos, params);
        trainingData = updateDataSet(trainingData, new_chromo, params);
//        timeEnd = time(nullptr);
        timeEnd = clock();
        //updatingTime = difftime(timeEnd, timeRunning);
        updatingTime = double(timeEnd - timeRunning) / 1000;
        meanUpdatingTime += updatingTime;
//        if (logging)
//            cout << "Time: " << runningTime << "; " << updatingTime << endl;
    }

    if (logging){
        cout << "Mean time:" << endl;
        cout << "Running\tUpdating" << endl;
        cout << meanRunningTime / ext_iter << "\t" << meanUpdatingTime / ext_iter << endl;
    }
    timeStart = clock();
    save_pictures(trainingData, chromos, params, file_names, false);
    timeEnd = clock();
    cout << "Saving time: " << double(timeEnd - timeStart) / 1000 << endl;

    freeDataSet(trainingData);
    // TODO: create a chromosome clearing method
    freeParameters(params);
    return 0;
}


int main() {
    const int ext_iter = 5;
    const int int_iter = 200;
    const int width = 256;
    const int height = 256;
    const int numImages = 5;
    int error;
    string path2srcimages = string("../images2gpfl_256");

    cout << "Number of external iterations: " << ext_iter << endl;
    cout << "Number of internal iteratons: " << int_iter <<  endl;

    time_t start_time, end_time;

    start_time = time(nullptr);

    error = learn_features(int_iter, ext_iter, width, height, numImages, path2srcimages, true);
    if (error!=0){
        cout << "An error has occurred!" << endl;
    }

    end_time = time(nullptr);

    cout << "Total time: " << difftime(end_time, start_time) << endl;    

    return 0;
}
