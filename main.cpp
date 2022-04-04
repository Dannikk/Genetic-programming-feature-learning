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
#include "CustomNodeFunctions.h"

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
        vector<string> file_names, bool show_new_pictures = false, string prefix = string()) {

    int width = getWidth(params);
    int height = getHeight(params);

    auto* imageArray = new double[width*height];
    double res;
    string name;
    int numRes = getImageResolution(params);
    int numImages = getNumImages(params);

    if (!prefix.empty()) {
        prefix += string("_");
    }

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

        name = string("Result/") + prefix + file_names[i] + ".png";
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


int learn_features(char* func_set, int max_int_iter, int ext_iter, const int width, const int height,
                   const int numImages, const string& path2srcimages, bool logging = false, int savingImageFreq = 0,
                   bool saveCromosomes = true) {
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

    int numThreads = 8;
    int numGens = max_int_iter;
    double targetFitness = 0.000001;
    int updateFrequency = 200;

//    time_t timeStart, timeRunning, timeEnd;
    unsigned int timeStart, timeRunning, timeEnd;
    double runningTime, updatingTime, savingTime;
    double meanRunningTime = 0, meanUpdatingTime = 0;

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);
    //addNodeFunction(params, "add,sub,mul,div,sin,pow,exp,1,sig,tanh");
    //addNodeFunction(params, "add,sub,mul,div,sin,pow,exp,1,0,sq,sqrt,step,sig,tanh");
    //addNodeFunction(params, "add,sub,mul,div,sig,1,min,max,not,xor,and,or,wire,step");
    addNodeFunction(params, func_set);
    addCustomNodeFunction(params, min, "min", -1);
    addCustomNodeFunction(params, max, "max", -1);

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
    printParameters(params);
    cout << "Number of threads: " << numThreads << endl;

    if (logging)
        cout << "- Dataset loading from images starts:" << endl;

    trainingData = loadDataSetFromImages(path2srcimages, file_names, numImages, width, height, true);
    if (trainingData == nullptr) {
        return 1;
    }

    if (logging)
        cout << "- Dataset loading from images completed!" << endl;

    FILE *fp = fopen("coeffs.txt", "w");

    for (int i=0; i < numImages; i++) {
        fprintf(fp, "%s\t\t\t\t", file_names[i].c_str());
    }
    fprintf(fp, "\n");

    string file_name;
    char *fname;
    for (int ei = 0; ei < ext_iter; ei++) {
        cout << "\nIteration " << ei << endl;
        timeStart = clock();
        new_chromo = runCGP(params, trainingData, numGens);
        timeRunning = clock();
        runningTime = double(timeRunning - timeStart);
        meanRunningTime += runningTime;

        trainingData = updateDataSet(trainingData, new_chromo, params);
        removeInactiveNodes(new_chromo);
        if (saveCromosomes) {
            file_name = string("Chromosomes/") + "chromo" + to_string(ei+1) + ".txt";
            fname = strcpy((char*)malloc(file_name.length()+1), file_name.c_str());
            saveChromosome(new_chromo, fname);
            delete [] fname;

            for (int i=0; i < numImages; i++) {
                fprintf(fp, "%.15f %.15f\t", getA(new_chromo, i), getB(new_chromo, i));
            }
            fprintf(fp, "\n");
            file_name = string("Graphs/") + "graph" + to_string(ei+1) + ".dot";
            fname = strcpy((char*)malloc(file_name.length()+1), file_name.c_str());
            saveChromosomeDot(new_chromo, 1, fname);
            delete [] fname;
        }
        timeEnd = clock();
        chromos.push_back(new_chromo);
        updatingTime = double(timeEnd - timeRunning);
        meanUpdatingTime += updatingTime;

        if (savingImageFreq != 0 && (ei + 1) % savingImageFreq == 0) {
            timeStart = clock();
            save_pictures(trainingData, chromos, params, file_names, false, to_string(ei+1));
            timeEnd = clock();
            savingTime += double(timeEnd - timeStart);
        }
/*        if (logging)
            cout << "Time: " << runningTime << "; " << updatingTime << endl;*/
    }
    fclose(fp);

    if (savingImageFreq == 0 || ext_iter % savingImageFreq != 0) {
        cout << "Save final predicts" << endl;
        timeStart = clock();
        save_pictures(trainingData, chromos, params, file_names, false, to_string(ext_iter));
        timeEnd = clock();
        savingTime += double(timeEnd - timeStart);
    }

    if (logging){
        cout << "Mean time:" << endl;
        cout << "\tRunning\tUpdating" << endl;
        cout << "\t" << meanRunningTime / ext_iter / 1000 << "\t" << meanUpdatingTime / ext_iter / 1000 << endl;
    }
    cout << "Saving time: " << savingTime / 1000 << endl;

    freeDataSet(trainingData);
    // TODO: create a chromosome clearing method
    freeParameters(params);
    return 0;
}


int main() {
    char func_set[] = "add,sub,mul,div,pow,exp,sig,1,0,wire,sin,cos";
    //char func_set[] = "add,sub,mul,div,sq,sqrt,tanh,1,0,wire,sin,cos";
    const int ext_iter = 12;
    const int int_iter = 10;
    const int width = 256;
    const int height = 256;
    const int numImages = 5;
    const int savingImageFrequency = 5;
    int error;
    string path2srcimages = string("../images2gpfl_256");

    cout << "Number of external iterations: " << ext_iter << endl;
    cout << "Number of internal iterations: " << int_iter <<  endl;
    cout << "Frequency of image saving: " << savingImageFrequency << endl;

    time_t start_time, end_time;

    start_time = time(nullptr);

    error = learn_features(func_set, int_iter, ext_iter, width, height, numImages,
                           path2srcimages, true, savingImageFrequency, true);
    if (error!=0) {
        cout << "An error has occurred!" << endl;
    }

    end_time = time(nullptr);

    int hours = 0, minutes = 0, seconds = 0;
    int total = difftime(end_time, start_time);

    cout << "Total time: " << total << " seconds" << endl;
    cout << "\tor " << total / 3600 << " hours " << total % 3600 / 60 << " minutes " << total % 3600 << " seconds" << endl;

    return 0;
}
