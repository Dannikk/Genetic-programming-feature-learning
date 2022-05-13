//
// Created by DNA on 25.3.22.
//

#include "Image4CGP.h"
#include <string>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
//#include <Windows.h>
#include <filesystem>
//#include <char.h>

using namespace std;
using namespace cv;
namespace fs = std::filesystem;


struct dataSet* loadDataSetFromImages(const string& sourcePath, vector<string>& file_names,
        int numImages, int width, int height, bool logging, bool transpose){
    vector<string> imageDirs = vector<string>();

    /*string source = R"(C:\Users\nikit\CLionProjects\GPFL\images2gpfl\)";
    char *src2find = strdup(R"(C:\Users\nikit\CLionProjects\GPFL\images2gpfl\*)");*/

    int numRes = width * height;

    double* outputs;
    double* inputs;
    struct dataSet* dataSet2return;

    Mat image;

    /*WIN32_FIND_DATA FindFileData;
    HANDLE hf;*/

    for (const auto & entry: fs::directory_iterator(sourcePath)) {
//        cout << "______" << entry.path().string() << endl;
        imageDirs.emplace_back(entry.path().string());
	    file_names.emplace_back(entry.path().stem().string());
    }


/*    hf=FindFirstFile(src2find, &FindFileData);

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
    }*/

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
        if (numImages == 0) {
            cout << "Images not found!" << endl;
            return nullptr;
        }
    }

    outputs = new double[numImages * numRes];
    inputs = new double[2 * numImages * numRes];

    for (int i = 0; i < numImages; i++) {
        /*if (logging){
            cout << "Reading: " << imageDirs[i] << " ..." << endl;
        }*/
        image = imread(imageDirs[i], IMREAD_GRAYSCALE);

        if (transpose){
            image = image.t();
        }

        /*std::string greyArrWindow = String("Image") + to_string(i);
        cv::namedWindow(greyArrWindow, cv::WINDOW_NORMAL);
        cv::imshow(greyArrWindow, image);
        waitKey(0);*/

        uchar *arr = image.isContinuous() ? image.data : image.clone().data;

        for (int k = 0; k < numRes; k++) {
            inputs[2*i * numRes + 2 * k] = double(k / width);
            inputs[2*i * numRes + 2 * k + 1] = double(k % width);
            outputs[i * numRes + k] = double(arr[k]) / 255.0;
        }
    }

    dataSet2return = initialiseDataSetFromArrays(2, 1, numImages * numRes, inputs, outputs);
    cout << "Recording finished!" << endl;

    delete[] inputs;
    delete[] outputs;

    return dataSet2return;
}

struct dataSet* updateDataSet(dataSet* oldDataSet, chromosome* chromo, parameters* params, [[maybe_unused]] bool logging){
    int numImages = getNumImages(params);
    int numRes = getImageResolution(params);
    int width = getWidth(params);
    struct dataSet* dataSet2return;

    auto* outputs = new double[numImages * numRes];
    auto* inputs = new double[2 * numImages * numRes];

    for (int i = 0; i < numImages; i++) {
        double a = getA(chromo, i);
        double b = getB(chromo, i);

        for (int p = 0; p < numRes; p++) {
            inputs[2*i * numRes + 2 * p] = double(p / width);
            inputs[2*i * numRes + 2 * p + 1] = double(p % width);

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
