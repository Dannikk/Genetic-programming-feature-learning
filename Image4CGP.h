//
// Created by nikit on 25.3.22.
//

#include "cgp.h"
#include <string>

#ifndef GPFL_IMAGE4CGP_H
#define GPFL_IMAGE4CGP_H

using namespace std;

struct dataSet* loadDataSetFromImages(const string& sourcePath, int numImages, int width, int height, bool logging=false);


struct dataSet* updateDataSet(dataSet* oldDataSet, chromosome* chromo, parameters* params, [[maybe_unused]] bool logging = false);



#endif //GPFL_IMAGE4CGP_H
