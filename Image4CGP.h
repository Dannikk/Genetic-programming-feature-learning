//
// Created by nikit on 25.3.22.
//

#include "cgp.h"

#ifndef GPFL_IMAGE4CGP_H
#define GPFL_IMAGE4CGP_H


struct dataSet* loadDataSetFromImages(char* sourceDir, int numImages, int width, int height, bool logging=false);


struct dataSet* updateDataSet(dataSet* oldDataSet, chromosome* chromo, parameters* params, [[maybe_unused]] bool logging = false);



#endif //GPFL_IMAGE4CGP_H
