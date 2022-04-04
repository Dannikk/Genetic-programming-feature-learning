//
// Created by nikit on 2.4.22.
//

#include "CustomNodeFunctions.h"

double min(const int numInputs, const double *inputs, const double *connectionWeights) {

    int i;
    double currentMin= inputs[0];

    for (i = 1; i < numInputs; i++) {
        if (inputs[i] < currentMin)
            currentMin = inputs[i];
    }

    return currentMin;
}

double max(const int numInputs, const double *inputs, const double *connectionWeights) {

    int i;
    double currentMax= inputs[0];

    for (i = 1; i < numInputs; i++) {
        if (inputs[i] > currentMax)
            currentMax = inputs[i];
    }

    return currentMax;
}
