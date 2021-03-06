//
// Created by Konstantin Neureither on 13.10.20.
//
#include <vector>
#include "inc/getCosmicSIDtracks.h"

int main(int argc, char *argv[]) {
    /**
     * This routine is just a basic a test for the method
     * getCosmicSIDtracks(...)
     *
     * TODO Move to test/ folder
     */

    std::vector <std::vector<unsigned int>> result1;
    std::vector <std::vector<unsigned int>> result2;
    int cosmic_testing_dataset = 30;
    int centralTPcount = 400;
    float spWZratio = 1;
    int mode = 0;


    //for testing: First call returns the data directly from calculation and writes file
    result1 = getCosmicSIDtracks(cosmic_testing_dataset, centralTPcount, spWZratio, mode);

    //second call reads the data from the file
    result2 = getCosmicSIDtracks(cosmic_testing_dataset, centralTPcount, spWZratio, mode);

    //if called the first time (no file exists), this checks if the data writing and reading worked correctly.
    assert(result1 == result2);
}
