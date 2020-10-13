//
// Created by Konstantin Neureither on 13.10.20.
//
#include <vector>

std::vector<std::vector<unsigned int>> produceCosmicSIDtracks(int cosmic_testing_dataset, int centralTPcount, float spWZratio, int mode);

int main(int argc, char *argv[]) {
    std::vector <std::vector<unsigned int>> result1;
    std::vector <std::vector<unsigned int>> result2;
    int cosmic_testing_dataset = 30;
    int centralTPcount = 400;
    float spWZratio = 1;
    int mode = 0;


    result1 = produceCosmicSIDtracks(cosmic_testing_dataset, centralTPcount, spWZratio, mode);

    result2 = produceCosmicSIDtracks(cosmic_testing_dataset, centralTPcount, spWZratio, mode);

    assert(result1 == result2);
}
