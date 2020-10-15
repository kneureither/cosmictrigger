//
// Created by Konstantin Neureither on 06.08.20.
//

#ifndef COSMICTRIGGER_COSMICTEMPLATEEVAL_H
#define COSMICTRIGGER_COSMICTEMPLATEEVAL_H

#include <vector>

void cosmicTemplatesEfficiency(const int dataset, unsigned int centralTPcount, float spWZratio, float stopping_efficiency);
std::vector<std::vector<unsigned int>> produceCosmicSIDtracks(int cosmic_testing_dataset, int centralTPcount, float spWZratio, int mode);

#endif //COSMICTRIGGER_COSMICTEMPLATEEVAL_H


