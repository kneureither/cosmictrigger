//
// Created by Konstantin Neureither on 09.01.21.
//

#ifndef COSMICTRIGGER_GETCOSMICSIDTRACKS_H
#define COSMICTRIGGER_GETCOSMICSIDTRACKS_H

#include <vector>

std::vector<std::vector<unsigned int>> getCosmicSIDtracks(int cosmic_testing_dataset, int centralTPcount, float spWZratio, int mode);

#endif //COSMICTRIGGER_GETCOSMICSIDTRACKS_H