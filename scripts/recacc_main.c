//
// Created by Konstantin Neureither on 29.05.20.
//
#include <iostream>

int main(int argc, char *argv[]) {
    int run;
    if(argc < 2) {
        std::cout << "ERR : No argument! Usage Mu3eRecAcc <run number> " << std::endl;
        exit(0);
    } else {
        run = atoi(argv[1]);
    }

    std::cout << "Running reconstruction_accuracy() for run " << run << "..." << std::endl;

    void reconstruction_accuracy(int);
    reconstruction_accuracy(run);
}
