//
// Created by Konstantin Neureither on 29.05.20.
//
#include <iostream>

int main(int argc, char *argv[]) {
    int run;
    int filter;
    if(argc < 2) {
        std::cout << "ERROR: Error in argument! Usage: Mu3eRecAcc <run number> <filter arg> " << std::endl;
        std::cout << "\t Find filter options in filter_usage.json" << std::endl;
        exit(0);
    } else {
        run = atoi(argv[1]);
        filter = atoi(argv[2]);
    }

    std::cout << "Running reconstruction_accuracy() for run " << run << "..." << std::endl;

    void reconstruction_accuracy(int, int);
    reconstruction_accuracy(run, filter);

    return 0;
}
