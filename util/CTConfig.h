//
// Created by Konstantin Neureither on 17.09.20.
//

#ifndef COSMICTRIGGER_CTCONFIG_H
#define COSMICTRIGGER_CTCONFIG_H


class CTConfig {
public:
    std::vector<int> spcounts;
    std::vector<float> spratios;
    std::vector<float> stopping_efficiencies;

    int dataset;
    int bg_run;
    int pe_mode;

    //TDB build
    int combination_id;
    int max_cosmic_entries;
    bool write_tdb_file;
    bool make_plot;

    //CosPatPlots
    std::vector<int> cycle_plotting_order;


    //TODO Also handle plotting files
private:

};


#endif //COSMICTRIGGER_CTCONFIG_H
