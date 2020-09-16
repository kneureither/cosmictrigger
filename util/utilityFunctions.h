#ifndef COSMICTRIGGER_UTILITY_FUNCTIONS_H
#define COSMICTRIGGER_UTILITY_FUNCTIONS_H

#include <sstream>
#include <string>
#include <iostream>
#include <experimental/filesystem>
#include "custom_types.h"

using std::endl;
using std::cout;
namespace fs = std::experimental::filesystem;


// Configuration variables


// pad an int number
static std::string get_padded_string(int number, int n, char c) {
    std::stringstream ss;
    ss << number;
    std::string str = ss.str();
    str.insert(str.begin(), n-str.length(), c);
    return str;
}

template <class T>
static std::string get_string(T number) {
    std::stringstream ss;
    ss << number;
    return ss.str();
}

static PXID process_pixel_id(unsigned int pixel_id) {
    PXID pixid;
    pixid.sensor        = (pixel_id)>>16;
    pixid.column        = ((pixel_id)>>8)&0xFF;
    pixid.row           = ((pixel_id))&0xFF;
    pixid.columnaddress = pixid.column + (pixid.sensor << 8);

    return pixid;
}

template <class T>
static float vector_mean(std::vector<T> vec) {
    float sum = 0.0;
    for(int i = 0; i < vec.size(); i++) {
        sum += vec[i] / (float) vec.size();
    }
    return sum;
}

template <typename T>
static int sgn(const T val) {
    return (T(0) < val) - (val < T(0));
}


static void check_create_directory(std::string path) {
    if (fs::exists(path.c_str())) {
        cout << "STATUS : path exists: " << path << endl;
    } else {
        fs::create_directory(path);
        cout << "STATUS : no such path - created: " << path << endl;
    }
}

static void print_status_bar(int entry, int max_entries, std::string label, std::string info) {
    const int STATUS_BAR_LEN = 50;
    float prog_perc = entry /  (float) max_entries;
    std::string prog_bar_fill((int) (STATUS_BAR_LEN * prog_perc), '=');
    std::string prog_bar_empty((int) (STATUS_BAR_LEN * (1-prog_perc)), ' ');
    std::cout << "\r(STATUS) : " << label <<  " [" << prog_bar_fill << ">" << prog_bar_empty << "] ";

    if(prog_bar_empty != ""){
        std::cout << entry /  (float) max_entries * 100 << "% | " << info << std::flush;
    } else {
        std::cout << 100 << "% | " << info << std::flush;
    }
}

static std::string getfileidtag(int mydataset, int mode, int wbins, int zbins) {
    return "dataset" + get_string(mydataset) + "_mode" + get_string(mode) + "wBins" + get_string(wbins) + "zBins" + get_string(zbins);
}

static std::string getfileidtag(int mydataset, int mode, int wbins, int zbins, float stopping_efficiency) {
    return "dataset" + get_string(mydataset) + "_mode" + get_string(mode) + "wBins" + get_string(wbins) + "zBins" + get_string(zbins) + "_maxeff" + get_string(stopping_efficiency);
}

static std::string getfileidtag(int mode, int wbins, int zbins) {
    return "mode" + get_string(mode) + "wBins" + get_string(wbins) + "zBins" + get_string(zbins);
}


#endif //COSMICTRIGGER_UTILITY_FUNCTIONS_H