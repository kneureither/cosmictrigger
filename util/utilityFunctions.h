#ifndef COSMICTRIGGER_UTILITY_FUNCTIONS_H
#define COSMICTRIGGER_UTILITY_FUNCTIONS_H

#endif //COSMICTRIGGER_UTILITY_FUNCTIONS_H

#include <string.h>
#include <iostream>
#include <experimental/filesystem>
#include "custom_types.h"

using std::endl;
using std::cout;
namespace fs = std::experimental::filesystem;


// pad an int number
std::string get_padded_string(int number, int n, char c) {
    std::stringstream ss;
    ss << number;
    std::string str = ss.str();
    str.insert(str.begin(), n-str.length(), c);
    return str;
}

template <class T>
std::string get_string(T number) {
    std::stringstream ss;
    ss << number;
    return ss.str();
}

PXID process_pixel_id(unsigned int pixel_id) {
    PXID pixid;
    pixid.sensor        = (pixel_id)>>16;
    pixid.column        = ((pixel_id)>>8)&0xFF;
    pixid.row           = ((pixel_id))&0xFF;
    pixid.columnaddress = pixid.column + (pixid.sensor << 8);

    return pixid;
}

template <class T>
float vector_mean(std::vector<T> vec) {
    float sum = 0.0;
    for(int i = 0; i < vec.size(); i++) {
        sum += vec[i] / (float) vec.size();
    }
    return sum;
}

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


void check_create_directory(std::string path) {
    if (fs::exists(path.c_str())) {
        cout << "STATUS : path exists: " << path << endl;
    } else {
        fs::create_directory(path);
        cout << "STATUS : no such path - created: " << path << endl;
    }
}

int get_charge_from_type(int type) {
    if (type == 4) {
        //mu-
        return 1;
    } else if (type == 3) {
        //mu+
        return -1;
    } else {
        //not a muon
        return 1;
    }
}