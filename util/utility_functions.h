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

unsigned int get_layer(int sid) {
    if(0 <= sid && sid < 1024) {
        return 0;
    } else if (1024 <= sid && sid < 2048) {
        return 1;
    } else if ((2000 <= sid && sid < 3000) || (10000 <= sid && sid < 11500) || (14000 <= sid && sid < 15200)) {
        return 2;
    } else if ((3000 <= sid && sid < 4000) || (11500 <= sid && sid < 12500) || (15200 <= sid && sid < 16500)) {
        return 3;
    } else {
        return 100;
    }
}