#ifndef __YODF5_CONST_KEYNAMES_H__
#define __YODF5_CONST_KEYNAMES_H__

#include <string>
#include <vector>
using namespace std;

namespace YODF5 {
    const string BINID = "binids";
    const string VARIATIONS = "variations";
    const string SUMW = "sumw";
    const string YVAL = "yval";

    const string H1D = "Histo1D";
    const string H2D = "Histo2D";
    const string P1D = "Profile1D";
    const string S1D = "Scatter1D";
    const string S2D = "Scatter2D";
    const string CNT = "Counter";
    const vector<string> ALL_DATA {
        H1D, H2D, P1D, S1D, S2D    
    };
    const vector<string> ALL_HIST {
        H1D, H2D, P1D
    };

    const double INVALID_NUMBER = -999999;
};
#endif