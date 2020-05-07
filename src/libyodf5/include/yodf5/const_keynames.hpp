#ifndef __YODF5_CONST_KEYNAMES_H__
#define __YODF5_CONST_KEYNAMES_H__

#include <string>
#include <vector>

namespace YODF5 {
    const std::string BINID = "binids";
    const std::string VARIATIONS = "variations";
    const std::string SUMW = "sumw";
    const std::string YVAL = "yval";

    const std::string H1D = "Histo1D";
    const std::string H2D = "Histo2D";
    const std::string P1D = "Profile1D";
    const std::string S1D = "Scatter1D";
    const std::string S2D = "Scatter2D";
    const std::string CNT = "Counter";
    const std::vector<std::string> ALL_DATA {
        H1D, H2D, P1D, S1D, S2D    
    };
    const std::vector<std::string> ALL_HIST {
        H1D, H2D, P1D
    };

    const double INVALID_NUMBER = -999999;
};
#endif