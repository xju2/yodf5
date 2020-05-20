#ifndef __YODF5_CONST_KEYNAMES_H__
#define __YODF5_CONST_KEYNAMES_H__

#include <string>
#include <vector>

namespace YODF5 {
    const std::string BINID = "binids";
    const std::string VARIATIONS = "variations";

    const std::string SUMW = "sumw";
    const std::string SUMW2 = "sumw2";
    const std::string SUMWX = "sumwx";
    const std::string SUMWX2 = "sumwx2";
    const std::string SUMWY = "sumwy";
    const std::string SUMWY2 = "sumwy2";
    const std::string SUMWXY = "sumwxy";
    const std::string NUMENTRIES = "numEntries";
    const std::string XVAL = "xval";
    const std::string XERRM = "xerr-";
    const std::string XERRP = "xerr+";
    const std::string YVAL = "yval";
    const std::string YERRM = "yerr-";
    const std::string YERRP = "yerr+";
    const std::string XMIN = "xmin";
    const std::string XMAX = "xmax";
    const std::string YMIN = "ymin";
    const std::string YMAX = "ymax";
    const std::vector<std::string> ALL_COLUMNS {
        SUMW, SUMW2, SUMWX, SUMWX2, SUMWY, SUMWY2, SUMWXY,
        NUMENTRIES, XVAL, XERRM, XERRP,
        YVAL, YERRM, YERRP,
        XMIN, XMAX, YMIN, YMAX
    };

    const std::string NA = "N/A";
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