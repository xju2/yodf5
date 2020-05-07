#ifndef __YODF5_H5UTILS_H__
#define __YODF5_H5UTILS_H__
/**
 * This provides util functions that explore the H5 file contents 
 * for estimating the theory systematic uncertainties
 * */

#include <vector>
#include <string>

#include "xtensor/xarray.hpp"

namespace YODF5{
    namespace H5Utils {

        std::vector<float> get_pdf_sys(
                const xt::xarray<float>& all_obs_values,
                const std::string& pdf_name,
                const std::vector<int>& variation_idx);

        // sumw2 for Histo1D, Histo2D and Profile1D
        // yval for Scatter1D
        xt::xarray<double> get_obs_values(std::string& file_name);
        
        int get_pdf_id(const std::string& variation);
        int find_nominal_pdfid(const std::string& filename);
      
        // std::vector<std::string> > get_variations(std::string& filename);
    };
};
#endif