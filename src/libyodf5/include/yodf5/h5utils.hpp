#ifndef __YODF5_H5UTILS_H__
#define __YODF5_H5UTILS_H__
/**
 * This provides util functions that explore the H5 file contents 
 * for estimating the theory systematic uncertainties
 * */

#include <vector>
#include <string>
#include <unordered_map>

#include "YODA/YODA.h" // include all analysis objects
#include "const_keynames.hpp"

#include "xtensor/xarray.hpp"

#define H5_USE_XTENSOR
#include <highfive/H5Easy.hpp>

namespace YODF5{
    namespace H5Utils {
        void split(std::vector<std::string>& results, 
            std::string input, const std::string& delimiter);

        std::vector<float> get_pdf_sys(
                const xt::xarray<float>& all_obs_values,
                const std::string& pdf_name,
                const std::vector<int>& variation_idx);

        // sumw2 for Histo1D, Histo2D and Profile1D
        // yval for Scatter1D
        xt::xarray<double> get_obs_values(std::string& file_name);
        
        int get_pdf_id(const std::string& variation);
        int find_nominal_pdfid(const std::string& filename);
      
        void writeYoda2H5(std::vector<YODA::AnalysisObject*> const &, std::string const & outh5_name, int compressio=1);

        // functions for converting yoda to h5.
        template <typename T> 
        void print_shape(xt::xarray<T>& a){
            const auto& s = a.shape();
            std::copy(s.cbegin(), s.cend(), std::ostream_iterator<double>(std::cout, " "));
        }

        template <typename T>
        xt::xarray<double> dbn0ToArray(T* dbn){
            xt::xarray<double> array({3}, 0);
            array(0) = dbn->sumW();
            array(1) = dbn->sumW2();
            array(2) = dbn->numEntries();
            return std::move(array);
        }

        std::string replaceall(std::string const & in_str, std::string const & source, std::string const & target)
        {
            std::string::size_type pos;
            std::string results(in_str);
            while ( (pos = results.find(source)) != std::string::npos) {
                results.replace(pos, source.size(), target);
            }
            return results;
        }

        int num_bins(YODA::AnalysisObject* ao){
            std::string aotype = ao->type();
            int nbins = 1;
            if(aotype == YODF5::H1D) nbins = dynamic_cast<YODA::Histo1D*>(ao)->numBins();
            else if(aotype == YODF5::H2D) nbins = dynamic_cast<YODA::Histo2D*>(ao)->numBins();
            else if(aotype == YODF5::P1D) nbins = dynamic_cast<YODA::Profile1D*>(ao)->numBins();
            else if(aotype == YODF5::S1D) nbins = dynamic_cast<YODA::Scatter1D*>(ao)->numPoints();
            else if(aotype == YODF5::S2D) nbins = dynamic_cast<YODA::Scatter2D*>(ao)->numPoints();
            else {
                std::cerr << "Not known type: " << aotype << ". Number of bins set to 1." << std::endl;
            }
            return nbins;
        }

        // convert YODA::AnalysisObject to xt::xarray
        // for each bin in Histo1D
        xt::xarray<double> dbn1ToArray(YODA::HistoBin1D const & dbn);
        // for over/underflow bins
        xt::xarray<double> dbn1ToArray(YODA::Dbn1D const & dbn);
        // for each bin in profile1D
        xt::xarray<double> dbn1ToArray(YODA::ProfileBin1D const & dbn);
        // for each bin in Hist2D
        xt::xarray<double> Hdbn2ToArray(YODA::HistoBin2D const & dbn);
        // for over/underflow bins, 2D
        xt::xarray<double> Hdbn2ToArray(YODA::Dbn2D const & dbn);
        // for profile 1D total-bin
        xt::xarray<double> dbn2ToArray(YODA::Dbn2D const & dbn);

        xt::xarray<double> point1ToArray(YODA::Point1D const & pnt);
        xt::xarray<double> point2ToArray(YODA::Point2D const & pnt);
        
        std::vector<std::string> get_ao_names(std::vector<YODA::AnalysisObject*> const & aos, std::string const & aotype);
        std::vector<std::string> get_variations(std::vector<YODA::AnalysisObject*> const & aos, std::string const & hname);

        std::vector<std::string> mk_binids(std::vector<YODA::AnalysisObject*> const & aos);
        void  fill_dataset(
            H5Easy::File& file,
            std::unordered_map<std::string, YODA::AnalysisObject*>& aomap,
            std::string const & hname,
            std::vector<std::string> const & binids,
            std::vector<std::string> const & variations
        );


    };
};
#endif