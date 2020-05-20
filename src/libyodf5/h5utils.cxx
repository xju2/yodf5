#include "yodf5/h5utils.hpp"

#include "yodf5/const_keynames.hpp"
#include "yodf5/lhapdf_lookup.hpp"

#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xaxis_slice_iterator.hpp"
#include "xtensor/xadapt.hpp"

#define H5_USE_XTENSOR
#include <highfive/H5Easy.hpp>
#include <highfive/H5Group.hpp>

#include "LHAPDF/LHAPDF.h"

namespace YODF5{
    namespace H5Utils {

    void split(std::vector<std::string>& results, std::string input, const std::string&  delimiter){
        std::string::size_type pos;
        while ( (pos = input.find(delimiter)) != std::string::npos) {
            results.push_back(input.substr(0, pos));
            input.erase(0, pos+delimiter.length());
        }
        results.push_back(input);
    }

    int get_pdf_id(const std::string& variation) {
        std::vector<std::string> items;
        items.clear();
        split(items, variation, "_");
        std::string& item = items[items.size()-1];
        return std::stoi(item.substr(3, item.size()-3));
    }

    std::vector<float> get_pdf_sys(
            const xt::xarray<float>& all_obs_values,
            const std::string& pdf_name,
            const std::vector<int>& variation_idx
        ){
            const LHAPDF::PDFSet this_pdf(pdf_name);
            std::vector<float> pdf_sys;
            int cl = 68;

            auto get_pdf_sys_fn = [&](const xt::xarray<double>& array)
            {
                std::vector<double> varies;
                for(size_t idx: variation_idx) {
                    varies.push_back(array(idx));
                }
                const LHAPDF::PDFUncertainty err = this_pdf.uncertainty(varies, cl);
                pdf_sys.push_back(err.errsymm);
            };
            auto iter = axis_slice_begin(all_obs_values, 1);
            auto end = axis_slice_end(all_obs_values, 1);
            std::for_each(iter, end, get_pdf_sys_fn);
            // cout << "internal pdf size: " << pdf_sys.size() << endl;
            // assert(pdf_sys.size() == all_obs_values.shape(0));

            return pdf_sys;
        }

        xt::xarray<double> get_obs_values(std::string& file_name){
            H5Easy::File file(file_name, H5Easy::File::ReadOnly);
            auto sumw = H5Easy::load<xt::xarray<double>>(file, YODF5::SUMW);
            std::cout << sumw << std::endl;
            auto yval = H5Easy::load<xt::xarray<double>>(file, YODF5::YVAL);
            long unsigned int n_variations = sumw.shape(1);
            long unsigned int n_bins = sumw.shape(0);
            auto obs_values = xt::xarray<double>::from_shape({n_bins, n_variations});
            obs_values.fill(YODF5::INVALID_NUMBER);
            
            // the measured values of histograms (1D, 2D and profile 1D) are from sumw
            // and the measured values of scatters from yval;
            for (const auto& hist_keyname: YODF5::ALL_HIST) {
                HighFive::Group g = file.getGroup(hist_keyname);
                const std::vector<std::string> hist_names = g.listObjectNames();
                for (const auto& hist_name: hist_names){
                    auto bin_idx = H5Easy::load<std::vector<int>>(file, hist_keyname+"/"+hist_name);
                    for(int bin: bin_idx) {
                        for(int ivar = 0; ivar < n_variations; ivar++){
                            obs_values(bin, ivar) = sumw(bin, ivar);
                        }
                    }
                }
            }

            HighFive::Group g = file.getGroup(YODF5::S1D);
            const std::vector<std::string> scatter_names = g.listObjectNames();
            for(const auto& scatter_name: scatter_names) {
                auto bin_idx = H5Easy::load<std::vector<int>>(file, YODF5::S1D+"/"+scatter_name);
                for(int bin: bin_idx) {
                    for(int ivar = 0; ivar < n_variations; ivar++){
                        obs_values(bin, ivar) = yval(bin, ivar);
                    }
                }
            }
            // cout << obs_values << endl;
            return obs_values;
        }

        int find_nominal_pdfid(const std::string& filename){
            H5Easy::File file(filename, H5Easy::File::ReadOnly);
            auto variations = H5Easy::load<std::vector<std::string>>(file, YODF5::VARIATIONS);
            int nominal_pdf_id = -1;
            for(const auto& variation: variations) {
                if (variation.find("MUR1_MUF0.5") != std::string::npos) {
                    nominal_pdf_id = get_pdf_id(variation);
                    break;
                }
            }
            if(nominal_pdf_id < 0) {
                throw(std::range_error("Cannot find nominal PDF ID"));
            }
            return nominal_pdf_id;
        }
    };
};