#include <iostream>
#include <string>
#include <vector>

#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xaxis_slice_iterator.hpp"
#include "xtensor/xadapt.hpp"

#define H5_USE_XTENSOR
#include <highfive/H5Easy.hpp>
#include <highfive/H5Group.hpp>

#include "const_keynames.hpp"
#include "lhapdf_lookup.hpp"

#include "LHAPDF/LHAPDF.h"


// #include <bits/stdc++.h> 
// #include <boost/algorithm/string.hpp> 

using namespace std;

void print_shape(const xt::xarray<double>& a) {
    const auto& s = a.shape();
    std::copy(s.cbegin(), s.cend(), std::ostream_iterator<double>(std::cout, " "));
}

vector<float> get_pdf_sys(
    const xt::xarray<float>& all_obs_values,
    const string& pdf_name,
    const vector<int>& variation_idx
){
    const LHAPDF::PDFSet this_pdf(pdf_name);
    vector<float> pdf_sys;
    int cl = 68;

    auto get_pdf_sys_fn = [&](const xt::xarray<double>& array)
    {
        vector<double> varies;
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

xt::xarray<double> get_obs_values(string& file_name){
    H5Easy::File file(file_name, H5Easy::File::ReadOnly);
    auto sumw = H5Easy::load<xt::xarray<double>>(file, YODF5::SUMW);
    cout << sumw << endl;
    auto yval = H5Easy::load<xt::xarray<double>>(file, YODF5::YVAL);
    long unsigned int n_variations = sumw.shape(1);
    long unsigned int n_bins = sumw.shape(0);
    auto obs_values = xt::xarray<double>::from_shape({n_bins, n_variations});
    obs_values.fill(YODF5::INVALID_NUMBER);
    
    // the measured values of histograms (1D, 2D and profile 1D) are from sumw
    // and the measured values of scatters from yval;
    for (const auto& hist_keyname: YODF5::ALL_HIST) {
        HighFive::Group g = file.getGroup(hist_keyname);
        const vector<string> hist_names = g.listObjectNames();
        for (const auto& hist_name: hist_names){
            auto bin_idx = H5Easy::load<vector<int>>(file, hist_keyname+"/"+hist_name);
            for(int bin: bin_idx) {
                for(int ivar = 0; ivar < n_variations; ivar++){
                    obs_values(bin, ivar) = sumw(bin, ivar);
                }
            }
        }
    }

    HighFive::Group g = file.getGroup(YODF5::S1D);
    const vector<string> scatter_names = g.listObjectNames();
    for(const auto& scatter_name: scatter_names) {
        auto bin_idx = H5Easy::load<vector<int>>(file, YODF5::S1D+"/"+scatter_name);
        for(int bin: bin_idx) {
            for(int ivar = 0; ivar < n_variations; ivar++){
                obs_values(bin, ivar) = yval(bin, ivar);
            }
        }
    }
    // cout << obs_values << endl;
    return obs_values;
}

std::vector<std::string> getListOfVariations(std::string const & fname) {
    H5Easy::File file(fname, H5Easy::File::ReadOnly);
    auto vars = H5Easy::load<std::vector<std::string> >(file, "variations");
    return vars;
}

void split(vector<string>& results, string input, string delimiter){
    std::string::size_type pos;
    while ( (pos = input.find(delimiter)) != std::string::npos) {
        results.push_back(input.substr(0, pos));
        input.erase(0, pos+delimiter.length());
    }
    results.push_back(input);
}
int get_pdf_id(std::string& variation) {
    vector<string> items;
    items.clear();
    split(items, variation, "_");
    string& item = items[items.size()-1];
    return stoi(item.substr(3, item.size()-3));
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << argv[0] << " file-name" << std::endl;
        return 1;
    }
    std::string file_name(argv[1]);

    // read the h5 file
    H5Easy::File file(file_name, H5Easy::File::ReadOnly);
    std::vector<std::string> variations = H5Easy::load<std::vector<std::string> >(file, YODF5::VARIATIONS);
    std::cout << "total " << variations.size() << " variations." << std::endl;

    auto binids = H5Easy::load<std::vector<std::string>>(file, YODF5::BINID);
    std::cout << "total " << binids.size() << " bins." << std::endl;

    std::string nominal_pdf_name;
    int nominal_pdf_id = -1;
    std::map<std::string, std::size_t> pdfs_variations;
    std::map<std::string, std::vector<int> > pdfs_variation_idx;
    for(std::vector<std::string>::size_type idx = 0; idx != variations.size(); idx++) {
        std::string& variation = variations[idx];
        if (variation.find("MUR1_MUF0.5") != std::string::npos) {
            nominal_pdf_id = get_pdf_id(variation);
            try {
                nominal_pdf_name = YODF5::pset_name(nominal_pdf_id);
            } catch(out_of_range& e){
                cerr << "Cannot find pdf ID" << nominal_pdf_id << endl;
            }
        }
        if (variation.find("MUR1_MUF1") != string::npos) {
            string pdf_name = YODF5::pset_name(get_pdf_id(variation));
            pdfs_variations[pdf_name] ++;
            pdfs_variation_idx[pdf_name].push_back(idx);
        }
    }
    cout << "nominal pdf: " << nominal_pdf_name << endl;
    int n_external_pdfs = 0;
    for (const auto& pair : pdfs_variations) {
        cout << pair.first << ": " << pair.second << " variations" << endl;
        if (pair.first.find(nominal_pdf_name) == string::npos) {
            n_external_pdfs ++;
        }
    }
    cout <<  n_external_pdfs << " external pdfs" << endl;

    // for(const auto& key_name: YODF5::ALL_DATA){
    //     auto values
    // }
    auto value = H5Easy::load<xt::xarray<double>>(file, YODF5::SUMW);
    long unsigned int n_bins = value.shape(0);
    // cout << value << endl;
    // auto shape = value.shape();
    // cout << shape[0] << " " << shape[1] << endl;
    // auto sumw = H5Easy::load<xt::x>(file, YODF5::SUMW);

    auto obs_values = get_obs_values(file_name);

    auto pdf_sys = xt::xarray<double>::from_shape({n_bins, 3});
    // pdf_sys.fill(YODF5::INVALID_NUMBER);

    auto& nominal_idx = pdfs_variation_idx[nominal_pdf_name];
    vector<float> internal_pdf_sys = get_pdf_sys(obs_values, nominal_pdf_name, nominal_idx);

    // evaluate pdf systematic uncertainties
    int inorm = -1;
    int idx = 0;
    auto external_pdf_values = xt::xarray<double>::from_shape({n_bins, (long unsigned int) n_external_pdfs+1});
    external_pdf_values.fill(0.);
    for (const auto& pair: pdfs_variations){
        const string& pdf_name = pair.first;
        size_t n_varies = pair.second;
        if(pdf_name == nominal_pdf_name) {
            inorm = idx;
            // continue;
        }
        if(n_varies > 1) {
            vector<float> this_pdf_sys = get_pdf_sys(obs_values, pdf_name, pdfs_variation_idx[pdf_name]);
            vector<size_t> shape = {this_pdf_sys.size()};
            xt::col(external_pdf_values, idx) = xt::adapt(this_pdf_sys, shape);
        } else {
            xt::col(external_pdf_values, idx) = xt::col(obs_values, pdfs_variation_idx[pdf_name][0]);
        }
        idx ++;
    }
    vector<size_t> sys_shape = {n_bins, 1};
    // auto nominal_pdf_sys = xt::col(external_pdf_values, inorm);
    // nominal_pdf_sys.reshape({nbins, 1});
    // auto external_diff = external_pdf_values - xt::col(external_pdf_values, inorm).reshape({n_bins, 1});
    auto external_diff = external_pdf_values - xt::adapt(internal_pdf_sys, sys_shape);
    auto external_pdf_sys = xt::amax(external_diff, 1);

    sys_shape = {n_bins};
    xt::col(pdf_sys, 0) = xt::adapt(internal_pdf_sys, sys_shape);
    xt::col(pdf_sys, 1) = external_pdf_sys;
    xt::col(pdf_sys, 2) = xt::amax(xt::view(pdf_sys, xt::all(), xt::range(0, 2)), 1);

    H5Easy::File outfile("Rivet_pdfsys_cpp.h5", H5Easy::File::Create);
    H5Easy::dump(outfile, "pdfsys", pdf_sys);

    return 0;
}