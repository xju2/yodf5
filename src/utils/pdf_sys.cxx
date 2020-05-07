#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include <fstream>

#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xaxis_slice_iterator.hpp"
#include "xtensor/xadapt.hpp"

#define H5_USE_XTENSOR
#include <highfive/H5Easy.hpp>
#include <highfive/H5Group.hpp>

#include "yodf5/const_keynames.hpp"
#include "yodf5/lhapdf_lookup.hpp"
#include "yodf5/h5utils.hpp"

#include "LHAPDF/LHAPDF.h"


// #include <bits/stdc++.h> 
// #include <boost/algorithm/string.hpp> 

void print_shape(const xt::xarray<double>& a) {
    const auto& s = a.shape();
    std::copy(s.cbegin(), s.cend(), std::ostream_iterator<double>(std::cout, " "));
}

int main(int argc, char** argv) {
    std::string filename = "Rivet.h5";
    std::string output_filename = "Rivet_sys.h5";
    bool help = false;
    int opt;
    while ((opt = getopt(argc, argv, "hf:o:")) != -1) {
        switch(opt) {
            case 'f':
                filename = optarg;
                break;
            case 'o':
                output_filename = optarg;
                break;
            case 'h':
                help = true;
            default:
                fprintf(stderr, "Usage: %s [-h] [-f FILENAME] [-o OUTPUT_FILENAME]\n", argv[0]);
                if (help) {
                    printf("    -f FILENAME : input HDF5 file. Default is \"Rivet.h5\"");
                    printf("    -o OUTPUT_FILENAME : output HDF5 file. Default is \"Rivet_sys.h5\"");
                    printf("    -h HELP : print help info");
                }
                exit(EXIT_FAILURE);
        }
    }

    std::ifstream f(filename);
    if(!f.good()) {
        std::cerr << "input file \"" << filename << "\" does not exist\n";
        exit(EXIT_FAILURE);
    }

    // read the h5 file
    H5Easy::File file(filename, H5Easy::File::ReadOnly);
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
            nominal_pdf_id = YODF5::H5Utils::get_pdf_id(variation);
            try {
                nominal_pdf_name = YODF5::pset_name(nominal_pdf_id);
            } catch(std::out_of_range& e){
                std::cerr << "Cannot find pdf ID" << nominal_pdf_id << std::endl;
            }
        }
        if (variation.find("MUR1_MUF1") != std::string::npos) {
            std::string pdf_name = YODF5::pset_name(YODF5::H5Utils::get_pdf_id(variation));
            pdfs_variations[pdf_name] ++;
            pdfs_variation_idx[pdf_name].push_back(idx);
        }
    }
    std::cout << "nominal pdf: " << nominal_pdf_name << std::endl;
    int n_external_pdfs = 0;
    for (const auto& pair : pdfs_variations) {
        std::cout << pair.first << ": " << pair.second << " variations" << std::endl;
        if (pair.first.find(nominal_pdf_name) == std::string::npos) {
            n_external_pdfs ++;
        }
    }
    std::cout <<  n_external_pdfs << " external pdfs" << std::endl;

    // for(const auto& key_name: YODF5::ALL_DATA){
    //     auto values
    // }
    auto value = H5Easy::load<xt::xarray<double>>(file, YODF5::SUMW);
    long unsigned int n_bins = value.shape(0);
    // cout << value << endl;
    // auto shape = value.shape();
    // cout << shape[0] << " " << shape[1] << endl;
    // auto sumw = H5Easy::load<xt::x>(file, YODF5::SUMW);

    auto obs_values = YODF5::H5Utils::get_obs_values(filename);

    auto pdf_sys = xt::xarray<double>::from_shape({n_bins, 3});
    // pdf_sys.fill(YODF5::INVALID_NUMBER);

    auto& nominal_idx = pdfs_variation_idx[nominal_pdf_name];
    std::vector<float> internal_pdf_sys = YODF5::H5Utils::get_pdf_sys(obs_values, nominal_pdf_name, nominal_idx);

    // evaluate pdf systematic uncertainties
    int inorm = -1;
    int idx = 0;
    auto external_pdf_values = xt::xarray<double>::from_shape({n_bins, (long unsigned int) n_external_pdfs+1});
    external_pdf_values.fill(0.);
    for (const auto& pair: pdfs_variations){
        const std::string& pdf_name = pair.first;
        size_t n_varies = pair.second;
        if(pdf_name == nominal_pdf_name) {
            inorm = idx;
            // continue;
        }
        if(n_varies > 1) {
            std::vector<float> this_pdf_sys = YODF5::H5Utils::get_pdf_sys(obs_values, pdf_name, pdfs_variation_idx[pdf_name]);
            std::vector<size_t> shape = {this_pdf_sys.size()};
            xt::col(external_pdf_values, idx) = xt::adapt(this_pdf_sys, shape);
        } else {
            xt::col(external_pdf_values, idx) = xt::col(obs_values, pdfs_variation_idx[pdf_name][0]);
        }
        idx ++;
    }
    std::vector<size_t> sys_shape = {n_bins, 1};
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