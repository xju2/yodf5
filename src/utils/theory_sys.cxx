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
                    printf("    -f FILENAME : input HDF5 file. Default is \"Rivet.h5\"\n");
                    printf("    -o OUTPUT_FILENAME : output HDF5 file. Default is \"Rivet_sys.h5\"\n");
                    printf("    -h HELP : print help info\n");
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

    const std::string scale_sys_name("QCDScales");

    int nominal_pdf_id = YODF5::H5Utils::find_nominal_pdfid(filename);
    std::string nominal_pdf_name = YODF5::pset_name(nominal_pdf_id);
    std::map<std::string, std::size_t> pdfs_variations;
    std::map<std::string, std::vector<int> > sys_variation_idx;
    for(std::vector<std::string>::size_type idx = 0; idx != variations.size(); idx++) {
        std::string& variation = variations[idx];
        if (variation.find("MUR1_MUF1") != std::string::npos) {
            std::string pdf_name = YODF5::pset_name(YODF5::H5Utils::get_pdf_id(variation));
            pdfs_variations[pdf_name] ++;
            sys_variation_idx[pdf_name].push_back(idx);
        }
        if (variation.find(nominal_pdf_name)){
            sys_variation_idx[scale_sys_name].push_back(idx);
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

    auto value = H5Easy::load<xt::xarray<double>>(file, YODF5::SUMW);
    long unsigned int n_bins = value.shape(0);

    auto obs_values = YODF5::H5Utils::get_obs_values(filename);

    auto pdf_sys = xt::xarray<double>::from_shape({n_bins, 3});

    auto& nominal_idx = sys_variation_idx[nominal_pdf_name];
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
        }
        if(n_varies > 1) {
            std::vector<float> this_pdf_sys = YODF5::H5Utils::get_pdf_sys(obs_values, pdf_name, sys_variation_idx[pdf_name]);
            std::vector<size_t> shape = {this_pdf_sys.size()};
            xt::col(external_pdf_values, idx) = xt::adapt(this_pdf_sys, shape);
        } else {
            xt::col(external_pdf_values, idx) = xt::col(obs_values, sys_variation_idx[pdf_name][0]);
        }
        idx ++;
    }
    std::vector<size_t> sys_shape_ext = {n_bins, 1};
    auto external_diff = external_pdf_values - xt::adapt(internal_pdf_sys, sys_shape_ext);
    auto external_pdf_sys = xt::amax(external_diff, 1);

    std::vector<size_t> sys_shape = {n_bins};
    xt::col(pdf_sys, 0) = xt::adapt(internal_pdf_sys, sys_shape);
    xt::col(pdf_sys, 1) = external_pdf_sys;
    xt::col(pdf_sys, 2) = xt::amax(xt::view(pdf_sys, xt::all(), xt::range(0, 2)), 1);

    // QCD scale SYS
    size_t n_scale_variations = sys_variation_idx[scale_sys_name].size();
    auto scale_values = xt::xarray<double>::from_shape({n_bins, n_scale_variations});
    long unsigned int scale_idx = 0;
    for(const auto idx: sys_variation_idx[scale_sys_name]){
        xt::col(scale_values, scale_idx) = xt::col(obs_values, idx);
        scale_idx += 1;
    }
    auto scale_diff = scale_values - xt::adapt(internal_pdf_sys, sys_shape_ext);
    // auto scale_diff = xt::col(obs_values, sys_variation_idx[scale_sys_name]) - xt::adapt(internal_pdf_sys, sys_shape_ext);
    auto scale_sys = xt::xarray<double>::from_shape({n_bins, 1});
    xt::col(scale_sys, 0) = xt::amax(scale_diff, 1);

    H5Easy::File outfile("Rivet_pdfsys_cpp.h5", H5Easy::File::Create);
    H5Easy::dump(outfile, "pdf_sys", pdf_sys);
    H5Easy::dump(outfile, "scale_sys", scale_sys);

    return 0;
}