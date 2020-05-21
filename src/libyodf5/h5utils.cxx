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
            const std::vector<int>& variation_idx)
    {
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

    xt::xarray<double> get_obs_values(std::string& file_name)
    {
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

    // for each bin in Histo1D
    xt::xarray<double> dbn1ToArray(YODA::HistoBin1D const & dbn)
    {
        xt::xarray<double> array({7}, 0);
        array(0) = dbn.sumW();
        array(1) = dbn.sumW2();
        array(2) = dbn.sumWX();
        array(3) = dbn.sumWX2();
        array(4) = dbn.numEntries();
        array(5) = dbn.xMin();
        array(6) = dbn.xMax();
        return std::move(array);
    }

    // for over/underflow bins
    xt::xarray<double> dbn1ToArray(YODA::Dbn1D const & dbn) {
        xt::xarray<double> array({7}, 0);
        array(0) = dbn.sumW();
        array(1) = dbn.sumW2();
        array(2) = dbn.sumWX();
        array(3) = dbn.sumWX2();
        array(4) = dbn.numEntries();
        // array(5) = 0.;
        // array(6) = 0.;
        return std::move(array);
    }

    // for each bin in profile1D
    xt::xarray<double> dbn1ToArray(YODA::ProfileBin1D const & dbn)
    {
        xt::xarray<double> array({9}, 0);
        array(0) = dbn.sumW();
        array(1) = dbn.sumW2();
        array(2) = dbn.sumWX();
        array(3) = dbn.sumWX2();
        array(4) = dbn.sumWY();
        array(5) = dbn.sumWY2();
        array(6) = dbn.numEntries();
        array(7) = dbn.xMin();
        array(8) = dbn.xMax();
        return std::move(array);
    }

    // for each bin in Hist2D
    xt::xarray<double> Hdbn2ToArray(YODA::HistoBin2D const & dbn)
    {
        xt::xarray<double> array({12}, 0);
        array(0) = dbn.sumW();
        array(1) = dbn.sumW2();
        array(2) = dbn.sumWX();
        array(3) = dbn.sumWX2();
        array(4) = dbn.sumWY();
        array(5) = dbn.sumWY2();
        array(6) = dbn.sumWXY();
        array(7) = dbn.numEntries();
        array(8) = dbn.xMin();
        array(9) = dbn.xMax();
        array(10) = dbn.yMin();
        array(11) = dbn.yMax();
        return std::move(array);
    }

    // for over/underflow bins, 2D
    xt::xarray<double> Hdbn2ToArray(YODA::Dbn2D const & dbn)
    {
        xt::xarray<double> array({12}, 0);
        array(0) = dbn.sumW();
        array(1) = dbn.sumW2();
        array(2) = dbn.sumWX();
        array(3) = dbn.sumWX2();
        array(4) = dbn.sumWY();
        array(5) = dbn.sumWY2();
        array(6) = dbn.sumWXY();
        array(7) = dbn.numEntries();
        // array(8) = 0.;
        // array(9) = 0.;
        // array(10) = 0.;
        // array(11) = 0.;
        return std::move(array);
    }

    // for profile 1D 
    xt::xarray<double> dbn2ToArray(YODA::Dbn2D const & dbn)
    {
        xt::xarray<double> array({9}, 0);
        array(0) = dbn.sumW();
        array(1) = dbn.sumW2();
        array(2) = dbn.sumWX();
        array(3) = dbn.sumWX2();
        array(4) = dbn.sumWY();
        array(5) = dbn.sumWY2();
        array(6) = dbn.numEntries();
        // array(7) = 0.; // default value, no need
        // array(8) = 0.;
        return std::move(array);
    }

    xt::xarray<double> point1ToArray(YODA::Point1D const & pnt)
    {
        xt::xarray<double> array({3}, 0);
        array(0) = pnt.val(1);
        array(1) = pnt.errMinus(1);
        array(2) = pnt.errPlus(1);
        return std::move(array);
    }

    xt::xarray<double> point2ToArray(YODA::Point2D const & pnt)
    {
        xt::xarray<double> array({6}, 0);
        array(0) = pnt.val(1);
        array(1) = pnt.errMinus(1);
        array(2) = pnt.errPlus(1);
        array(3) = pnt.val(2);
        array(4) = pnt.errMinus(2);
        array(5) = pnt.errPlus(2);
        return std::move(array);
    }

    std::vector<std::string> get_ao_names(const std::vector<YODA::AnalysisObject*>& aos, std::string const & aotype)
    {
        std::vector<std::string> aonames;
        for (auto ao : aos) 
        {
            if(aotype != YODF5::NA && ao->type() != aotype) {
                continue;
            }
            if (! ao->path().ends_with( ']')) {
                aonames.push_back(ao->path());
            }
        }
        sort(aonames.begin(), aonames.end());
        return aonames;
    }


    std::vector<std::string> get_variations(const std::vector<YODA::AnalysisObject*>& aos, const std::string & hname)
    {
        std::vector<std::string> vars= {""};
        for (auto ao : aos) 
        {
            if (ao->path().starts_with(hname+"[") && ao->path().ends_with( ']'))
            {
                std::string s(ao->path());
                s.erase(s.begin(), s.begin()+hname.length()+1);
                s.erase(s.end()-1);
                vars.push_back(s);
            }
        }
        sort(vars.begin(), vars.end());
        return vars;
    }

    std::vector<std::string> mk_binids(std::vector<YODA::AnalysisObject*> const & aos) 
    {
        std::vector<std::string> binids;
        for (auto ao : aos) {
            auto hname = ao->path();
            if(hname[hname.size()-1] == ']') continue; // ignore systematic variations

            std::string base = replaceall(hname, "/", "|");
            std::vector<std::string> suffixes;
            auto aotype = ao->type();
            if (aotype == YODF5::H1D || aotype == YODF5::H2D || aotype == YODF5::P1D) {
                // for histo-1d, 2d and profile 1d, add total, overflow and underflow bins
                suffixes.push_back("T"); // total
                suffixes.push_back("O"); // overflow
                suffixes.push_back("U"); // underflow
            }
            if (aotype == YODF5::CNT) {
                suffixes.push_back("0");
            } else {
                int n_bins =num_bins(ao);
                for (int nb = 0; nb < n_bins; ++nb){
                    suffixes.push_back(std::to_string(nb));
                }
            }
            for(auto const & suffix: suffixes) {
                binids.push_back(base+"#"+suffix);
            }
        }
        return binids;
    }

    void  fill_dataset(
        H5Easy::File& file,
        std::unordered_map<std::string, YODA::AnalysisObject*>& aomap,
        std::string const & hname,
        std::vector<std::string> const & binids,
        std::vector<std::string> const & variations
        )
    {
        std::string _hname = replaceall(hname, "/", "|");
        // std::cout << hname << " --> " << _hname << std::endl;

        // all results are saved in a matrix (nbins, n_variations)
        // bin_idx saves the row index of this variable.
        std::vector<int> bin_idx;
        for(int i = 0; i < (int)binids.size(); i++){
            std::string binid_basename = binids.at(i).substr(0, binids.at(i).find("#"));
            if(binid_basename == _hname){
                bin_idx.push_back(i);
                // std::cout <<" find: " << binids.at(i) << " for " << _hname << std::endl;
            }
        }
        auto ao = aomap[hname];
        std::string aotype = ao->type();
        // H5Easy::dump(file, aotype+"/"+hname, bin_idx); // save the index
        //
        // std::cout << bin_idx.size() << " " << variations.size() << std::endl;

        long unsigned int n_bins = bin_idx.size();
        long unsigned int n_variations = variations.size();
        long unsigned int n_fields;

        std::vector<std::string> field_keys;

        if(aotype == YODF5::H1D){
            field_keys.assign({
                YODF5::SUMW, YODF5::SUMW2, YODF5::SUMWX, YODF5::SUMWX2, 
                YODF5::NUMENTRIES, YODF5::XMIN, YODF5::XMAX});
        } else if (aotype == YODF5::H2D) {
            field_keys.assign({
                YODF5::SUMW, YODF5::SUMW2, YODF5::SUMWX, YODF5::SUMWX2,
                YODF5::SUMWY, YODF5::SUMWY2,
                YODF5::SUMWXY,
                YODF5::NUMENTRIES, YODF5::XMIN, YODF5::XMAX,
                YODF5::YMIN, YODF5::YMAX});
        } else if (aotype == YODF5::P1D) {
            field_keys.assign({
                YODF5::SUMW, YODF5::SUMW2, YODF5::SUMWX, YODF5::SUMWX2,
                YODF5::SUMWY, YODF5::SUMWY2,
                YODF5::NUMENTRIES, YODF5::XMIN, YODF5::XMAX});
        } else if (aotype == YODF5::S1D) {
            field_keys.assign({
                YODF5::XVAL, YODF5::XERRM, YODF5::XERRP
            });
        } else if (aotype == YODF5::S2D) {
            field_keys.assign({
                YODF5::XVAL, YODF5::XERRM, YODF5::XERRP,
                YODF5::YVAL, YODF5::YERRM, YODF5::YERRP
            });
        } else if (aotype == YODF5::CNT) {
            field_keys.assign({YODF5::SUMW, YODF5::SUMW2, YODF5::NUMENTRIES});
        } else {
            fprintf(stderr, "ERROR: type %s is unknown\n", aotype);
            return ;
        }
        // for(auto& field: field_keys){
        //     std::cout << field << " ";
        // }
        // std::cout << std::endl;

        // std::cout << "bins: " << n_bins << " " << n_variations << std::endl;
        n_fields = field_keys.size();
        xt::xarray<double> temp({n_bins, n_variations, n_fields}, 0);

        std::vector<std::string> hist_ids;
        for(const auto& v: variations) {
            char name[521];
            if (v != ""){
                sprintf(name, "%s[%s]", hname.c_str(), v.c_str());
            } else {
                sprintf(name, "%s", hname.c_str());
            }
            std::string name_(name);
            hist_ids.push_back(name_);
        }

        // loop over all variations, starting from nominal
        // the output is a 2D array: [n_bins, n_variations]
        // n_bins for Histo1D includes overflow and underflow
        for(int col = 0; col < (int)hist_ids.size(); col++) {
            auto& hist_name = hist_ids.at(col);
            YODA::AnalysisObject* ao = nullptr;
            try {
                ao = aomap.at(hist_name);
            } catch(const std::out_of_range& oot){
                printf("hist %s is missing\n", hist_name.c_str());
                continue;
            }
            std::string this_type = ao->type();

            // printf("checking hist %s with type %s\n", hist_name.c_str(), this_type.c_str());

            if(this_type == YODF5::CNT){
                YODA::Counter* ao_pt = dynamic_cast<YODA::Counter*>(ao);
                xt::view(temp, 0, col, xt::all())  = dbn0ToArray<YODA::Counter>(ao_pt);
            } else if (this_type == YODF5::H1D) {
                YODA::Histo1D* ao_pt = dynamic_cast<YODA::Histo1D*>(ao);
                xt::view(temp, 0, col, xt::all()) = dbn1ToArray(ao_pt->totalDbn());
                xt::view(temp, 1, col, xt::all()) = dbn1ToArray(ao_pt->overflow());
                xt::view(temp, 2, col, xt::all()) = dbn1ToArray(ao_pt->underflow());
                // loop over all bins
                int islice = 3;
                for(int ibin = 0; ibin < n_bins-3; ibin++) {
                    xt::view(temp, islice++, col, xt::all()) = dbn1ToArray(ao_pt->bin(ibin));
                }
            } else if (this_type == YODF5::H2D) {
                YODA::Histo2D* ao_pt = dynamic_cast<YODA::Histo2D*>(ao);
                xt::view(temp, 0, col, xt::all()) = Hdbn2ToArray(ao_pt->totalDbn());
                xt::view(temp, 1, col, xt::all()) = 0.;
                xt::view(temp, 2, col, xt::all()) = 0.;
                // loop over all bins
                int islice = 3;
                for(int ibin = 0; ibin < n_bins-3; ibin++) {
                    xt::view(temp, islice++, col, xt::all()) = Hdbn2ToArray(ao_pt->bin(ibin));
                }
            } else if (this_type == YODF5::P1D) {
                YODA::Profile1D* ao_pt = dynamic_cast<YODA::Profile1D*>(ao);
                xt::view(temp, 0, col, xt::all()) = dbn2ToArray(ao_pt->totalDbn());
                xt::view(temp, 1, col, xt::all()) = dbn2ToArray(ao_pt->overflow());
                xt::view(temp, 2, col, xt::all()) = dbn2ToArray(ao_pt->underflow());
                // loop over all bins
                int islice = 3;
                for(int ibin = 0; ibin < n_bins-3; ibin++) {
                    xt::view(temp, islice++, col, xt::all()) = dbn1ToArray(ao_pt->bin(ibin));
                }

            } else if (this_type == YODF5::S1D) {
                auto ao_pt = dynamic_cast<YODA::Scatter1D*>(ao);
                for(int ip = 0; ip < n_bins; ip++) {
                    xt::view(temp, ip, col, xt::all()) = point1ToArray(ao_pt->point(ip));
                }
            } else if (this_type == YODF5::S2D) {
                auto ao_pt = dynamic_cast<YODA::Scatter2D*>(ao);
                for(int ip = 0; ip < n_bins; ip++) {
                    xt::view(temp, ip, col, xt::all()) = point2ToArray(ao_pt->point(ip));
                }
            } else {
                fprintf(stderr, "ERROR: type %s is unknown\n", this_type.c_str());
                return ;
            }
        }

        for(int ifield = 0; ifield < (int) field_keys.size(); ifield++) {
            auto& field = field_keys.at(ifield);
            // std::cout <<"updating: " << ifield << " " << field << std::endl;
            auto data = H5Easy::load<xt::xarray<double> > (file, field);
            auto v1 = xt::view(data, xt::keep(bin_idx), xt::all());
            auto v2 = xt::view(temp, xt::all(), xt::all(), xt::range(ifield, ifield+1));
            xt::xarray<double> v22(v2);
            v22.reshape({n_bins, n_variations});
            v1 = v22;
            // std::cout << std::endl;
            H5Easy::dump(file, field, data, H5Easy::DumpMode::Overwrite);
        }
    }

    void writeYoda2H5(std::vector<YODA::AnalysisObject*> const & aos, std::string const & outh5_name, int /*not used*/)
    {
        std::unordered_map<std::string, YODA::AnalysisObject*> aomap;
        for (auto ao : aos) aomap.insert({ao->path(), ao});

        // Determines sorting order --- in a robust application all hnames must be known
        auto const hnames = get_ao_names(aos, YODF5::NA);
        // for(const auto& h1_name: h1d_names) {
        //     std::cout << h1_name << std::endl;
        // }
        std::vector<std::string> const binids    = mk_binids(aos);
        std::cout << binids[0] << std::endl;
        auto const vars      = get_variations(aos, hnames[0]);

        std::cout << vars.size() << " variations" << std::endl;
        std::cout << binids.size() << " bins" << std::endl;


        H5Easy::File file(outh5_name, H5Easy::File::ReadWrite | H5Easy::File::Create| H5Easy::File::Truncate);
        H5Easy::dump(file, YODF5::BINID, binids);
        H5Easy::dump(file, YODF5::VARIATIONS, vars);
        std::cout << "hist name: " << hnames[0] << std::endl;
        std::cout << "binids : " << binids.size() << std::endl;

        std::vector<size_t> data_shape = {binids.size(), vars.size()};
        for(auto& column_name: YODF5::ALL_COLUMNS) {
            H5Easy::dump(file, column_name, xt::xarray<double>::from_shape(data_shape));
        }
        for(auto& hname: hnames) {
            fill_dataset(file, aomap, hname, binids, vars);
        }
    }

    };
};