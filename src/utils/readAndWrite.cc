#include "YODA/ReaderYODA.h"
#include "YODA/YODA.h" // include all analysis objects

#include <iostream>
#include <regex>
#include <unordered_map>
#include <functional>
#include <unistd.h>
#include <fstream>

#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"

#define H5_USE_XTENSOR
#include <highfive/H5Easy.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

#include "yodf5/const_keynames.hpp"
// #include "yodf5/h5utils.hpp"

using namespace YODA;
using namespace std;
// using namespace Eigen;
using namespace HighFive;


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

// for each bin in Histo1D
xt::xarray<double> dbn1ToArray(YODA::HistoBin1D const & dbn){
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
xt::xarray<double> Hdbn2ToArray(YODA::HistoBin2D const & dbn){
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

xt::xarray<double> point1ToArray(YODA::Point1D const & pnt) {
	xt::xarray<double> array({3}, 0);
	array(0) = pnt.val(1);
	array(1) = pnt.errMinus(1);
	array(2) = pnt.errPlus(1);
	return std::move(array);
}

xt::xarray<double> point2ToArray(YODA::Point2D const & pnt) {
	xt::xarray<double> array({6}, 0);
	array(0) = pnt.val(1);
	array(1) = pnt.errMinus(1);
	array(2) = pnt.errPlus(1);
	array(3) = pnt.val(2);
	array(4) = pnt.errMinus(2);
	array(5) = pnt.errPlus(2);
	return std::move(array);
}

unordered_map< string, unordered_map< string, int> > myyodaindex(istream& stream_) {

    unordered_map<string, unordered_map< string, int>> hmap;
    hmap.insert({"Histo1D",   unordered_map< string, int>()});
    hmap.insert({"Histo2D",   unordered_map< string, int>()});
    hmap.insert({"Profile1D", unordered_map< string, int>()});
    hmap.insert({"Profile2D", unordered_map< string, int>()});
    hmap.insert({"Scatter1D", unordered_map< string, int>()});
    hmap.insert({"Scatter2D", unordered_map< string, int>()});
    hmap.insert({"Scatter3D", unordered_map< string, int>()});
    hmap.insert({"Counter",   unordered_map< string, int>()});

    //#ifdef HAVE_LIBZ
    //// NB. zstr auto-detects if file is deflated or plain-text
    //zstr::istream stream(stream_);
    //#else
    istream& stream = stream_;
    //#endif
    enum Context
    {
        NONE, SCATTER1D, SCATTER2D, SCATTER3D, COUNTER, HISTO1D, HISTO2D,
        PROFILE1D, PROFILE2D
    };


    /// State of the parser: line number, line, parser context, and pointer(s) to the object currently being assembled
    unsigned int nline = 0;
    string s;
    Context context = NONE;
    std::string curpath="";
    int nbins=0;

    // Loop over all lines of the input file
    bool in_anns = false;
    string fmt = "1";
    while (std::getline(stream, s))
    {
        nline += 1;
        if (!in_anns)
        {
            Utils::itrim(s);
            if (s.empty()) continue;
            if (s.find("#") == 0 && s.find("BEGIN") == string::npos && s.find("END") == string::npos) continue;
        }
        if (context == NONE)
        {
            if (s.find("BEGIN ") == string::npos)
            {
                stringstream ss;
                ss << "Unexpected line in YODA format parsing when BEGIN expected: '" << s << "' on line " << nline;
                throw ReadError(ss.str());
            }
            while (s.find("#") == 0) s = Utils::trim(s.substr(1));
            vector<string> parts;
            istringstream iss(s); string tmp;
            while (iss >> tmp) parts.push_back(tmp);

            if (parts.size() < 2 || parts[0] != "BEGIN")
            {
                stringstream ss;
                ss << "Unexpected BEGIN line structure when BEGIN expected: '" << s << "' on line " << nline;
                throw ReadError(ss.str());
            }
            const string ctxstr = parts[1];
            curpath = (parts.size() >= 3) ? parts[2] : "";
            nbins=0;
            if      (Utils::startswith(ctxstr, "YODA_COUNTER"))   {context = COUNTER;  }
            else if (Utils::startswith(ctxstr, "YODA_SCATTER1D")) {context = SCATTER1D;}
            else if (Utils::startswith(ctxstr, "YODA_SCATTER2D")) {context = SCATTER2D;}
            else if (Utils::startswith(ctxstr, "YODA_SCATTER3D")) {context = SCATTER3D;}
            else if (Utils::startswith(ctxstr, "YODA_HISTO1D"))   {context = HISTO1D;  }
            else if (Utils::startswith(ctxstr, "YODA_HISTO2D"))   {context = HISTO2D;  }
            else if (Utils::startswith(ctxstr, "YODA_PROFILE1D")) {context = PROFILE1D;}
            else if (Utils::startswith(ctxstr, "YODA_PROFILE2D")) {context = PROFILE2D;}
            const size_t vpos = ctxstr.find_last_of("V");
            fmt = vpos != string::npos ? ctxstr.substr(vpos+1) : "1";
            if (fmt != "1") in_anns = true;
        }
        else
        { //< not a BEGIN line
            if (s.find("BEGIN ") != string::npos) throw ReadError("Unexpected BEGIN line in YODA format parsing before ending current BEGIN..END block");
            // FINISHING THE CURRENT CONTEXT
            if (s.find("END ") != string::npos)
            {
                if      (context == HISTO1D)   {hmap["Histo1D"]  .insert({curpath, nbins});}
                else if (context == HISTO2D)   {hmap["Histo2D"]  .insert({curpath, nbins});}
                else if (context == PROFILE1D) {hmap["Profile1D"].insert({curpath, nbins});}
                else if (context == PROFILE2D) {hmap["Profile2D"].insert({curpath, nbins});}
                else if (context == SCATTER1D) {hmap["Scatter1D"].insert({curpath, nbins});}
                else if (context == SCATTER2D) {hmap["Scatter2D"].insert({curpath, nbins});}
                else if (context == SCATTER3D) {hmap["Scatter3D"].insert({curpath, nbins});}
                else if (context == COUNTER)   {hmap["Counter"]  .insert({curpath, nbins});}
                in_anns = false;
                context = NONE;
                continue;
            }
            // ANNOTATIONS PARSING
            if (fmt == "1")
            {
                // First convert to one-key-per-line YAML syntax
                const size_t ieq = s.find("=");
                if (ieq != string::npos) s.replace(ieq, 1, ": ");
                // Special-case treatment for syntax clashes
                const size_t icost = s.find(": *");
                if (icost != string::npos)
                {
                    s.replace(icost, 1, ": '*");
                    s += "'";
                }
                // Store reformatted annotation
                const size_t ico = s.find(":");
                if (ico != string::npos) continue;
            }
            else if (in_anns)
            {
                if (s == "---") in_anns = false;
                continue;
            }

            if ( (context == HISTO1D) || (context==HISTO2D) || (context==PROFILE1D) || (context==PROFILE2D))
            {
                if (s.find("Total") != string::npos || s.find("Underflow") != string::npos || s.find("Overflow") != string::npos) continue;
                nbins++;
            }
            else if ( (context==SCATTER1D)  || (context==SCATTER2D) || (context==SCATTER3D))
            {
                nbins++;
            }
        }
    }
    return hmap;
}



// pre C++20
inline bool ends_with(string const & value, string const & ending)
{
    if (ending.size() > value.size()) return false;
    return equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline vector<string> get_ao_names(const vector<AnalysisObject*>& aos, string const & aotype)
{
    vector<string> aonames;
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


inline vector<string> get_variations(const vector<AnalysisObject*>& aos, const string & hname)
{
    vector<string> vars= {""};
    for (auto ao : aos) 
    {
        if (ao->path().starts_with(hname+"[") && ao->path().ends_with( ']'))
        {
            string s(ao->path());
            s.erase(s.begin(), s.begin()+hname.length()+1);
            s.erase(s.end()-1);
            vars.push_back(s);
        }
    }
    sort(vars.begin(), vars.end());
    return vars;
}

unordered_map<string, int> mk_nbinmap(unordered_map<string, AnalysisObject* > aomap, vector<string> const & hnames)
{
    unordered_map<string, int> nbin_map;
    for (auto hn : hnames) {
        nbin_map.insert({hn, dynamic_cast<Histo1D*>(aomap[hn])->numBins()});
    }
    return nbin_map;
}

string replaceall(const string& in_str, const string& source, const string& target) {
    std::string::size_type pos;
    std::string results(in_str);
    while ( (pos = results.find(source)) != std::string::npos) {
        results.replace(pos, source.size(), target);
    }
    return results;
}

int num_bins(AnalysisObject* ao){
    string aotype = ao->type();
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

inline vector<string> mk_binids(vector<AnalysisObject*> const & aos) 
{
    vector<string> binids;
    for (auto ao : aos) {
        auto hname = ao->path();
        if(hname[hname.size()-1] == ']') continue; // ignore systematic variations

        std::string base = replaceall(hname, "/", "|");
        vector<string> suffixes;
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
            int n_bins = num_bins(ao);
            for (int nb = 0; nb < n_bins; ++nb){
                suffixes.push_back(to_string(nb));
            }
        }
        for(auto const & suffix: suffixes) {
            binids.push_back(base+"#"+suffix);
        }
    }
    return binids;
}



inline vector<string> mk_binids(unordered_map<string, AnalysisObject* > aomap, vector<string> const & hnames) 
{
    vector<string> binids;
    for (auto hn : hnames) 
    {
        binids.push_back(hn+"#T");
        binids.push_back(hn+"#O");
        binids.push_back(hn+"#U");
        for (int nb=0;nb<dynamic_cast<Histo1D*>(aomap[hn])->numBins();++nb) binids.push_back(hn+"#" + to_string(nb));
    }
    return binids;
}


inline vector<double> get_edge(unordered_map<string, AnalysisObject* > aomap, vector<string> const & hnames, double (HistoBin1D::*function)() const) 
{
    vector<double> edges;
    for (auto hn : hnames) 
    {
        edges.push_back(0.0);
        edges.push_back(0.0);
        edges.push_back(0.0);
        auto h = dynamic_cast<Histo1D*>(aomap[hn]);
        for (int nb=0;nb< h->numBins();++nb) edges.push_back((h->bin(nb).*function)());
    }
    return edges;
}

void  fill_dataset(
        H5Easy::File& file,
        unordered_map<std::string, AnalysisObject*>& aomap,
        std::string const & hname,
        std::vector<std::string> const & binids,
        std::vector<std::string> const & variations
        ){
    std::string _hname = replaceall(hname, "/", "|");
	cout << hname << " --> " << _hname << endl;

    // all results are saved in a matrix (nbins, n_variations)
    // bin_idx saves the row index of this variable.
    std::vector<int> bin_idx;
    for(int i = 0; i < (int)binids.size(); i++){
		string binid_basename = binids.at(i).substr(0, binids.at(i).find("#"));
        if(binid_basename == _hname){
            bin_idx.push_back(i);
			cout <<" find: " << binids.at(i) << " for " << _hname << endl;
        }
    }
    auto ao = aomap[hname];
    std::string aotype = ao->type();
    // H5Easy::dump(file, aotype+"/"+hname, bin_idx); // save the index
	//
	cout << bin_idx.size() << " " << variations.size() << endl;

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
	for(auto& field: field_keys){
		cout << field << " ";
	}
	cout << endl;

	cout << "bins: " << n_bins << " " << n_variations << endl;
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
		string name_(name);
        hist_ids.push_back(name_);
    }

	// loop over all variations, starting from nominal
	// the output is a 2D array: [n_bins, n_variations]
	// n_bins for Histo1D includes overflow and underflow
    for(int col = 0; col < (int)hist_ids.size(); col++) {
        auto& hist_name = hist_ids.at(col);
		AnalysisObject* ao = nullptr;
		try {
        	ao = aomap.at(hist_name);
		} catch(const std::out_of_range& oot){
			printf("hist %s is missing\n", hist_name.c_str());
			continue;
		}
        std::string this_type = ao->type();

		printf("checking hist %s with type %s\n", hist_name.c_str(), this_type.c_str());

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
			cout << "HERE1 " << endl;
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
		cout <<"updating: " << ifield << " " << field << endl;
		auto data = H5Easy::load<xt::xarray<double> > (file, field);
		auto v1 = xt::view(data, xt::keep(bin_idx), xt::all());
		auto v2 = xt::view(temp, xt::all(), xt::all(), xt::range(ifield, ifield+1));
		xt::xarray<double> v22(v2);
		v22.reshape({n_bins, n_variations});
		/***
		// xt::view(data, xt::keep(bin_idx), xt::all()) = xt::view(temp, xt::all(), xt::all(), xt::range(ifield, ifield+1));
		cout << "v1: " ;
		for(int id=0; id < (int) v1.dimension(); id++) {
			cout  << v1.shape(id) << ", ";
		}
		cout << endl;
		cout << "v2: " ;
		for(int id=0; id < (int) v2.dimension(); id++) {
			cout  << v2.shape(id) << ", ";
		}
		**/
		v1 = v22;
		cout << endl;
		H5Easy::dump(file, field, data, H5Easy::DumpMode::Overwrite);
	}
}


 

int main(int argc, char** argv)
{
    int compression = 1;
    std::string filename = "Rivet.yoda";
    std::string outname = "examples.h5";
    bool help = false;
    int opt;
    while ((opt = getopt(argc, argv, "hf:c:o:")) != -1) {
        switch(opt) {
            case 'f':
                filename = optarg;
                break;
            case 'c':
                compression = atoi(optarg);
                break;
            case 'o':
                outname = optarg;
                break;
            case 'h':
                help = true;
            default:
                fprintf(stderr, "Usage: %s [-h] [-f FILENAME] [-c COMPRESSION]\n", argv[0]);
                if (help) {
                    printf("    -f FILENAME: input yoda file. Default is \"Rivet.yoda\"\n");
                    printf("    -c COMPRESSION: compression level. Default is 1\n");
                    printf("    -h HELP: print help info\n");
                }
                exit(EXIT_FAILURE);
        }
    }
    auto aos = ReaderYODA::create().read(filename);
    unordered_map<string, AnalysisObject*> aomap;
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


    std::ifstream instream;
    instream.open(filename);
    // This is a map of AOtype : {hname: nbins} --- runs 10x faster than read -- good for initial survey when
    // running over many files
    auto hmap = myyodaindex(instream);
    instream.close();
    auto nbinmap = hmap["Histo1D"];


    H5Easy::File file(outname, H5Easy::File::ReadWrite | H5Easy::File::Create| H5Easy::File::Truncate);
    H5Easy::dump(file, YODF5::BINID, binids);
    H5Easy::dump(file, YODF5::VARIATIONS, vars);
    std::cout << "hist name: " << hnames[0] << std::endl;
    std::cout << "binids : " << binids.size() << std::endl;

	vector<size_t> data_shape = {binids.size(), vars.size()};
	for(auto& column_name: YODF5::ALL_COLUMNS) {
		H5Easy::dump(file, column_name, xt::xarray<double>::from_shape(data_shape));
	}
	for(auto& hname: hnames) {
		fill_dataset(file, aomap, hname, binids, vars);
	}

    return 0;
}
