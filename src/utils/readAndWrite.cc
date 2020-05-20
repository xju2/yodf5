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
    // H5Easy::dump(file, aotype+"/"+hname, bin_idx);
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
    using arr2d_t = xt::xarray<double>; 
    std::map<std::string, arr2d_t > data_map;
    for(auto field_key: field_keys) {
        data_map[field_key] = xt::xarray<double>::from_shape({n_bins, n_variations});
    }
	print_shape(data_map[YODF5::SUMW]);
	cout << data_map[YODF5::SUMW] << endl;
	data_map[YODF5::SUMW](0, 0) = 12.0;
	cout << data_map[YODF5::SUMW](0, 0) << endl;


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
    cout << "HERE? " << endl;
    for(int col = 0; col < (int)hist_ids.size(); col++) {
        auto& hist_name = hist_ids.at(col);
		printf("checking hist %s\n", hist_name.c_str());
		AnalysisObject* ao = nullptr;
		try {
        	ao = aomap.at(hist_name);
		} catch(const std::out_of_range& oot){
			printf("hist %s is missing\n", hist_name.c_str());
			continue;
		}
        std::string this_type = ao->type();

        if(this_type == YODF5::CNT){
			printf("processing %s %s\n", this_type.c_str(), ao->path().c_str());
			YODA::Counter* ao_pt = dynamic_cast<YODA::Counter*>(ao);
			if(ao_pt == 0) {
				printf("ERROR\n");
				continue;
			}
			printf("%.2f\n", ao_pt->sumW());
			printf("%.2f %.2f %.2f\n", ao_pt->sumW(), ao_pt->sumW2(), ao_pt->numEntries());
            data_map[YODF5::SUMW](0, col) = 1.0;
			printf("HERE\n");
            data_map[YODF5::SUMWY](0, col) = ao_pt->sumW2();
            data_map[YODF5::NUMENTRIES](0, col) = ao_pt->numEntries();
			printf("done with %s\n", this_type.c_str());
        } else if (this_type == YODF5::H1D) {
        } else if (this_type == YODF5::H2D) {
            n_fields = 12;
        } else if (this_type == YODF5::P1D) {
            n_fields = 9;
        } else if (this_type == YODF5::S1D) {
            n_fields = 3;
        } else if (this_type == YODF5::S2D) {
        } else if (this_type == YODF5::CNT) {
            n_fields = 3;
        } else {
            fprintf(stderr, "ERROR: type %s is unknown\n", this_type.c_str());
            return ;
        }
    }

}


 

int main(int argc, char** argv)
{
    int compression = 1;
    std::string filename = "Rivet.yoda";
    std::string outname = "examples.h5";
    bool help = false;
    int opt;
    while ((opt = getopt(argc, argv, "hf:c:")) != -1) {
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
    // HighFive::File file(outname, H5Easy::File::ReadWrite | H5Easy::File::Create| H5Easy::File::Truncate);
    // DataSet dataset = file.createDataSet<std::string>(YODF5::BINID, DataSpace::From(binids));
    // dataset.write(binids);
    // dataset = file.createDataSet<std::string>(YODF5::VARIATIONS, DataSpace::From(vars));
    // dataset.write(vars);
    std::cout << "hist name: " << hnames[0] << std::endl;
    std::cout << "binids : " << binids.size() << std::endl;
	fill_dataset(file, aomap, hnames[0], binids, vars);
    return 0;
    
    // file.create_group("Histo1D")
    // file.create_group("Histo2D")
    // file.create_group("Profile1D")
    // file.create_group("Counter")
    // file.create_group("Scatter1D")
    // file.create_group("Scatter2D")


    // DataSetCreateProps props;
    // props.add(Deflate(compression));
    // // This is the same chunking strategy h5py applies seemingly
    // props.add(Chunking(std::vector<hsize_t>{
    //             size_t(ceil(binids.size()/8.)),
    //             size_t(ceil(vars.size()/8.)),
    //             size_t(ceil(nFiles/8.))
    //             }));
    
    // //DataSetAccessProps cacheConfig;
    // //cacheConfig.add(Caching(1024, 1024, 0.5));

    // // Initial DS extends
    // size_t const NB = binids.size();
    // size_t const NV = vars.size();
    // // size_t const NF = nFiles;


    // // Dataset exttensible in 3rd dim, i.e. nfiles
    // file.createDataSet<double>("/Histo1D/sumW",       DataSpace( { NB, NV},{ NB, NV}), props);//, cacheConfig );
    // file.createDataSet<double>("/Histo1D/sumW2",      DataSpace( { NB, NV},{ NB, NV}), props);//, cacheConfig );
    // file.createDataSet<double>("/Histo1D/sumWX",      DataSpace( { NB, NV},{ NB, NV}), props);//, cacheConfig );
    // file.createDataSet<double>("/Histo1D/sumWX2",     DataSpace( { NB, NV},{ NB, NV}), props);//, cacheConfig );
    // file.createDataSet<double>("/Histo1D/numEntries", DataSpace( { NB, NV},{ NB, NV}), props);//, cacheConfig );

    // // Dummy here, we write the data of the first input file as many times as there are command line arguments
    // for (size_t fidx=0; fidx<NF; ++fidx)
    // {
    //     file.getDataSet("/Histo1D/sumW").select(       {0, 0, fidx}, {NB, NV, 1}).write(get_h1d_field(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::sumW,       &HistoBin1D::sumW));
    //     file.getDataSet("/Histo1D/sumW2").select(      {0, 0, fidx}, {NB, NV, 1}).write(get_h1d_field(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::sumW2,      &HistoBin1D::sumW2));
    //     file.getDataSet("/Histo1D/sumWX").select(      {0, 0, fidx}, {NB, NV, 1}).write(get_h1d_field(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::sumWX,      &HistoBin1D::sumWX));
    //     file.getDataSet("/Histo1D/sumWX2").select(     {0, 0, fidx}, {NB, NV, 1}).write(get_h1d_field(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::sumWX2,     &HistoBin1D::sumWX2));
    //     file.getDataSet("/Histo1D/numEntries").select( {0, 0, fidx}, {NB, NV, 1}).write(get_h1d_field(aomap, h1d_names, vars, NB, nbinmap, &Dbn1D::numEntries, &HistoBin1D::numEntries));
    // }
    
    // file.createDataSet("Histo1D/xMin",   get_edge(aomap, h1d_names, &HistoBin1D::xMin));
    // file.createDataSet("Histo1D/xMax",   get_edge(aomap, h1d_names, &HistoBin1D::xMax));
    
    // H5Easy::dump(file, "Histo1D/binids", binids);
    // H5Easy::dump(file, "Histo1D/names", h1d_names);


    

    return 0;
}
