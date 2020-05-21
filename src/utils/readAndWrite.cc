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
#include "yodf5/h5utils.hpp"

using namespace YODA;
using namespace std;
// using namespace Eigen;
using namespace HighFive;



unordered_map< string, unordered_map< string, int> > myyodaindex(istream& stream_) 
{

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






unordered_map<string, int> mk_nbinmap(unordered_map<string, AnalysisObject* > aomap, vector<string> const & hnames)
{
    unordered_map<string, int> nbin_map;
    for (auto hn : hnames) {
        nbin_map.insert({hn, dynamic_cast<Histo1D*>(aomap[hn])->numBins()});
    }
    return nbin_map;
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
    YODF5::H5Utils::writeYoda2H5(aos, outname);

    return 0;
}
