#include <boost/filesystem.hpp>

#include "data_tools/transforms.h"
#include "data_tools/benchmark.h"

#include "data_tools/navi_data.h"
#include "data_tools/std_data.h"
#include "data_tools/gsf_data.h"
#include <data_tools/csv_data.h>
#include <data_tools/xtf_data.h>
#include <data_tools/all_data.h>

#include "submaps_tools/submaps.hpp"
#include "submaps_tools/cxxopts.hpp"

#include "registration/utils_visualization.hpp"

using namespace std;

int main(int argc, char** argv) {
    string folder_str;
    string file_str;
    string type;

    cxxopts::Options options("example_reader", "Reads different mbes and sss file formats and saves them to a common format");
    options.add_options()
      ("help", "Print help")
      ("folder", "Input folder containing mbes files", cxxopts::value(folder_str))
      ("file", "Output file", cxxopts::value(file_str))
      ("type", "Type of data to read, options: all, xtf, navi, gsf", cxxopts::value(type));

    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        cout << options.help({ "", "Group" }) << endl;
        exit(0);
    }

    if (result.count("folder") == 0) {
        cout << "Please provide folder containing mbes or sss files..." << endl;
        exit(0);
    }
    if (result.count("type") == 0) {
        cout << "Please provide input type, options: all, xtf, navi, gsf" << endl;
        exit(0);
    }


    boost::filesystem::path folder(folder_str);
    boost::filesystem::path pings_path = folder / "mbes_pings.cereal";
    boost::filesystem::path submaps_path = folder / "submaps.cereal";

    // if side scan, read and save, then return
    if (type == "xtf") {
        cout << "Not much to do by now" << endl;
        return 0;
    }

    // otherwise we have multibeam data, read and save
    std_data::mbes_ping::PingsT std_pings;
    if (type == "gsf") {
        gsf_data::gsf_mbes_ping::PingsT pings = std_data::parse_folder<gsf_data::gsf_mbes_ping>(folder);
        std::stable_sort(pings.begin(), pings.end(), [](const gsf_data::gsf_mbes_ping& ping1, const gsf_data::gsf_mbes_ping& ping2) {
            return ping1.time_stamp_ < ping2.time_stamp_;
        });
        std_pings = gsf_data::convert_pings(pings);
    }
    else if (type == "all") {
        all_data::all_nav_entry::EntriesT entries = std_data::parse_folder<all_data::all_nav_entry>(folder);
        all_data::all_mbes_ping::PingsT pings = std_data::parse_folder<all_data::all_mbes_ping>(folder);
        all_data::all_nav_attitude::EntriesT attitudes = std_data::parse_folder<all_data::all_nav_attitude>(folder);
        std_pings = all_data::convert_matched_entries(pings, entries);
        std_pings = all_data::match_attitude(std_pings, attitudes);
    }
    else if (type == "navi") {
        std_pings = std_data::parse_folder<std_data::mbes_ping>(folder / "Pings");
        std_data::nav_entry::EntriesT entries = std_data::parse_folder<std_data::nav_entry>(folder / "NavUTM");
        navi_data::match_timestamps(std_pings, entries);
    }
    else {
        cout << "Type " << type << " is not supported!" << endl;
        return 0;
    }

    // write to disk
    std_data::write_data<std_data::mbes_ping::PingsT>(std_pings, pings_path);


    return 0;
}
