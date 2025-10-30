#include <ROOT/RDataFrame.hxx>
#include <TString.h>
#include <vector>
#include <iostream>

int MergeFile() {
    ROOT::EnableImplicitMT();
   
    std::set<int> exclude = {237, 270, 276, 277, 278, 356, 357, 367};

    std::vector<std::string> files;
    for (int i = 1; i <= 500; ++i) {
        if (exclude.count(i)) continue;
        std::string filename = "/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/PVStudy/output_" + std::to_string(i) + ".root";
        files.push_back(filename);
    }


    ROOT::RDataFrame df("PVStudy/Events", files);

    df.Snapshot("Events", "/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/PVStudy_noAP.root");

    return 0;
}
