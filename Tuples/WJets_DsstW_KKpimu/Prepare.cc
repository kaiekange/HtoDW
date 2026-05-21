#include <map>
#include <vector>
#include <algorithm>

Int_t Prepare()
{
    TFile *fin = TFile::Open("root://maite.iihe.ac.be//store/user/kakang/H2DW/Analysis/MC/2022/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/WJets_DsstW_KKpimu/260511_105314/0000/WJets_DsstW_KKpimu.root");
    TTree *tree = fin->Get<TTree>("DsstWKKpimu/Events");

    // Pass 1: collect entries passing isolation cuts, grouped by event_id
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("event_id", 1);

    ULong64_t event_id;
    tree->SetBranchAddress("event_id", &event_id);

    std::map<ULong64_t, std::vector<Long64_t>> evtEntries;
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        evtEntries[event_id].push_back(i);
    }

    // Randomly select one candidate per event
    TRandom3 rng(42);
    std::vector<Long64_t> selected;
    selected.reserve(evtEntries.size());
    for (auto& [evtId, entries] : evtEntries) {
        int idx = (int)rng.Uniform(0, entries.size());
        selected.push_back(entries[idx]);
    }
    std::sort(selected.begin(), selected.end());

    // Pass 2: write selected entries
    tree->SetBranchStatus("*", 1);
    tree->ResetBranchAddresses();

    TFile *fout = new TFile("/eos/home-k/kakang/H2DW/Analysis/MVA/DsstW_KKpimu/Background.root", "RECREATE");
    TTree *tout = tree->CloneTree(0);
    for (Long64_t entry : selected) {
        tree->GetEntry(entry);
        tout->Fill();
    }

    fout->Write();
    fout->Close();
    fin->Close();

    return 0;
}
