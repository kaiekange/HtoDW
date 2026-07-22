#include <map>
#include <vector>
#include <algorithm>

Int_t Prepare()
{
    TFile *fin = TFile::Open("/eos/home-k/kakang/H2DW/Analysis/MC/2022/WJets_H_phipimu/WJets_H_phipimu_bestgamma.root");
    TTree *tree = fin->Get<TTree>("Hphipimu/Events");

    // Pass 1: group all entries by event_id
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("event_id", 1);

    ULong64_t event_id;
    tree->SetBranchAddress("event_id", &event_id);

    std::map<ULong64_t, std::vector<Long64_t>> evtEntries;
    Long64_t nEntries = tree->GetEntries();
    printf("Total entries: %lld\n", nEntries);

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        evtEntries[event_id].push_back(i);
    }
    printf("Unique events: %zu\n", evtEntries.size());

    // Randomly select one candidate per event
    TRandom3 rng(42);
    std::vector<Long64_t> selected;
    selected.reserve(evtEntries.size());
    for (auto& [evtId, entries] : evtEntries) {
        int idx = (int)rng.Uniform(0, entries.size());
        selected.push_back(entries[idx]);
    }
    std::sort(selected.begin(), selected.end());
    printf("Selected: %zu\n", selected.size());

    // Pass 2: write selected entries to EOS
    tree->SetBranchStatus("*", 1);
    tree->ResetBranchAddresses();

    TFile *fout = new TFile("/eos/home-k/kakang/H2DW/Analysis/MVA/phipimu/WJets_Background.root", "RECREATE");
    TTree *tout = tree->CloneTree(0);

    for (Long64_t entry : selected) {
        tree->GetEntry(entry);
        tout->Fill();
    }

    fout->Write();
    fout->Close();
    fin->Close();

    printf("Done. Output: /eos/home-k/kakang/H2DW/Analysis/MVA/phipimu/WJets_Background.root\n");
    return 0;
}
