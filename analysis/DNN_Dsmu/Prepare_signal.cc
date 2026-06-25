#include <map>
#include <vector>
#include <algorithm>

// Filter rows with Ds_match && mu_match, pick one per event,
// fill into tout (must have same branch structure).
void ProcessFile(const char *fname, const char *tree_name,
                 TTree *tout, TRandom3 &rng, Long64_t &n_selected)
{
    TFile *fin = TFile::Open(fname);
    TTree *tin = fin->Get<TTree>(tree_name);

    // Pass 1: collect qualifying entries grouped by event_id
    tin->SetBranchStatus("*", 0);
    tin->SetBranchStatus("event_id", 1);
    tin->SetBranchStatus("Ds_match", 1);
    tin->SetBranchStatus("mu_match", 1);

    ULong64_t event_id;
    Bool_t    Ds_match, mu_match;
    tin->SetBranchAddress("event_id", &event_id);
    tin->SetBranchAddress("Ds_match", &Ds_match);
    tin->SetBranchAddress("mu_match", &mu_match);

    std::map<ULong64_t, std::vector<Long64_t>> evtEntries;
    Long64_t nEntries = tin->GetEntries();
    printf("  Input entries: %lld\n", nEntries);

    for (Long64_t i = 0; i < nEntries; i++) {
        tin->GetEntry(i);
        if (Ds_match && mu_match)
            evtEntries[event_id].push_back(i);
    }
    printf("  Events with Ds_match && mu_match: %zu\n", evtEntries.size());

    // Randomly pick one candidate per event
    std::vector<Long64_t> selected;
    selected.reserve(evtEntries.size());
    for (auto &[evtId, entries] : evtEntries) {
        int idx = (int)rng.Uniform(0, entries.size());
        selected.push_back(entries[idx]);
    }
    std::sort(selected.begin(), selected.end());

    // Pass 2: point tin's branches at tout's memory, fill selected entries
    tin->SetBranchStatus("*", 1);
    tin->ResetBranchAddresses();
    tin->CopyAddresses(tout);   // tin reads into tout's branch buffers

    for (Long64_t entry : selected) {
        tin->GetEntry(entry);
        tout->Fill();
    }

    n_selected += (Long64_t)selected.size();
    tin->ResetBranchAddresses();
    fin->Close();
}

Int_t Prepare_signal()
{
    const char *dsstw_file = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Signal_DsstW_H_KKpimu/Signal_DsstW_H_KKpimu_bestgamma.root";
    const char *dsw_file   = "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Signal_H_KKpimu/Signal_H_KKpimu_bestgamma.root";
    const char *out_path   = "/eos/home-k/kakang/H2DW/Analysis/MVA/Dsmu/Signal.root";

    TRandom3 rng(42);

    // Create output file, clone branch structure from DsstW (same as DsW)
    TFile *f_templ = TFile::Open(dsstw_file);
    TTree *t_templ = f_templ->Get<TTree>("HKKpimuDsstW/Events");
    TFile *fout    = new TFile(out_path, "RECREATE");
    TTree *tout    = t_templ->CloneTree(0);
    tout->SetName("Events");
    f_templ->Close();   // closes template; tout still lives in fout

    Long64_t n_total = 0;

    printf("Processing DsstW signal...\n");
    ProcessFile(dsstw_file, "HKKpimuDsstW/Events", tout, rng, n_total);

    printf("Processing DsW signal...\n");
    ProcessFile(dsw_file,   "HKKpimuDsW/Events",   tout, rng, n_total);

    fout->cd();
    fout->Write();
    fout->Close();

    printf("\nTotal selected: %lld\n", n_total);
    printf("Output: %s\n", out_path);
    return 0;
}
