// Keep one row per event: the highest-pt Ds candidate, and within it the
// highest-pt gamma. A Ds candidate with no photon is a single row
// (has_gamma = false, gamma_pt = sentinel) and is kept as-is if it wins on Ds_pt.
//
// A Lumi_reweight branch (constant per sample, equivalent-luminosity weight) is
// added to the output. The Info tree is copied unchanged.
//
// Rows are stored grouped by event and, within an event, by Ds candidate, so a
// single pass suffices: track the best gamma row inside the current candidate,
// commit it against the event's best Ds_pt at each candidate boundary, and emit
// the event's winner at each event boundary.
//
// Every cfg has used the module label Hphipimu since 2026-08-18, so subdir and outsubdir
// are now the same string for new tuples. Both parameters are kept because tuples produced
// before that date still hold the signal trees under HphipimuDsW / HphipimuDsstW.
//
// Usage:
//   root -l -b -q 'best_cand_tuple.cc("INPUT.root", "OUTPUT.root", "Hphipimu", 6.3308e-06, "Hphipimu")'

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <iostream>
#include <limits>
#include <vector>

void best_cand_tuple(const char* input, const char* output, const char* subdir, Double_t lumi_reweight, const char* outsubdir)
{
    TFile* inputFile = TFile::Open(input);
    if (inputFile == nullptr || inputFile->IsZombie())
    {
        std::cerr << "Cannot open input file: " << input << std::endl;
        return;
    }

    TTree* eventsIn = (TTree*) inputFile->Get(TString(subdir) + "/Events");
    TTree* infoIn   = (TTree*) inputFile->Get(TString(subdir) + "/Info");
    if (eventsIn == nullptr)
    {
        std::cerr << "No Events tree found under " << subdir << std::endl;
        return;
    }

    Long64_t nRows = eventsIn->GetEntries();
    std::cout << "Events rows: " << nRows << std::endl;
    if (nRows == 0)
        return;

    // ---------------------------------------------------------------------
    // First pass: read only the key columns, Ds_pt and gamma_pt, and flag the
    // single row that survives for each event.
    // ---------------------------------------------------------------------
    UInt_t    runId   = 0;
    UInt_t    lumiId  = 0;
    ULong64_t eventId = 0;
    Int_t     candId  = 0;
    Float_t   dsPt    = 0;
    Float_t   gammaPt = 0;

    eventsIn->SetBranchStatus("*", false);
    eventsIn->SetBranchStatus("run_id",   true);
    eventsIn->SetBranchStatus("lumi_id",  true);
    eventsIn->SetBranchStatus("event_id", true);
    eventsIn->SetBranchStatus("cand_id",  true);
    eventsIn->SetBranchStatus("Ds_pt",    true);
    eventsIn->SetBranchStatus("gamma_pt", true);

    eventsIn->SetBranchAddress("run_id",   &runId);
    eventsIn->SetBranchAddress("lumi_id",  &lumiId);
    eventsIn->SetBranchAddress("event_id", &eventId);
    eventsIn->SetBranchAddress("cand_id",  &candId);
    eventsIn->SetBranchAddress("Ds_pt",    &dsPt);
    eventsIn->SetBranchAddress("gamma_pt", &gammaPt);

    std::vector<bool> keepRow(nRows, false);

    const Float_t lowest = -std::numeric_limits<Float_t>::max();

    // current Ds candidate
    Long64_t candBestRow     = 0;   // row of this candidate with the highest gamma_pt
    Float_t  candBestGammaPt = 0;
    Float_t  candDsPt        = 0;

    // current event
    Long64_t evSelRow    = -1;      // row that wins the event so far
    Float_t  evBestDsPt  = lowest;
    Long64_t nEvents     = 0;

    UInt_t    prevRun   = 0;
    UInt_t    prevLumi  = 0;
    ULong64_t prevEvent = 0;
    Int_t     prevCand  = 0;

    for (Long64_t row = 0; row < nRows; ++row)
    {
        eventsIn->GetEntry(row);

        bool sameEvent = (row > 0)
                      && runId   == prevRun
                      && lumiId  == prevLumi
                      && eventId == prevEvent;
        bool sameCand = sameEvent && candId == prevCand;

        if (row > 0 && !sameCand)
        {
            // The previous Ds candidate is complete: let it compete for the event.
            if (candDsPt > evBestDsPt)
            {
                evBestDsPt = candDsPt;
                evSelRow   = candBestRow;
            }

            if (!sameEvent)
            {
                // The previous event is complete: emit its winner.
                keepRow[evSelRow] = true;
                ++nEvents;
                evBestDsPt = lowest;
                evSelRow   = -1;
            }
        }

        if (row == 0 || !sameCand)
        {
            // Start tracking a new Ds candidate. Ds_pt is constant across its rows.
            candBestRow     = row;
            candBestGammaPt = gammaPt;
            candDsPt        = dsPt;
        }
        else if (gammaPt > candBestGammaPt)
        {
            candBestGammaPt = gammaPt;
            candBestRow     = row;
        }

        prevRun   = runId;
        prevLumi  = lumiId;
        prevEvent = eventId;
        prevCand  = candId;
    }

    // Close the final candidate and the final event.
    if (candDsPt > evBestDsPt)
    {
        evBestDsPt = candDsPt;
        evSelRow   = candBestRow;
    }
    keepRow[evSelRow] = true;
    ++nEvents;

    std::cout << "Selected rows (one per event): " << nEvents << std::endl;

    // ---------------------------------------------------------------------
    // Second pass: copy every branch of the flagged rows, plus Lumi_reweight.
    // ---------------------------------------------------------------------
    eventsIn->SetBranchStatus("*", true);

    TFile* outputFile = TFile::Open(output, "RECREATE");
    TDirectory* outputDir = outputFile->mkdir(outsubdir);
    outputDir->cd();

    TTree* eventsOut = eventsIn->CloneTree(0);

    Double_t lumiReweight = lumi_reweight;
    eventsOut->Branch("Lumi_reweight", &lumiReweight, "Lumi_reweight/D");

    for (Long64_t row = 0; row < nRows; ++row)
    {
        if (!keepRow[row])
            continue;
        eventsIn->GetEntry(row);
        eventsOut->Fill();
    }
    eventsOut->Write("", TObject::kOverwrite);

    if (infoIn != nullptr)
    {
        TTree* infoOut = infoIn->CloneTree(-1, "fast");
        infoOut->Write("", TObject::kOverwrite);
    }

    outputFile->Close();
    inputFile->Close();

    std::cout << "Written: " << output
              << "  (Events = " << nEvents
              << ", Lumi_reweight = " << lumi_reweight << ")" << std::endl;
}
