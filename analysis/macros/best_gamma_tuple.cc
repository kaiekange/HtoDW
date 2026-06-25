// For each Ds candidate (identified by run / lumi / event / cand_id), keep only
// the row with the highest gamma_pt. A Ds candidate with no photon is a single
// row (gamma_pt = sentinel) and is kept unchanged. The output has one row per
// Ds candidate.
//
// The rows of one Ds candidate are stored contiguously, so we just walk the
// tree in order and, whenever the (run, lumi, event, cand_id) key changes,
// remember which row of the previous group had the largest gamma_pt.
//
// Usage:
//   root -l -b -q 'analysis/macros/best_gamma_tuple.cc("INPUT.root", "OUTPUT.root", "HKKpimu")'
// The third argument is the in-file directory that holds the Events and Info trees.

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <iostream>
#include <vector>

void best_gamma_tuple(const char* input, const char* output, const char* subdir)
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
    // First pass: read only the key columns and gamma_pt. For each group of
    // rows that share the same (run, lumi, event, cand_id), find the row with
    // the highest gamma_pt and flag it for keeping.
    // ---------------------------------------------------------------------
    UInt_t    runId    = 0;
    UInt_t    lumiId   = 0;
    ULong64_t eventId  = 0;
    Int_t     candId   = 0;
    Float_t   gammaPt  = 0;

    eventsIn->SetBranchStatus("*", false);
    eventsIn->SetBranchStatus("run_id",   true);
    eventsIn->SetBranchStatus("lumi_id",  true);
    eventsIn->SetBranchStatus("event_id", true);
    eventsIn->SetBranchStatus("cand_id",  true);
    eventsIn->SetBranchStatus("gamma_pt", true);

    eventsIn->SetBranchAddress("run_id",   &runId);
    eventsIn->SetBranchAddress("lumi_id",  &lumiId);
    eventsIn->SetBranchAddress("event_id", &eventId);
    eventsIn->SetBranchAddress("cand_id",  &candId);
    eventsIn->SetBranchAddress("gamma_pt", &gammaPt);

    std::vector<bool> keepRow(nRows, false);

    UInt_t    prevRun     = 0;
    UInt_t    prevLumi    = 0;
    ULong64_t prevEvent   = 0;
    Int_t     prevCand    = 0;
    Long64_t  bestRow     = 0;
    Float_t   bestGammaPt = 0;

    for (Long64_t row = 0; row < nRows; ++row)
    {
        eventsIn->GetEntry(row);

        bool sameGroup = (row > 0)
                      && runId   == prevRun
                      && lumiId  == prevLumi
                      && eventId == prevEvent
                      && candId  == prevCand;

        if (sameGroup)
        {
            // Still inside the current Ds candidate: update the best row.
            if (gammaPt > bestGammaPt)
            {
                bestGammaPt = gammaPt;
                bestRow     = row;
            }
        }
        else
        {
            // A new Ds candidate starts here. Commit the best row of the
            // previous one, then start tracking this group.
            if (row > 0)
                keepRow[bestRow] = true;

            bestGammaPt = gammaPt;
            bestRow     = row;
        }

        prevRun   = runId;
        prevLumi  = lumiId;
        prevEvent = eventId;
        prevCand  = candId;
    }

    // Commit the best row of the very last group.
    keepRow[bestRow] = true;

    Long64_t nSelected = 0;
    for (Long64_t row = 0; row < nRows; ++row)
        if (keepRow[row])
            ++nSelected;
    std::cout << "Selected rows (one per Ds cand): " << nSelected << std::endl;

    // ---------------------------------------------------------------------
    // Second pass: copy every branch, but write out only the flagged rows.
    // ---------------------------------------------------------------------
    eventsIn->SetBranchStatus("*", true);

    TFile* outputFile = TFile::Open(output, "RECREATE");
    TDirectory* outputDir = outputFile->mkdir(subdir);
    outputDir->cd();

    TTree* eventsOut = eventsIn->CloneTree(0);
    for (Long64_t row = 0; row < nRows; ++row)
    {
        if (!keepRow[row])
            continue;
        eventsIn->GetEntry(row);
        eventsOut->Fill();
    }
    eventsOut->Write("", TObject::kOverwrite);

    // The Info tree (bookkeeping) is copied unchanged.
    if (infoIn != nullptr)
    {
        TTree* infoOut = infoIn->CloneTree(-1, "fast");
        infoOut->Write("", TObject::kOverwrite);
    }

    outputFile->Close();
    inputFile->Close();

    std::cout << "Written: " << output << "  (Events = " << nSelected << ")" << std::endl;
}
