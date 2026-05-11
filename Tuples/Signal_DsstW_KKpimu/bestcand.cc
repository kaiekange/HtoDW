#include "TFile.h"
#include "TTree.h"

Int_t bestcand()
{
    ROOT::EnableImplicitMT();
    TFile *myfile = TFile::Open("/eos/home-k/kakang/H2DW/Analysis/MC/2022/RecoAll_DsW_KKpimu/RecoAll_DsW_KKpimu.root");
    TTree *mytree = (TTree *)myfile->Get("RecoAllDsWKKpimu/Events");

    Int_t nGen = 225277;

    Int_t nEvents = mytree->GetEntries();
    Int_t nMatched = mytree->GetEntries("event_match");

    Int_t nMatch = 0;
    Int_t nDsMatch = 0;
    Int_t nDsNoMatch = 0;
    Int_t nEvtMatch = 0;
    Int_t nEvtNoMatch = 0;

    std::vector<Bool_t> *event_match = nullptr;
    std::vector<Bool_t> *Ds_match = nullptr;
    mytree->SetBranchAddress("event_match", &event_match);
    mytree->SetBranchAddress("Ds_match", &Ds_match);

    std::vector<Float_t> *Ds_FDxysig = nullptr;
    std::vector<Float_t> *Ds_IsoR04_RelIso = nullptr;
    std::vector<Float_t> *dR_Kp_Km = nullptr;
    std::vector<Float_t> *dR_Kp_pi = nullptr;
    std::vector<Float_t> *dR_Km_pi = nullptr;
    std::vector<Float_t> *dR_Kp_Ds = nullptr;
    std::vector<Float_t> *dR_Km_Ds = nullptr;
    std::vector<Float_t> *dR_pi_Ds = nullptr;
    mytree->SetBranchAddress("Ds_FDxysig", &Ds_FDxysig);
    mytree->SetBranchAddress("Ds_IsoR04_RelIso", &Ds_IsoR04_RelIso);
    mytree->SetBranchAddress("dR_Kp_Km", &dR_Kp_Km);
    mytree->SetBranchAddress("dR_Kp_pi", &dR_Kp_pi);
    mytree->SetBranchAddress("dR_Km_pi", &dR_Km_pi);
    mytree->SetBranchAddress("dR_Kp_Ds", &dR_Kp_Ds);
    mytree->SetBranchAddress("dR_Km_Ds", &dR_Km_Ds);
    mytree->SetBranchAddress("dR_pi_Ds", &dR_pi_Ds);
    
    std::vector<Float_t> *Ds_pt = nullptr;
    mytree->SetBranchAddress("Ds_pt", &Ds_pt);

    for (Int_t i = 0; i < nEvents; i++)
    {
        mytree->GetEntry(i);

        Int_t vecsize = event_match->size();
        if (vecsize < 1)
            continue;

        Float_t best_Ds_pt = -1.;
        Int_t idx_bestcand = -1;

        for (Int_t j = 0; j < vecsize; j++)
        {
            if(Ds_FDxysig->at(j) < 1)
                continue;
            if(Ds_IsoR04_RelIso->at(j) > 1)
                continue;
            if(dR_Kp_Km->at(j) > 0.4)
                continue;
            if(dR_Kp_pi->at(j) > 0.5)
                continue;
            if(dR_Km_pi->at(j) > 0.5)
                continue;
            if(dR_Kp_Ds->at(j) > 0.25)
                continue;
            if(dR_Km_Ds->at(j) > 0.25)
                continue;
            if(dR_pi_Ds->at(j) > 0.4)
                continue;

            if(event_match->at(j))
                nMatch++;

            if (Ds_pt->at(j) > best_Ds_pt)
            {
                best_Ds_pt = Ds_pt->at(j);
                idx_bestcand = j;
            }
        }
        if(idx_bestcand < 0)
            continue;

        if (Ds_match->at(idx_bestcand))
            nDsMatch++;
        else
            nDsNoMatch++;
        if (event_match->at(idx_bestcand))
            nEvtMatch++;
        else
            nEvtNoMatch++;
    }

    std::cout << "nMatch: " << nMatch << ", nDsMatch: " << nDsMatch << ", nDsNoMatch: " << nDsNoMatch << ", nEvtMatch: " << nEvtMatch << ", nEvtNoMatch: " << nEvtNoMatch << std::endl;

    return 0;
}