#include "TFile.h"
#include "TTree.h"

Int_t PrepareSB()
{
    TFile *myfile = TFile::Open("/eos/home-k/kakang/H2DW/Analysis/MC/2022/RecoAll_DsstW_KKpimu/RecoAll_DsstW_KKpimu.root");
    TTree *mytree = (TTree *)myfile->Get("RecoAllDsstWKKpimu/Events");

    Int_t nEvents = mytree->GetEntries();

    std::vector<Bool_t> *Ds_match_vec = nullptr;
    std::vector<Bool_t> *Kp_match_vec = nullptr;
    std::vector<Bool_t> *Km_match_vec = nullptr;
    std::vector<Bool_t> *pi_match_vec = nullptr;

    std::vector<Float_t> *Kp_eta_vec = nullptr;
    std::vector<Float_t> *Kp_phi_vec = nullptr;
    std::vector<Float_t> *Kp_pt_vec = nullptr;
    std::vector<Float_t> *Kp_p_vec = nullptr;
    
    std::vector<Float_t> *Km_eta_vec = nullptr;
    std::vector<Float_t> *Km_phi_vec = nullptr;
    std::vector<Float_t> *Km_pt_vec = nullptr;
    std::vector<Float_t> *Km_p_vec = nullptr;

    std::vector<Float_t> *pi_eta_vec = nullptr;
    std::vector<Float_t> *pi_phi_vec = nullptr;
    std::vector<Float_t> *pi_pt_vec = nullptr;
    std::vector<Float_t> *pi_p_vec = nullptr;

    std::vector<Float_t> *KK_eta_vec = nullptr;
    std::vector<Float_t> *KK_phi_vec = nullptr;
    std::vector<Float_t> *KK_pt_vec = nullptr;
    std::vector<Float_t> *KK_p_vec = nullptr;

    std::vector<Float_t> *Kppi_eta_vec = nullptr;
    std::vector<Float_t> *Kppi_phi_vec = nullptr;
    std::vector<Float_t> *Kppi_pt_vec = nullptr;
    std::vector<Float_t> *Kppi_p_vec = nullptr;

    std::vector<Float_t> *Kmpi_eta_vec = nullptr;
    std::vector<Float_t> *Kmpi_phi_vec = nullptr;
    std::vector<Float_t> *Kmpi_pt_vec = nullptr;
    std::vector<Float_t> *Kmpi_p_vec = nullptr;

    std::vector<Float_t> *Ds_eta_vec = nullptr;
    std::vector<Float_t> *Ds_phi_vec = nullptr;
    std::vector<Float_t> *Ds_pt_vec = nullptr;
    std::vector<Float_t> *Ds_p_vec = nullptr;

    std::vector<Float_t> *Ds_prob_vec = nullptr;
    std::vector<Float_t> *dira_vec = nullptr;

    std::vector<Float_t> *dR_Kp_Km_vec = nullptr;
    std::vector<Float_t> *dR_Kp_pi_vec = nullptr;
    std::vector<Float_t> *dR_Km_pi_vec = nullptr;
    std::vector<Float_t> *dR_Kp_Ds_vec = nullptr;
    std::vector<Float_t> *dR_Km_Ds_vec = nullptr;
    std::vector<Float_t> *dR_pi_Ds_vec = nullptr;

    std::vector<Float_t> *Kp_ipsig_vec = nullptr;
    std::vector<Float_t> *Kp_ipxysig_vec = nullptr;
    std::vector<Float_t> *Km_ipsig_vec = nullptr;
    std::vector<Float_t> *Km_ipxysig_vec = nullptr;
    std::vector<Float_t> *pi_ipsig_vec = nullptr;
    std::vector<Float_t> *pi_ipxysig_vec = nullptr;
    std::vector<Float_t> *Ds_ipsig_vec = nullptr;
    std::vector<Float_t> *Ds_ipxysig_vec = nullptr;

    std::vector<Float_t> *Ds_FDsig_vec = nullptr;
    std::vector<Float_t> *Ds_FDxysig_vec = nullptr;
    std::vector<Float_t> *Ds_IsoR04_RelIso_vec = nullptr;

    mytree->SetBranchAddress("Ds_match", &Ds_match_vec);
    mytree->SetBranchAddress("Kp_match", &Kp_match_vec);
    mytree->SetBranchAddress("Km_match", &Km_match_vec);
    mytree->SetBranchAddress("pi_match", &pi_match_vec);

    mytree->SetBranchAddress("Kp_eta", &Kp_eta_vec);
    mytree->SetBranchAddress("Kp_phi", &Kp_phi_vec);
    mytree->SetBranchAddress("Kp_pt", &Kp_pt_vec);
    mytree->SetBranchAddress("Kp_p", &Kp_p_vec);
    
    mytree->SetBranchAddress("Km_eta", &Km_eta_vec);
    mytree->SetBranchAddress("Km_phi", &Km_phi_vec);
    mytree->SetBranchAddress("Km_pt", &Km_pt_vec);
    mytree->SetBranchAddress("Km_p", &Km_p_vec);

    mytree->SetBranchAddress("pi_eta", &pi_eta_vec);
    mytree->SetBranchAddress("pi_phi", &pi_phi_vec);
    mytree->SetBranchAddress("pi_pt", &pi_pt_vec);
    mytree->SetBranchAddress("pi_p", &pi_p_vec);

    mytree->SetBranchAddress("KK_eta", &KK_eta_vec);
    mytree->SetBranchAddress("KK_phi", &KK_phi_vec);
    mytree->SetBranchAddress("KK_pt", &KK_pt_vec);
    mytree->SetBranchAddress("KK_p", &KK_p_vec);

    mytree->SetBranchAddress("Kppi_eta", &Kppi_eta_vec);
    mytree->SetBranchAddress("Kppi_phi", &Kppi_phi_vec);
    mytree->SetBranchAddress("Kppi_pt", &Kppi_pt_vec);
    mytree->SetBranchAddress("Kppi_p", &Kppi_p_vec);

    mytree->SetBranchAddress("Kmpi_eta", &Kmpi_eta_vec);
    mytree->SetBranchAddress("Kmpi_phi", &Kmpi_phi_vec);
    mytree->SetBranchAddress("Kmpi_pt", &Kmpi_pt_vec);
    mytree->SetBranchAddress("Kmpi_p", &Kmpi_p_vec);

    mytree->SetBranchAddress("Ds_eta", &Ds_eta_vec);
    mytree->SetBranchAddress("Ds_phi", &Ds_phi_vec);
    mytree->SetBranchAddress("Ds_pt", &Ds_pt_vec);
    mytree->SetBranchAddress("Ds_p", &Ds_p_vec);

    mytree->SetBranchAddress("Ds_prob", &Ds_prob_vec);
    mytree->SetBranchAddress("dira", &dira_vec);

    mytree->SetBranchAddress("dR_Kp_Km", &dR_Kp_Km_vec);
    mytree->SetBranchAddress("dR_Kp_pi", &dR_Kp_pi_vec);
    mytree->SetBranchAddress("dR_Km_pi", &dR_Km_pi_vec);
    mytree->SetBranchAddress("dR_Kp_Ds", &dR_Kp_Ds_vec);
    mytree->SetBranchAddress("dR_Km_Ds", &dR_Km_Ds_vec);
    mytree->SetBranchAddress("dR_pi_Ds", &dR_pi_Ds_vec);

    mytree->SetBranchAddress("Kp_ipsig", &Kp_ipsig_vec);
    mytree->SetBranchAddress("Kp_ipxysig", &Kp_ipxysig_vec);
    mytree->SetBranchAddress("Km_ipsig", &Km_ipsig_vec);
    mytree->SetBranchAddress("Km_ipxysig", &Km_ipxysig_vec);
    mytree->SetBranchAddress("pi_ipsig", &pi_ipsig_vec);
    mytree->SetBranchAddress("pi_ipxysig", &pi_ipxysig_vec);
    mytree->SetBranchAddress("Ds_ipsig", &Ds_ipsig_vec);
    mytree->SetBranchAddress("Ds_ipxysig", &Ds_ipxysig_vec);

    mytree->SetBranchAddress("Ds_FDsig", &Ds_FDsig_vec);
    mytree->SetBranchAddress("Ds_FDxysig", &Ds_FDxysig_vec);
    mytree->SetBranchAddress("Ds_IsoR04_RelIso", &Ds_IsoR04_RelIso_vec);

    Int_t event_id;
    Float_t Kp_eta;
    Float_t Kp_phi;
    Float_t Kp_pt;
    Float_t Kp_p;
    
    Float_t Km_eta;
    Float_t Km_phi;
    Float_t Km_pt;
    Float_t Km_p;

    Float_t pi_eta;
    Float_t pi_phi;
    Float_t pi_pt;
    Float_t pi_p;

    Float_t KK_eta;
    Float_t KK_phi;
    Float_t KK_pt;
    Float_t KK_p;

    Float_t Kppi_eta;
    Float_t Kppi_phi;
    Float_t Kppi_pt;
    Float_t Kppi_p;

    Float_t Kmpi_eta;
    Float_t Kmpi_phi;
    Float_t Kmpi_pt;
    Float_t Kmpi_p;

    Float_t Ds_eta;
    Float_t Ds_phi;
    Float_t Ds_pt;
    Float_t Ds_p;

    Float_t Ds_prob;
    Float_t dira;

    Float_t dR_Kp_Km;
    Float_t dR_Kp_pi;
    Float_t dR_Km_pi;
    Float_t dR_Kp_Ds;
    Float_t dR_Km_Ds;
    Float_t dR_pi_Ds;

    Float_t Kp_ipsig;
    Float_t Kp_ipxysig;
    Float_t Km_ipsig;
    Float_t Km_ipxysig;
    Float_t pi_ipsig;
    Float_t pi_ipxysig;
    Float_t Ds_ipsig;
    Float_t Ds_ipxysig;

    Float_t Ds_FDsig;
    Float_t Ds_FDxysig;
    Float_t Ds_IsoR04_RelIso;

    Bool_t isSig;
    Bool_t isBkg;

    TFile * MVAfile = new TFile("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/RecoAll_DsstW_KKpimu/MVAflatten.root", "RECREATE");
    TTree * MVAtree = new TTree("Events", "MVAtree");
    MVAtree->Branch("event_id", &event_id);
    MVAtree->Branch("Kp_eta", &Kp_eta);
    MVAtree->Branch("Kp_phi", &Kp_phi);
    MVAtree->Branch("Kp_pt", &Kp_pt);
    MVAtree->Branch("Kp_p", &Kp_p);
    
    MVAtree->Branch("Km_eta", &Km_eta);
    MVAtree->Branch("Km_phi", &Km_phi);
    MVAtree->Branch("Km_pt", &Km_pt);
    MVAtree->Branch("Km_p", &Km_p);

    MVAtree->Branch("pi_eta", &pi_eta);
    MVAtree->Branch("pi_phi", &pi_phi);
    MVAtree->Branch("pi_pt", &pi_pt);
    MVAtree->Branch("pi_p", &pi_p);

    MVAtree->Branch("KK_eta", &KK_eta);
    MVAtree->Branch("KK_phi", &KK_phi);
    MVAtree->Branch("KK_pt", &KK_pt);
    MVAtree->Branch("KK_p", &KK_p);

    MVAtree->Branch("Kppi_eta", &Kppi_eta);
    MVAtree->Branch("Kppi_phi", &Kppi_phi);
    MVAtree->Branch("Kppi_pt", &Kppi_pt);
    MVAtree->Branch("Kppi_p", &Kppi_p);

    MVAtree->Branch("Kmpi_eta", &Kmpi_eta);
    MVAtree->Branch("Kmpi_phi", &Kmpi_phi);
    MVAtree->Branch("Kmpi_pt", &Kmpi_pt);
    MVAtree->Branch("Kmpi_p", &Kmpi_p);

    MVAtree->Branch("Ds_eta", &Ds_eta);
    MVAtree->Branch("Ds_phi", &Ds_phi);
    MVAtree->Branch("Ds_pt", &Ds_pt);
    MVAtree->Branch("Ds_p", &Ds_p);

    MVAtree->Branch("Ds_prob", &Ds_prob);
    MVAtree->Branch("dira", &dira);

    MVAtree->Branch("dR_Kp_Km", &dR_Kp_Km);
    MVAtree->Branch("dR_Kp_pi", &dR_Kp_pi);
    MVAtree->Branch("dR_Km_pi", &dR_Km_pi);
    MVAtree->Branch("dR_Kp_Ds", &dR_Kp_Ds);
    MVAtree->Branch("dR_Km_Ds", &dR_Km_Ds);
    MVAtree->Branch("dR_pi_Ds", &dR_pi_Ds);

    MVAtree->Branch("Kp_ipsig", &Kp_ipsig);
    MVAtree->Branch("Kp_ipxysig", &Kp_ipxysig);
    MVAtree->Branch("Km_ipsig", &Km_ipsig);
    MVAtree->Branch("Km_ipxysig", &Km_ipxysig);
    MVAtree->Branch("pi_ipsig", &pi_ipsig);
    MVAtree->Branch("pi_ipxysig", &pi_ipxysig);
    MVAtree->Branch("Ds_ipsig", &Ds_ipsig);
    MVAtree->Branch("Ds_ipxysig", &Ds_ipxysig);

    MVAtree->Branch("Ds_FDsig", &Ds_FDsig);
    MVAtree->Branch("Ds_FDxysig", &Ds_FDxysig);
    MVAtree->Branch("Ds_IsoR04_RelIso", &Ds_IsoR04_RelIso);

    MVAtree->Branch("isSig", &isSig);
    MVAtree->Branch("isBkg", &isBkg);

    for (Int_t i = 0; i < nEvents; i++)
    {
        mytree->GetEntry(i);

        event_id = i;

        Int_t vecsize = Ds_match_vec->size();
        if (vecsize < 1)
            continue;

        for (Int_t j = 0; j < vecsize; j++)
        {
            Kp_eta = Kp_eta_vec->at(j);
            Kp_phi = Kp_phi_vec->at(j);
            Kp_pt = Kp_pt_vec->at(j);
            Kp_p = Kp_p_vec->at(j);
            
            Km_eta = Km_eta_vec->at(j);
            Km_phi = Km_phi_vec->at(j);
            Km_pt = Km_pt_vec->at(j);
            Km_p = Km_p_vec->at(j);

            pi_eta = pi_eta_vec->at(j);
            pi_phi = pi_phi_vec->at(j);
            pi_pt = pi_pt_vec->at(j);
            pi_p = pi_p_vec->at(j);

            KK_eta = KK_eta_vec->at(j);
            KK_phi = KK_phi_vec->at(j);
            KK_pt = KK_pt_vec->at(j);
            KK_p = KK_p_vec->at(j);

            Kppi_eta = Kppi_eta_vec->at(j);
            Kppi_phi = Kppi_phi_vec->at(j);
            Kppi_pt = Kppi_pt_vec->at(j);
            Kppi_p = Kppi_p_vec->at(j);

            Kmpi_eta = Kmpi_eta_vec->at(j);
            Kmpi_phi = Kmpi_phi_vec->at(j);
            Kmpi_pt = Kmpi_pt_vec->at(j);
            Kmpi_p = Kmpi_p_vec->at(j);

            Ds_eta = Ds_eta_vec->at(j);
            Ds_phi = Ds_phi_vec->at(j);
            Ds_pt = Ds_pt_vec->at(j);
            Ds_p = Ds_p_vec->at(j);

            Ds_prob = Ds_prob_vec->at(j);
            dira = dira_vec->at(j);

            dR_Kp_Km = dR_Kp_Km_vec->at(j);
            dR_Kp_pi = dR_Kp_pi_vec->at(j);
            dR_Km_pi = dR_Km_pi_vec->at(j);
            dR_Kp_Ds = dR_Kp_Ds_vec->at(j);
            dR_Km_Ds = dR_Km_Ds_vec->at(j);
            dR_pi_Ds = dR_pi_Ds_vec->at(j);

            Kp_ipsig = Kp_ipsig_vec->at(j);
            Kp_ipxysig = Kp_ipxysig_vec->at(j);
            Km_ipsig = Km_ipsig_vec->at(j);
            Km_ipxysig = Km_ipxysig_vec->at(j);
            pi_ipsig = pi_ipsig_vec->at(j);
            pi_ipxysig = pi_ipxysig_vec->at(j);
            Ds_ipsig = Ds_ipsig_vec->at(j);
            Ds_ipxysig = Ds_ipxysig_vec->at(j);

            Ds_FDsig = Ds_FDsig_vec->at(j);
            Ds_FDxysig = Ds_FDxysig_vec->at(j);
            Ds_IsoR04_RelIso = Ds_IsoR04_RelIso_vec->at(j);

            isSig = false;
            isBkg = false;
            isSig = Ds_match_vec->at(j);
            isBkg = (!Kp_match_vec->at(j)) && (!Km_match_vec->at(j)) && (!pi_match_vec->at(j));

            if(isSig != isBkg)
                MVAtree->Fill();
        }
    }

    MVAtree->Write();
	MVAfile->Close();

    return 0;
}