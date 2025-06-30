#include <TChain.h>

int DsVariable(){

    TChain * mychain = new TChain("PVStudy/Events");
    mychain->Add("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/PVStudy/PVStudy.root");

    TFile * outfile = new TFile("/pnfs/iihe/cms/store/user/kakang/Analysis/Simulation/20250417/2017UL/tuples/PVStudy/AddVar2.root", "RECREATE");
	TTree * outtree = new TTree("Events", "Events");
    /* TTree * outtree = mychain->CloneTree(0); */

    double Gen_H_ORIVX_X; mychain->SetBranchAddress("Gen_H_ORIVX_X", &Gen_H_ORIVX_X);
    double Gen_H_ORIVX_Y; mychain->SetBranchAddress("Gen_H_ORIVX_Y", &Gen_H_ORIVX_Y);
    double Gen_H_ORIVX_Z; mychain->SetBranchAddress("Gen_H_ORIVX_Z", &Gen_H_ORIVX_Z);
    double Gen_Kp_ORIVX_X; mychain->SetBranchAddress("Gen_Kp_ORIVX_X", &Gen_Kp_ORIVX_X);
    double Gen_Kp_ORIVX_Y; mychain->SetBranchAddress("Gen_Kp_ORIVX_Y", &Gen_Kp_ORIVX_Y);
    double Gen_Kp_ORIVX_Z; mychain->SetBranchAddress("Gen_Kp_ORIVX_Z", &Gen_Kp_ORIVX_Z);
    double Gen_Ds_PX; mychain->SetBranchAddress("Gen_Ds_PX", &Gen_Ds_PX);
    double Gen_Ds_PY; mychain->SetBranchAddress("Gen_Ds_PY", &Gen_Ds_PY);
    double Gen_Ds_PZ; mychain->SetBranchAddress("Gen_Ds_PZ", &Gen_Ds_PZ);
    double Gen_phi_PX; mychain->SetBranchAddress("Gen_phi_PX", &Gen_phi_PX);
    double Gen_phi_PY; mychain->SetBranchAddress("Gen_phi_PY", &Gen_phi_PY);
    double Gen_phi_PZ; mychain->SetBranchAddress("Gen_phi_PZ", &Gen_phi_PZ);
    double Gen_Kp_PX; mychain->SetBranchAddress("Gen_Kp_PX", &Gen_Kp_PX);
    double Gen_Kp_PY; mychain->SetBranchAddress("Gen_Kp_PY", &Gen_Kp_PY);
    double Gen_Kp_PZ; mychain->SetBranchAddress("Gen_Kp_PZ", &Gen_Kp_PZ);
    double Gen_Km_PX; mychain->SetBranchAddress("Gen_Km_PX", &Gen_Km_PX);
    double Gen_Km_PY; mychain->SetBranchAddress("Gen_Km_PY", &Gen_Km_PY);
    double Gen_Km_PZ; mychain->SetBranchAddress("Gen_Km_PZ", &Gen_Km_PZ);
    double Gen_pi_PX; mychain->SetBranchAddress("Gen_pi_PX", &Gen_pi_PX);
    double Gen_pi_PY; mychain->SetBranchAddress("Gen_pi_PY", &Gen_pi_PY);
    double Gen_pi_PZ; mychain->SetBranchAddress("Gen_pi_PZ", &Gen_pi_PZ);

    std::vector<bool> * PV_withBS_IsValid = nullptr; mychain->SetBranchAddress("PV_withBS_IsValid", &PV_withBS_IsValid);
    std::vector<bool> * PV_withBS_IsFake = nullptr; mychain->SetBranchAddress("PV_withBS_IsFake", &PV_withBS_IsFake);
    std::vector<double> * PV_withBS_X = nullptr; mychain->SetBranchAddress("PV_withBS_X", &PV_withBS_X);
    std::vector<double> * PV_withBS_Y = nullptr; mychain->SetBranchAddress("PV_withBS_Y", &PV_withBS_Y);
    std::vector<double> * PV_withBS_Z = nullptr; mychain->SetBranchAddress("PV_withBS_Z", &PV_withBS_Z);
    std::vector<double> * PV_withBS_XERR = nullptr; mychain->SetBranchAddress("PV_withBS_XERR", &PV_withBS_XERR);
    std::vector<double> * PV_withBS_YERR = nullptr; mychain->SetBranchAddress("PV_withBS_YERR", &PV_withBS_YERR);
    std::vector<double> * PV_withBS_ZERR = nullptr; mychain->SetBranchAddress("PV_withBS_ZERR", &PV_withBS_ZERR);
    std::vector<double> * match_DsFit_ENDVX_X = nullptr; mychain->SetBranchAddress("match_DsFit_ENDVX_X", &match_DsFit_ENDVX_X);
    std::vector<double> * match_DsFit_ENDVX_Y = nullptr; mychain->SetBranchAddress("match_DsFit_ENDVX_Y", &match_DsFit_ENDVX_Y);
    std::vector<double> * match_DsFit_ENDVX_Z = nullptr; mychain->SetBranchAddress("match_DsFit_ENDVX_Z", &match_DsFit_ENDVX_Z);
    std::vector<double> * match_DsFit_ENDVX_XERR = nullptr; mychain->SetBranchAddress("match_DsFit_ENDVX_XERR", &match_DsFit_ENDVX_XERR);
    std::vector<double> * match_DsFit_ENDVX_YERR = nullptr; mychain->SetBranchAddress("match_DsFit_ENDVX_YERR", &match_DsFit_ENDVX_YERR);
    std::vector<double> * match_DsFit_ENDVX_ZERR = nullptr; mychain->SetBranchAddress("match_DsFit_ENDVX_ZERR", &match_DsFit_ENDVX_ZERR);
    std::vector<double> * match_DsFit_Ds_PX = nullptr; mychain->SetBranchAddress("match_DsFit_Ds_PX", &match_DsFit_Ds_PX);
    std::vector<double> * match_DsFit_Ds_PY = nullptr; mychain->SetBranchAddress("match_DsFit_Ds_PY", &match_DsFit_Ds_PY);
    std::vector<double> * match_DsFit_Ds_PZ = nullptr; mychain->SetBranchAddress("match_DsFit_Ds_PZ", &match_DsFit_Ds_PZ);
    std::vector<double> * match_phiFit_phi_PX = nullptr; mychain->SetBranchAddress("match_phiFit_phi_PX", &match_phiFit_phi_PX);
    std::vector<double> * match_phiFit_phi_PY = nullptr; mychain->SetBranchAddress("match_phiFit_phi_PY", &match_phiFit_phi_PY);
    std::vector<double> * match_phiFit_phi_PZ = nullptr; mychain->SetBranchAddress("match_phiFit_phi_PZ", &match_phiFit_phi_PZ);
    std::vector<double> * match_Kp_PX = nullptr; mychain->SetBranchAddress("match_Kp_PX", &match_Kp_PX);
    std::vector<double> * match_Kp_PY = nullptr; mychain->SetBranchAddress("match_Kp_PY", &match_Kp_PY);
    std::vector<double> * match_Kp_PZ = nullptr; mychain->SetBranchAddress("match_Kp_PZ", &match_Kp_PZ);
    std::vector<double> * match_Km_PX = nullptr; mychain->SetBranchAddress("match_Km_PX", &match_Km_PX);
    std::vector<double> * match_Km_PY = nullptr; mychain->SetBranchAddress("match_Km_PY", &match_Km_PY);
    std::vector<double> * match_Km_PZ = nullptr; mychain->SetBranchAddress("match_Km_PZ", &match_Km_PZ);
    std::vector<double> * match_pi_PX = nullptr; mychain->SetBranchAddress("match_pi_PX", &match_pi_PX);
    std::vector<double> * match_pi_PY = nullptr; mychain->SetBranchAddress("match_pi_PY", &match_pi_PY);
    std::vector<double> * match_pi_PZ = nullptr; mychain->SetBranchAddress("match_pi_PZ", &match_pi_PZ);

    double Gen_Ds_dxy; outtree->Branch("Gen_Ds_dxy", &Gen_Ds_dxy);
    double Gen_Ds_dz; outtree->Branch("Gen_Ds_dz", &Gen_Ds_dz);
    double Gen_Ds_FD; outtree->Branch("Gen_Ds_FD", &Gen_Ds_FD);
    double Gen_Ds_DIRA; outtree->Branch("Gen_Ds_DIRA", &Gen_Ds_DIRA);
    double Gen_phi_IP; outtree->Branch("Gen_phi_IP", &Gen_phi_IP);
    double Gen_Kp_IP; outtree->Branch("Gen_Kp_IP", &Gen_Kp_IP);
    double Gen_Km_IP; outtree->Branch("Gen_Km_IP", &Gen_Km_IP);
    double Gen_pi_IP; outtree->Branch("Gen_pi_IP", &Gen_pi_IP);

    std::vector<double> * match_Ds_dxy = nullptr; outtree->Branch("match_Ds_dxy", &match_Ds_dxy);
    std::vector<double> * match_Ds_dz = nullptr; outtree->Branch("match_Ds_dz", &match_Ds_dz);
    std::vector<double> * match_Ds_FD = nullptr; outtree->Branch("match_Ds_FD", &match_Ds_FD);
    std::vector<double> * match_Ds_DIRA = nullptr; outtree->Branch("match_Ds_DIRA", &match_Ds_DIRA);
    std::vector<double> * match_phi_IP = nullptr; outtree->Branch("match_phi_IP", &match_phi_IP);
    std::vector<double> * match_Kp_IP = nullptr; outtree->Branch("match_Kp_IP", &match_Kp_IP);
    std::vector<double> * match_Km_IP = nullptr; outtree->Branch("match_Km_IP", &match_Km_IP);
    std::vector<double> * match_pi_IP = nullptr; outtree->Branch("match_pi_IP", &match_pi_IP);
    
    std::vector<double> * match_Ds_dxy_ERR = nullptr; outtree->Branch("match_Ds_dxy_ERR", &match_Ds_dxy_ERR);
    std::vector<double> * match_Ds_dz_ERR = nullptr; outtree->Branch("match_Ds_dz_ERR", &match_Ds_dz_ERR);
    std::vector<double> * match_Ds_FD_ERR = nullptr; outtree->Branch("match_Ds_FD_ERR", &match_Ds_FD_ERR);
    std::vector<double> * match_phi_IP_ERR = nullptr; outtree->Branch("match_phi_IP_ERR", &match_phi_IP_ERR);
    std::vector<double> * match_Kp_IP_ERR = nullptr; outtree->Branch("match_Kp_IP_ERR", &match_Kp_IP_ERR);
    std::vector<double> * match_Km_IP_ERR = nullptr; outtree->Branch("match_Km_IP_ERR", &match_Km_IP_ERR);
    std::vector<double> * match_pi_IP_ERR = nullptr; outtree->Branch("match_pi_IP_ERR", &match_pi_IP_ERR);

    std::vector<double> * match_Ds_FD_Chi2 = nullptr; outtree->Branch("match_Ds_FD_Chi2", &match_Ds_FD_Chi2);
    std::vector<double> * match_phi_IP_Chi2 = nullptr; outtree->Branch("match_phi_IP_Chi2", &match_phi_IP_Chi2);
    std::vector<double> * match_Kp_IP_Chi2 = nullptr; outtree->Branch("match_Kp_IP_Chi2", &match_Kp_IP_Chi2);
    std::vector<double> * match_Km_IP_Chi2 = nullptr; outtree->Branch("match_Km_IP_Chi2", &match_Km_IP_Chi2);
    std::vector<double> * match_pi_IP_Chi2 = nullptr; outtree->Branch("match_pi_IP_Chi2", &match_pi_IP_Chi2);

    int nentries = mychain->GetEntries();
    for(int i=0; i<nentries; i++){

        mychain->GetEntry(i);

        TVector3 Gen_PV(Gen_H_ORIVX_X, Gen_H_ORIVX_Y, Gen_H_ORIVX_Z); 
        TVector3 Gen_Ds_ENDVX(Gen_Kp_ORIVX_X, Gen_Kp_ORIVX_Y, Gen_Kp_ORIVX_Z);
        TVector3 Gen_Ds_FDirec = Gen_Ds_ENDVX - Gen_PV;
        TVector3 Gen_Ds_PDirec(Gen_Ds_PX, Gen_Ds_PY, Gen_Ds_PZ);
        TVector3 Gen_phi_PDirec(Gen_phi_PX, Gen_phi_PY, Gen_phi_PZ);
        TVector3 Gen_Kp_PDirec(Gen_Kp_PX, Gen_Kp_PY, Gen_Kp_PZ);
        TVector3 Gen_Km_PDirec(Gen_Km_PX, Gen_Km_PY, Gen_Km_PZ);
        TVector3 Gen_pi_PDirec(Gen_pi_PX, Gen_pi_PY, Gen_pi_PZ);

        Gen_Ds_dxy = Gen_Ds_FDirec.Perp(); 
        Gen_Ds_dz = Gen_Ds_FDirec.z(); 
        Gen_Ds_FD = Gen_Ds_FDirec.Mag();
        Gen_Ds_DIRA = Gen_Ds_FDirec.Dot(Gen_Ds_PDirec) / ( Gen_Ds_FDirec.Mag() * Gen_Ds_PDirec.Mag() );
        Gen_phi_IP = Gen_Ds_FDirec.Perp(Gen_phi_PDirec);
        Gen_Kp_IP = Gen_Ds_FDirec.Perp(Gen_Kp_PDirec);
        Gen_Km_IP = Gen_Ds_FDirec.Perp(Gen_Km_PDirec);
        Gen_pi_IP = Gen_Ds_FDirec.Perp(Gen_pi_PDirec);


        match_Ds_dxy->clear();
        match_Ds_dxy_ERR->clear();
        match_Ds_dz->clear();
        match_Ds_dz_ERR->clear();
        match_Ds_FD->clear();
        match_Ds_FD_ERR->clear();
        match_Ds_DIRA->clear();
        match_phi_IP->clear();
        match_Kp_IP->clear();
        match_Km_IP->clear();
        match_pi_IP->clear();
        if( PV_withBS_IsValid->at(0) && !(PV_withBS_IsFake->at(0)) && (match_Kp_PX->size()>0) ){

            TVector3 PV(PV_withBS_X->at(0), PV_withBS_Y->at(0), PV_withBS_Z->at(0));
            TVector3 PV_ERR(PV_withBS_XERR->at(0), PV_withBS_YERR->at(0), PV_withBS_ZERR->at(0));
            TVector3 Ds_ENDVX(match_DsFit_ENDVX_X->at(0), match_DsFit_ENDVX_Y->at(0), match_DsFit_ENDVX_Z->at(0));
            TVector3 Ds_ENDVX_ERR(match_DsFit_ENDVX_XERR->at(0), match_DsFit_ENDVX_YERR->at(0), match_DsFit_ENDVX_ZERR->at(0));

            TVector3 Ds_FDirec = Ds_ENDVX - PV;
            TVector3 Ds_PDirec(match_DsFit_Ds_PX->at(0), match_DsFit_Ds_PY->at(0), match_DsFit_Ds_PZ->at(0));
            TVector3 phi_PDirec(match_phiFit_phi_PX->at(0), match_phiFit_phi_PY->at(0), match_phiFit_phi_PZ->at(0));
            TVector3 Kp_PDirec(match_Kp_PX->at(0), match_Kp_PY->at(0), match_Kp_PZ->at(0));
            TVector3 Km_PDirec(match_Km_PX->at(0), match_Km_PY->at(0), match_Km_PZ->at(0));
            TVector3 pi_PDirec(match_pi_PX->at(0), match_pi_PY->at(0), match_pi_PZ->at(0));
            
            double Ds_dxy = Ds_FDirec.Perp();
            double Ds_dxy_ERR = std::sqrt( pow(Ds_FDirec.x(),2) * ( pow(PV_ERR.x(),2) + pow(Ds_ENDVX_ERR.x(),2) ) + pow(Ds_FDirec.y(),2) * ( pow(PV_ERR.y(),2) + pow(Ds_ENDVX_ERR.y(),2) ) ) / Ds_dxy;
            double Ds_dz = Ds_FDirec.z();
            double Ds_dz_ERR = std::sqrt( pow(PV_ERR.z(),2) + pow(Ds_ENDVX_ERR.z(),2) );
            double Ds_FD = Ds_FDirec.Mag();
            double Ds_FD_ERR = std::sqrt( pow(Ds_FDirec.x(),2) * ( pow(PV_ERR.x(),2) + pow(Ds_ENDVX_ERR.x(),2) ) + pow(Ds_FDirec.y(),2) * ( pow(PV_ERR.y(),2) + pow(Ds_ENDVX_ERR.y(),2) ) + pow(Ds_FDirec.z(),2) * ( pow(PV_ERR.z(),2) + pow(Ds_ENDVX_ERR.z(),2) ) ) / Ds_FD;
            double Ds_DIRA = Ds_FDirec.Dot(Ds_PDirec) / ( Ds_FDirec.Mag() * Ds_PDirec.Mag() );

            double phi_IP = Ds_FDirec.Perp(phi_PDirec);
            double Kp_IP = Ds_FDirec.Perp(Kp_PDirec);
            double Km_IP = Ds_FDirec.Perp(Km_PDirec);
            double pi_IP = Ds_FDirec.Perp(pi_PDirec);
            
            match_Ds_dxy->push_back(Ds_dxy);
            match_Ds_dxy_ERR->push_back(Ds_dxy_ERR);
            match_Ds_dz->push_back(Ds_dz);
            match_Ds_dz_ERR->push_back(Ds_dz_ERR);
            match_Ds_FD->push_back(Ds_FD);
            match_Ds_FD_ERR->push_back(Ds_FD_ERR);
            match_Ds_DIRA->push_back(Ds_DIRA);
            match_phi_IP->push_back(phi_IP);
            match_Kp_IP->push_back(Kp_IP);
            match_Km_IP->push_back(Km_IP);
            match_pi_IP->push_back(pi_IP);
        }
        outtree->Fill();
    }

    outtree->Write();
	delete outtree;
	outfile->Close();
	delete outfile;

    return 0;
}
