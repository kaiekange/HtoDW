#include <TChain.h>

const Float_t Mass_Ds = 1.96835;
const Float_t Mass_phi = 1.019460;
const Float_t Mass_K = 0.493677;
const Float_t Mass_pi = 0.13957039;
const Float_t Mass_mu = 0.1056583755;

Int_t sig_eff()
{

	TChain *mychain = new TChain("MatchRecoKineFit/Events");
	mychain->Add("/afs/cern.ch/work/k/kakang/H2DW/CMSSW_13_0_13/src/analysis/MatchRecoKineFit/MatchRecoKineFit.root");

	Bool_t trig_HLT_IsoMu24;
	Bool_t trig_HLT_Mu50;

	Float_t orig_Kp_pt;
	Float_t orig_Kp_px;
	Float_t orig_Kp_py;
	Float_t orig_Kp_pz;

	Float_t orig_Km_pt;
	Float_t orig_Km_px;
	Float_t orig_Km_py;
	Float_t orig_Km_pz;

	Float_t orig_pi_pt;
	Float_t orig_pi_px;
	Float_t orig_pi_py;
	Float_t orig_pi_pz;

	Float_t Ds_prob;

	Float_t phi_pt;
	Float_t phi_mass;

	Float_t Ds_pt;
	Float_t Ds_mass;

	Float_t dira;

	Float_t mu_pt;
	Bool_t mu_isTight;

	Int_t num_trig = 0;
	Int_t num_orig_Kp_pt = 0;
	Int_t num_orig_Km_pt = 0;
	Int_t num_orig_pi_pt = 0;
	Int_t num_orig_phi_pt = 0;
	Int_t num_orig_phi_m = 0;
	Int_t num_orig_Ds_pt = 0;
	Int_t num_orig_Ds_m = 0;
	Int_t num_Ds_prob = 0;
	Int_t num_phi_pt = 0;
	Int_t num_phi_mass = 0;
	Int_t num_Ds_pt = 0;
	Int_t num_Ds_mass = 0;
	Int_t num_dira = 0;
	Int_t num_mu_pt = 0;
	Int_t num_mu_isTight = 0;

	mychain->SetBranchAddress("trig_HLT_IsoMu24", &trig_HLT_IsoMu24);
	mychain->SetBranchAddress("trig_HLT_Mu50", &trig_HLT_Mu50);

	mychain->SetBranchAddress("orig_Kp_pt", &orig_Kp_pt);
	mychain->SetBranchAddress("orig_Kp_px", &orig_Kp_px);
	mychain->SetBranchAddress("orig_Kp_py", &orig_Kp_py);
	mychain->SetBranchAddress("orig_Kp_pz", &orig_Kp_pz);

	mychain->SetBranchAddress("orig_Km_pt", &orig_Km_pt);
	mychain->SetBranchAddress("orig_Km_px", &orig_Km_px);
	mychain->SetBranchAddress("orig_Km_py", &orig_Km_py);
	mychain->SetBranchAddress("orig_Km_pz", &orig_Km_pz);

	mychain->SetBranchAddress("orig_pi_pt", &orig_pi_pt);
	mychain->SetBranchAddress("orig_pi_px", &orig_pi_px);
	mychain->SetBranchAddress("orig_pi_py", &orig_pi_py);
	mychain->SetBranchAddress("orig_pi_pz", &orig_pi_pz);

	mychain->SetBranchAddress("Ds_prob", &Ds_prob);

	mychain->SetBranchAddress("phi_pt", &phi_pt);
	mychain->SetBranchAddress("phi_mass", &phi_mass);

	mychain->SetBranchAddress("Ds_pt", &Ds_pt);
	mychain->SetBranchAddress("Ds_mass", &Ds_mass);

	mychain->SetBranchAddress("dira", &dira);

	mychain->SetBranchAddress("mu_pt", &mu_pt);
	mychain->SetBranchAddress("mu_isTight", &mu_isTight);

	Int_t nevents = mychain->GetEntries();

	std::cout << nevents << std::endl;

	for (Int_t i = 0; i < nevents; i++)
	{
		mychain->GetEntry(i);

		TLorentzVector Kp_tlv, Km_tlv, pi_tlv, phi_tlv, Ds_tlv;

		Kp_tlv.SetXYZM(orig_Kp_px, orig_Kp_py, orig_Kp_pz, Mass_K);
		Km_tlv.SetXYZM(orig_Km_px, orig_Km_py, orig_Km_pz, Mass_K);
		pi_tlv.SetXYZM(orig_pi_px, orig_pi_py, orig_pi_pz, Mass_pi);
		phi_tlv = Kp_tlv + Km_tlv;
		Ds_tlv = Kp_tlv + Km_tlv + pi_tlv;

		if (!(trig_HLT_IsoMu24 || trig_HLT_Mu50))
			continue;
		num_trig++;

		if (orig_Kp_pt < 1)
			continue;
		num_orig_Kp_pt++;

		if (orig_Km_pt < 1)
			continue;
		num_orig_Km_pt++;

		if (orig_pi_pt < 0.5)
			continue;
		num_orig_pi_pt++;

		if (phi_tlv.Pt() < 2.5)
			continue;
		num_orig_phi_pt++;

		if (abs(phi_tlv.M() - Mass_phi) > 0.1)
			continue;
		num_orig_phi_m++;

		if (Ds_tlv.Pt() < 5)
			continue;
		num_orig_Ds_pt++;

		if (abs(Ds_tlv.M() - Mass_Ds) > 0.3)
			continue;
		num_orig_Ds_m++;

		if (Ds_prob < 0.01)
			continue;
		num_Ds_prob++;

		if (abs(phi_mass - Mass_phi) > 0.03)
			continue;
		num_phi_mass++;

		if (phi_pt < 3)
			continue;
		num_phi_pt++;

		if (abs(Ds_mass - Mass_Ds) > 0.15)
			continue;
		num_Ds_mass++;

		if (Ds_pt < 7)
			continue;
		num_Ds_pt++;

		if (dira < 0.9)
			continue;
		num_dira++;

		if (mu_pt < 24)
			continue;
		num_mu_pt++;

		if (!mu_isTight)
			continue;
		num_mu_isTight++;
	}



	std::cout << "num_trig " << num_trig << ", eff: " << Float_t(num_trig)/nevents*100 << ", eff/trig: " << Float_t(num_trig)/nevents*100 << std::endl;
	std::cout << "num_orig_Kp_pt " << num_orig_Kp_pt << ", eff: " << Float_t(num_orig_Kp_pt)/nevents*100 << ", eff/trig: " << Float_t(num_orig_Kp_pt)/nevents*100 << std::endl;
	std::cout << "num_orig_Km_pt " << num_orig_Km_pt << ", eff: " << Float_t(num_orig_Km_pt)/nevents*100 << ", eff/trig: " << Float_t(num_orig_Km_pt)/nevents*100 << std::endl;
	std::cout << "num_orig_pi_pt " << num_orig_pi_pt << ", eff: " << Float_t(num_orig_pi_pt)/nevents*100 << ", eff/trig: " << Float_t(num_orig_pi_pt)/nevents*100 << std::endl;
	std::cout << "num_orig_phi_pt " << num_orig_phi_pt << ", eff: " << Float_t(num_orig_phi_pt)/nevents*100 << ", eff/trig: " << Float_t(num_orig_phi_pt)/nevents*100 << std::endl;
	std::cout << "num_orig_phi_m " << num_orig_phi_m << ", eff: " << Float_t(num_orig_phi_m)/nevents*100 << ", eff/trig: " << Float_t(num_orig_phi_m)/nevents*100 << std::endl;
	std::cout << "num_orig_Ds_pt " << num_orig_Ds_pt << ", eff: " << Float_t(num_orig_Ds_pt)/nevents*100 << ", eff/trig: " << Float_t(num_orig_Ds_pt)/nevents*100 << std::endl;
	std::cout << "num_orig_Ds_m " << num_orig_Ds_m << ", eff: " << Float_t(num_orig_Ds_m)/nevents*100 << ", eff/trig: " << Float_t(num_orig_Ds_m)/nevents*100 << std::endl;
	std::cout << "num_Ds_prob " << num_Ds_prob << ", eff: " << Float_t(num_Ds_prob)/nevents*100 << ", eff/trig: " << Float_t(num_Ds_prob)/nevents*100 << std::endl;
	std::cout << "num_phi_pt " << num_phi_pt << ", eff: " << Float_t(num_phi_pt)/nevents*100 << ", eff/trig: " << Float_t(num_phi_pt)/nevents*100 << std::endl;
	std::cout << "num_phi_mass " << num_phi_mass << ", eff: " << Float_t(num_phi_mass)/nevents*100 << ", eff/trig: " << Float_t(num_phi_mass)/nevents*100 << std::endl;
	std::cout << "num_Ds_pt " << num_Ds_pt << ", eff: " << Float_t(num_Ds_pt)/nevents*100 << ", eff/trig: " << Float_t(num_Ds_pt)/nevents*100 << std::endl;
	std::cout << "num_Ds_mass " << num_Ds_mass << ", eff: " << Float_t(num_Ds_mass)/nevents*100 << ", eff/trig: " << Float_t(num_Ds_mass)/nevents*100 << std::endl;
	std::cout << "num_dira " << num_dira << ", eff: " << Float_t(num_dira)/nevents*100 << ", eff/trig: " << Float_t(num_dira)/nevents*100 << std::endl;
	std::cout << "num_mu_pt " << num_mu_pt << ", eff: " << Float_t(num_mu_pt)/nevents*100 << ", eff/trig: " << Float_t(num_mu_pt)/nevents*100 << std::endl;
	std::cout << "num_mu_isTight " << num_mu_isTight << ", eff: " << Float_t(num_mu_isTight)/nevents*100 << ", eff/trig: " << Float_t(num_mu_isTight)/nevents*100 << std::endl;

	std::cout << "\n\n\n" << std::endl;

	std::cout << "eff trig " << Float_t(num_trig)/nevents*100 << std::endl;
	std::cout << "eff orig_pi_pt " << Float_t(num_orig_pi_pt)/num_trig*100 << std::endl;
	std::cout << "eff orig_phi_m " << Float_t(num_orig_phi_m)/num_orig_pi_pt*100 << std::endl;
	std::cout << "eff orig_Ds_m " << Float_t(num_orig_Ds_m)/num_orig_phi_m*100 << std::endl;
	std::cout << "eff Ds_prob " << Float_t(num_Ds_prob)/num_orig_Ds_m*100 << std::endl;
	std::cout << "eff phi_mass " << Float_t(num_phi_mass)/num_Ds_prob*100 << std::endl;
	std::cout << "eff Ds_mass " << Float_t(num_Ds_mass)/num_phi_mass*100 << std::endl;
	std::cout << "eff dira " << Float_t(num_dira)/num_Ds_mass*100 << std::endl;
	std::cout << "eff mu_isTight " << Float_t(num_mu_isTight)/num_dira*100 << std::endl;

	return 0;
}