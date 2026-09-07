#ifndef H_PHIPIMU_GENTAG_TREE
#define H_PHIPIMU_GENTAG_TREE

#include <TTree.h>

#define null -777

class H_phipimu_GenTag_Tree
{
public:
	H_phipimu_GenTag_Tree(TTree *tree_);
	TTree *tree;

	void CreateBranches();
	void Init();
	void InitGamma();

	ULong64_t event_id;
	UInt_t run_id;
	UInt_t lumi_id;
	Int_t cand_id;

	// Event-level gen counts (last-copy only), for computing reconstruction efficiency:
	// how many true charm hadrons the event contained vs how many ended up as candidates.
	Int_t n_gen_Ds;        // |pdgId| == 431
	Int_t n_gen_D;         // |pdgId| == 411 or 421
	Int_t n_gen_Ds_phipi;  // Ds -> phi pi (phi present among the Ds's direct daughters)

	// trigger
	Bool_t trig_HLT_IsoMu24;
	Bool_t trig_HLT_Mu50;

	// PV
	Bool_t primvtx_isValid;
	Bool_t primvtx_isFake;
	Float_t primvtx_chi2;
	Int_t primvtx_ndof;
	Float_t primvtx_chi2ndof;
	Float_t primvtx_vx;
	Float_t primvtx_vy;
	Float_t primvtx_vz;
	Float_t primvtx_vxerr;
	Float_t primvtx_vyerr;
	Float_t primvtx_vzerr;

	// MET
	Float_t met_sig;
	Float_t met_phi;
	Float_t met_px;
	Float_t met_py;
	Float_t met_pt;
	Float_t met_sumEt;
	Float_t W_tmass;

	// jets
	Int_t num_jets;
	Int_t num_jets_pt20;
	Int_t num_jets_pt30;
	Int_t num_jets_pt40;
	Int_t num_jets_pt50;
	Int_t num_jets_pt60;
	Int_t num_jets_pt70;
	Int_t num_jets_pt80;
	Int_t num_jets_pt90;
	Int_t num_jets_pt100;
	Float_t jet_ht;         // scalar sum pT of tightJetId jets (0 if none)
	Float_t jet_ht_pt20;    // same, jets with pT>20 only
	Float_t jet_leadpt;     // highest-pT tightJetId jet (null if none)
	Float_t jet_maxbscore;  // max b-score over central (|eta|<2.5) tightJetId jets (null if none)
	// false = at least one jet landed in the jet veto map region. Before 2026-08-17 such
	// events were dropped by JvmAppliedEventFilter and never reached this tree, so this
	// branch is absent from every tuple produced earlier; the cut is now applied offline.
	Bool_t jet_vetomap_pass;

	// Per-event generator weight (1.0 on data). Double, not Float: the amcatnloFXFX
	// samples carry |w| ~ 2.6e5, which exceeds Float's ~7 significant digits. Only the
	// sign actually matters for those (weights are two-valued), but keeping the value
	// exact costs nothing and stays correct if a future sample has continuous weights.
	// Normalise with Info tree's sumGenWeights_, never with a raw event count.
	Double_t gen_weight;

	// lepton counts. ss/os = same/opposite charge as the pi (pi == -mu.charge() by
	// construction). os_* includes the best muon (it is opposite the pi).
	Int_t num_os_mu_loose;
	Int_t num_ss_mu_loose;
	Int_t num_os_mu_medium;
	Int_t num_ss_mu_medium;
	Int_t num_os_mu_tight;
	Int_t num_ss_mu_tight;
	Int_t num_os_e_cutloose;
	Int_t num_ss_e_cutloose;
	Int_t num_os_e_cutmedium;
	Int_t num_ss_e_cutmedium;
	Int_t num_os_e_cuttight;
	Int_t num_ss_e_cuttight;
	Int_t num_os_e_mvawp80;
	Int_t num_ss_e_mvawp80;
	Int_t num_os_e_mvawp90;
	Int_t num_ss_e_mvawp90;

	// muon kinematics
	Int_t mu_charge;
	Float_t mu_eta;
	Float_t mu_phi;
	Float_t mu_px;
	Float_t mu_py;
	Float_t mu_pz;
	Float_t mu_pt;
	Float_t mu_p;
	Bool_t mu_isHighPt;
	Bool_t mu_isLoose;
	Bool_t mu_isMedium;
	Bool_t mu_isSoft;
	Bool_t mu_isTight;
	Bool_t mu_isPF;
	Bool_t mu_isTracker;
	Bool_t mu_isGlobal;

	Float_t mu_ip;
	Float_t mu_iperr;
	Float_t mu_ipsig;
	Float_t mu_ipxy;
	Float_t mu_ipxyerr;
	Float_t mu_ipxysig;

	// muon isolation
	Float_t mu_IsoR03_sumChargedHadronPt;
	Float_t mu_IsoR03_sumNeutralHadronEt;
	Float_t mu_IsoR03_sumPhotonEt;
	Float_t mu_IsoR03_sumPUPt;
	Float_t mu_IsoR03_PFIso;
	Float_t mu_IsoR03_RelIso;

	Float_t mu_IsoR04_sumChargedHadronPt;
	Float_t mu_IsoR04_sumNeutralHadronEt;
	Float_t mu_IsoR04_sumPhotonEt;
	Float_t mu_IsoR04_sumPUPt;
	Float_t mu_IsoR04_PFIso;
	Float_t mu_IsoR04_RelIso;

	// refitted K+
	Int_t Kp_charge;
	Float_t Kp_eta;
	Float_t Kp_phi;
	Float_t Kp_px;
	Float_t Kp_py;
	Float_t Kp_pz;
	Float_t Kp_pt;
	Float_t Kp_p;
	Float_t Kp_mass;
	Int_t Kp_PVassoc;

	// refitted K-
	Int_t Km_charge;
	Float_t Km_eta;
	Float_t Km_phi;
	Float_t Km_px;
	Float_t Km_py;
	Float_t Km_pz;
	Float_t Km_pt;
	Float_t Km_p;
	Float_t Km_mass;
	Int_t Km_PVassoc;

	// refitted pi
	Int_t pi_charge;
	Float_t pi_eta;
	Float_t pi_phi;
	Float_t pi_px;
	Float_t pi_py;
	Float_t pi_pz;
	Float_t pi_pt;
	Float_t pi_p;
	Float_t pi_mass;
	Int_t pi_PVassoc;

	// phi = K+K-
	Int_t phi_charge;
	Float_t phi_eta;
	Float_t phi_phi;
	Float_t phi_px;
	Float_t phi_py;
	Float_t phi_pz;
	Float_t phi_pt;
	Float_t phi_p;
	Float_t phi_mass;

	// Ds
	Int_t Ds_charge;
	Float_t Ds_eta;
	Float_t Ds_phi;
	Float_t Ds_px;
	Float_t Ds_py;
	Float_t Ds_pz;
	Float_t Ds_pt;
	Float_t Ds_p;
	Float_t Ds_mass;

	// Ds + mu (always filled)
	Int_t Dsmu_charge;
	Float_t Dsmu_eta;
	Float_t Dsmu_phi;
	Float_t Dsmu_px;
	Float_t Dsmu_py;
	Float_t Dsmu_pz;
	Float_t Dsmu_pt;
	Float_t Dsmu_p;
	Float_t Dsmu_mass;
	Float_t Dsmu_Ht;

	// unconstrained gamma summary per Ds candidate (no iso/mass cut, only dM window [50,250] MeV)
	Int_t   n_gamma_window;
	Float_t gamma_uc_dM;
	Float_t gamma_uc_pt;
	Float_t gamma_uc_dR_Ds;

	// gamma (null if has_gamma == false)
	Bool_t has_gamma;
	Int_t gamma_charge;
	Float_t gamma_eta;
	Float_t gamma_phi;
	Float_t gamma_px;
	Float_t gamma_py;
	Float_t gamma_pz;
	Float_t gamma_pt;
	Float_t gamma_p;
	Float_t gamma_mass;
	Int_t gamma_PVassoc;

	// Ds* = Ds + gamma (null if has_gamma == false)
	Int_t Dsst_charge;
	Float_t Dsst_eta;
	Float_t Dsst_phi;
	Float_t Dsst_px;
	Float_t Dsst_py;
	Float_t Dsst_pz;
	Float_t Dsst_pt;
	Float_t Dsst_p;
	Float_t Dsst_mass;

	// Ds* + mu (null if has_gamma == false)
	Int_t Dsstmu_charge;
	Float_t Dsstmu_eta;
	Float_t Dsstmu_phi;
	Float_t Dsstmu_px;
	Float_t Dsstmu_py;
	Float_t Dsstmu_pz;
	Float_t Dsstmu_pt;
	Float_t Dsstmu_p;
	Float_t Dsstmu_mass;
	Float_t Dsstmu_Ht;

	// refitted vertex
	Float_t Ds_vx;
	Float_t Ds_vy;
	Float_t Ds_vz;
	Float_t Ds_vxerr;
	Float_t Ds_vyerr;
	Float_t Ds_vzerr;
	Float_t Ds_chi2;
	Int_t Ds_ndf;
	Float_t Ds_chi2ndf;
	Float_t Ds_prob;

	// angular distances (Ds system)
	Float_t dR_Kp_Km;
	Float_t dR_Kp_pi;
	Float_t dR_Kp_Ds;
	Float_t dR_Km_pi;
	Float_t dR_Km_Ds;
	Float_t dR_pi_Ds;
	Float_t dR_phi_pi;
	Float_t dR_phi_Ds;
	Float_t dR_Kp_mu;
	Float_t dR_Km_mu;
	Float_t dR_pi_mu;
	Float_t dR_phi_mu;
	Float_t dR_Ds_mu;

	// angular distances (gamma/Dsst, null if has_gamma == false)
	Float_t dR_Kp_gamma;
	Float_t dR_Kp_Dsst;
	Float_t dR_Km_gamma;
	Float_t dR_Km_Dsst;
	Float_t dR_pi_gamma;
	Float_t dR_pi_Dsst;
	Float_t dR_phi_gamma;
	Float_t dR_phi_Dsst;
	Float_t dR_Ds_gamma;
	Float_t dR_Ds_Dsst;
	Float_t dR_gamma_Dsst;
	Float_t dR_gamma_mu;
	Float_t dR_Dsst_mu;

	// IP
	Float_t Kp_ip;
	Float_t Kp_iperr;
	Float_t Kp_ipsig;
	Float_t Kp_ipxy;
	Float_t Kp_ipxyerr;
	Float_t Kp_ipxysig;
	Float_t Km_ip;
	Float_t Km_iperr;
	Float_t Km_ipsig;
	Float_t Km_ipxy;
	Float_t Km_ipxyerr;
	Float_t Km_ipxysig;
	Float_t pi_ip;
	Float_t pi_iperr;
	Float_t pi_ipsig;
	Float_t pi_ipxy;
	Float_t pi_ipxyerr;
	Float_t pi_ipxysig;
	Float_t Ds_ip;
	Float_t Ds_iperr;
	Float_t Ds_ipsig;
	Float_t Ds_ipxy;
	Float_t Ds_ipxyerr;
	Float_t Ds_ipxysig;

	// PV-SV
	Float_t dira;
	Float_t Ds_FDxy;
	Float_t Ds_FDxyerr;
	Float_t Ds_FDxysig;
	Float_t Ds_FD;
	Float_t Ds_FDerr;
	Float_t Ds_FDsig;

	// Ds isolation
	Float_t Ds_IsoR03_sumChargedHadronPt;
	Float_t Ds_IsoR03_sumNeutralHadronEt;
	Float_t Ds_IsoR03_sumPhotonEt;
	Float_t Ds_IsoR03_sumPUPt;
	Float_t Ds_IsoR03_PFIso;
	Float_t Ds_IsoR03_RelIso;
	Float_t Ds_IsoR04_sumChargedHadronPt;
	Float_t Ds_IsoR04_sumNeutralHadronEt;
	Float_t Ds_IsoR04_sumPhotonEt;
	Float_t Ds_IsoR04_sumPUPt;
	Float_t Ds_IsoR04_PFIso;
	Float_t Ds_IsoR04_RelIso;

	// jet nearest the Ds (null if no tightJetId jet in event)
	Float_t Ds_nearestjet_dR;
	Float_t Ds_nearestjet_bscore;

	// central tightJetId jets with b-score > loose/medium/tight WP, excluding jets within
	// dR03/dR04 of this Ds (avoids counting the jet the Ds itself is clustered into).
	// Per-candidate: recomputed for each Ds using this row's Ds_tlv.
	Int_t n_b_jets_loose_dR03;
	Int_t n_b_jets_loose_dR04;
	Int_t n_b_jets_medium_dR03;
	Int_t n_b_jets_medium_dR04;
	Int_t n_b_jets_tight_dR03;
	Int_t n_b_jets_tight_dR04;
	Int_t n_b_jets_loose_pt20_dR03;   // same, jets with pT>20 only
	Int_t n_b_jets_loose_pt20_dR04;
	Int_t n_b_jets_medium_pt20_dR03;
	Int_t n_b_jets_medium_pt20_dR04;
	Int_t n_b_jets_tight_pt20_dR03;
	Int_t n_b_jets_tight_pt20_dR04;

	// Ds* isolation (null if has_gamma == false)
	Float_t Dsst_IsoR03_sumChargedHadronPt;
	Float_t Dsst_IsoR03_sumNeutralHadronEt;
	Float_t Dsst_IsoR03_sumPhotonEt;
	Float_t Dsst_IsoR03_sumPUPt;
	Float_t Dsst_IsoR03_PFIso;
	Float_t Dsst_IsoR03_RelIso;
	Float_t Dsst_IsoR04_sumChargedHadronPt;
	Float_t Dsst_IsoR04_sumNeutralHadronEt;
	Float_t Dsst_IsoR04_sumPhotonEt;
	Float_t Dsst_IsoR04_sumPUPt;
	Float_t Dsst_IsoR04_PFIso;
	Float_t Dsst_IsoR04_RelIso;

	//==================================================//
	// Truth tagging of the three tracks (the point of this tuplizer)
	//==================================================//
	// Per-track truth, from packedPFCandidateToGenAssociation. That association is a
	// hit-level track->TrackingParticle match re-keyed into prunedGenParticles, so a null
	// ref does NOT mean "no truth information" -- it means the matched gen particle is not
	// a category prunedGenParticles keeps. After its "drop status == 1", only first-generation
	// daughters of c/b hadrons, phi(333) daughters, K_S/Lambda daughters, leptons, neutrinos
	// and a few others are re-kept, so an ordinary fragmentation pion is always null.
	// That is exactly the wanted behaviour here: the five modes of interest all have their
	// final-state tracks inside the kept set.
	//
	// Known blind spot, not a bug: a mode whose tracks sit >=2 generations below the charm
	// hadron without passing through a phi (e.g. Ds->K*K, K*->K pi -- K* is kept by a plain
	// "keep" with no "+", so its daughters are pruned) comes back null. Those modes do not
	// produce a KKpi peak at M_Ds, but it does mean this tagging is not a complete
	// decomposition of all Ds decays.
	//
	// The two *_key fields are indices into that event's prunedGenParticles. They are
	// COMPARABLE ONLY WITHIN ONE ROW -- they carry no meaning across rows or events. Their
	// whole purpose is to tell "the same Ds" from "two different Ds" offline.
	Int_t  Kp_gen_pdgId;      Int_t  Km_gen_pdgId;      Int_t  pi_gen_pdgId;
	Int_t  Kp_gen_mom_pdgId;  Int_t  Km_gen_mom_pdgId;  Int_t  pi_gen_mom_pdgId;
	Int_t  Kp_gen_mom_key;    Int_t  Km_gen_mom_key;    Int_t  pi_gen_mom_key;
	Int_t  Kp_gen_hf_pdgId;   Int_t  Km_gen_hf_pdgId;   Int_t  pi_gen_hf_pdgId;
	Int_t  Kp_gen_hf_key;     Int_t  Km_gen_hf_key;     Int_t  pi_gen_hf_key;
	Bool_t Kp_gen_hf_isMom;   Bool_t Km_gen_hf_isMom;   Bool_t pi_gen_hf_isMom;

	// Candidate-level summary, all derivable from the per-track fields above.
	Int_t  cand_gen_nmatched;        // tracks with a non-null association (0-3)
	Int_t  cand_gen_n_distinct_hf;   // distinct c/b ancestors among the three tracks (0-3)
	Bool_t cand_gen_KK_same_mom;     // the two K slots share one immediate mother (by key)
	Int_t  cand_gen_KK_mom_pdgId;    // that shared mother's pdgId (333 => a real phi)
	Int_t  cand_gen_common_hf_pdgId; // non-null only when all three share ONE c/b ancestor
	Bool_t cand_gen_all_direct;      // all three tracks' immediate mother IS that ancestor

	// How many of the three tracks are truly a K / truly a pi. The candidate loop assigns the
	// K,K,pi mass hypotheses blind, so these two counters ARE the misID bookkeeping:
	// 2K+1pi is the intended assignment, 1K+2pi is one track taken for a kaon, 0K+3pi is two.
	// Stored rather than only folded into cand_gen_mode so that any misID classification can
	// be redone offline -- the enum below cannot keep up with every (species, nK) combination.
	Int_t  cand_gen_ntrue_K;
	Int_t  cand_gen_ntrue_pi;

	// Convenience enum. The raw fields above are the primary product -- with the keys stored,
	// any regrouping can be redone offline without reproducing the sample.
	//   0 = no c/b ancestor on any track (pure fragmentation combinatorics)
	//   1 = Ds -> phi pi, phi -> KK
	//   2 = Ds -> K K pi, direct three-body
	//   3 = D  -> phi pi, phi -> KK
	//   4 = D  -> K K pi, direct three-body
	//   5 = D  -> K pi pi, one pi taken as a K (the misID reflection)
	//   6 = one common c/b ancestor but the topology is none of the named ones
	//   7 = >=2 distinct c/b ancestors: real charm tracks from different hadrons
	//   8 = partial -- some tracks have a c/b ancestor, the rest are null
	//   9 = Ds -> K pi pi, one pi taken as a K (the Ds-side twin of mode 5)
	//  10 = D or Ds -> pi pi pi, BOTH K slots misID. Measured to be the dominant
	//       previously-unnamed population (23 of 26 mode-6 rows in a 1J test run, mostly
	//       D0 -> pi pi pi), which is why it earns its own value rather than sitting in 6.
	// Species beyond D/Ds (Lambda_c, B hadrons) deliberately stay in mode 6: their species is
	// already in cand_gen_common_hf_pdgId, and naming every one of them does not scale.
	// Modes 7 and 8 are deliberately NOT folded into 0. Their tracks are genuine charm decay
	// products (displaced, isolated, signal-like to the MVA) but do not share an origin, so
	// M(KKpi) does not peak. They behave like combinatorics in a mass fit and like signal in
	// the classifier, which is precisely the region the SR closure problem lives in.
	Int_t cand_gen_mode;
};
#endif
