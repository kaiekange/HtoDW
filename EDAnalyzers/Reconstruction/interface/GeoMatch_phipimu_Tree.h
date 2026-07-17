#ifndef GEOMATCH_PHIPIMU_TREE
#define GEOMATCH_PHIPIMU_TREE

#include <TTree.h>

#define null -777

class GeoMatch_phipimu_Tree
{
public:
	GeoMatch_phipimu_Tree(TTree *tree_);
	TTree *tree;

	void CreateBranches();
	void Init();

	ULong64_t event_id;
	UInt_t run_id;
	UInt_t lumi_id;

	// Gen PV / SV
	Float_t Gen_PV_vx;
	Float_t Gen_PV_vy;
	Float_t Gen_PV_vz;
	Float_t Gen_SV_vx;
	Float_t Gen_SV_vy;
	Float_t Gen_SV_vz;

	// Gen H
	Int_t Gen_H_charge;
	Float_t Gen_H_eta;
	Float_t Gen_H_phi;
	Float_t Gen_H_px;
	Float_t Gen_H_py;
	Float_t Gen_H_pz;
	Float_t Gen_H_pt;
	Float_t Gen_H_p;
	Float_t Gen_H_mass;

	// Gen Ds* (DsstW channel only, else null)
	Int_t Gen_Dsst_charge;
	Float_t Gen_Dsst_eta;
	Float_t Gen_Dsst_phi;
	Float_t Gen_Dsst_px;
	Float_t Gen_Dsst_py;
	Float_t Gen_Dsst_pz;
	Float_t Gen_Dsst_pt;
	Float_t Gen_Dsst_p;
	Float_t Gen_Dsst_mass;

	// Gen Ds
	Int_t Gen_Ds_charge;
	Float_t Gen_Ds_eta;
	Float_t Gen_Ds_phi;
	Float_t Gen_Ds_px;
	Float_t Gen_Ds_py;
	Float_t Gen_Ds_pz;
	Float_t Gen_Ds_pt;
	Float_t Gen_Ds_p;
	Float_t Gen_Ds_mass;

	// Gen gamma (DsstW channel only, else null)
	Int_t Gen_gamma_charge;
	Float_t Gen_gamma_eta;
	Float_t Gen_gamma_phi;
	Float_t Gen_gamma_px;
	Float_t Gen_gamma_py;
	Float_t Gen_gamma_pz;
	Float_t Gen_gamma_pt;
	Float_t Gen_gamma_p;
	Float_t Gen_gamma_mass;

	// Gen W
	Int_t Gen_W_charge;
	Float_t Gen_W_eta;
	Float_t Gen_W_phi;
	Float_t Gen_W_px;
	Float_t Gen_W_py;
	Float_t Gen_W_pz;
	Float_t Gen_W_pt;
	Float_t Gen_W_p;
	Float_t Gen_W_mass;

	// Gen K+
	Int_t Gen_Kp_charge;
	Float_t Gen_Kp_eta;
	Float_t Gen_Kp_phi;
	Float_t Gen_Kp_px;
	Float_t Gen_Kp_py;
	Float_t Gen_Kp_pz;
	Float_t Gen_Kp_pt;
	Float_t Gen_Kp_p;
	Float_t Gen_Kp_mass;

	// Gen K-
	Int_t Gen_Km_charge;
	Float_t Gen_Km_eta;
	Float_t Gen_Km_phi;
	Float_t Gen_Km_px;
	Float_t Gen_Km_py;
	Float_t Gen_Km_pz;
	Float_t Gen_Km_pt;
	Float_t Gen_Km_p;
	Float_t Gen_Km_mass;

	// Gen pi
	Int_t Gen_pi_charge;
	Float_t Gen_pi_eta;
	Float_t Gen_pi_phi;
	Float_t Gen_pi_px;
	Float_t Gen_pi_py;
	Float_t Gen_pi_pz;
	Float_t Gen_pi_pt;
	Float_t Gen_pi_p;
	Float_t Gen_pi_mass;

	// Gen mu
	Int_t Gen_mu_charge;
	Float_t Gen_mu_eta;
	Float_t Gen_mu_phi;
	Float_t Gen_mu_px;
	Float_t Gen_mu_py;
	Float_t Gen_mu_pz;
	Float_t Gen_mu_pt;
	Float_t Gen_mu_p;
	Float_t Gen_mu_mass;

	// Gen Ds+mu
	Int_t Gen_Dsmu_charge;
	Float_t Gen_Dsmu_eta;
	Float_t Gen_Dsmu_phi;
	Float_t Gen_Dsmu_px;
	Float_t Gen_Dsmu_py;
	Float_t Gen_Dsmu_pz;
	Float_t Gen_Dsmu_pt;
	Float_t Gen_Dsmu_p;
	Float_t Gen_Dsmu_mass;

	// Gen Ds*+mu (DsstW channel only, else null)
	Int_t Gen_Dsstmu_charge;
	Float_t Gen_Dsstmu_eta;
	Float_t Gen_Dsstmu_phi;
	Float_t Gen_Dsstmu_px;
	Float_t Gen_Dsstmu_py;
	Float_t Gen_Dsstmu_pz;
	Float_t Gen_Dsstmu_pt;
	Float_t Gen_Dsstmu_p;
	Float_t Gen_Dsstmu_mass;

	// Gen dR (Ds daughters)
	Float_t Gen_dR_Kp_Km;
	Float_t Gen_dR_Kp_pi;
	Float_t Gen_dR_Kp_Ds;
	Float_t Gen_dR_Km_pi;
	Float_t Gen_dR_Km_Ds;
	Float_t Gen_dR_pi_Ds;

	// Gen dR (gamma/Dsst, DsstW channel only, else null)
	Float_t Gen_dR_Kp_gamma;
	Float_t Gen_dR_Kp_Dsst;
	Float_t Gen_dR_Km_gamma;
	Float_t Gen_dR_Km_Dsst;
	Float_t Gen_dR_pi_gamma;
	Float_t Gen_dR_pi_Dsst;
	Float_t Gen_dR_Ds_gamma;
	Float_t Gen_dR_Ds_Dsst;
	Float_t Gen_dR_gamma_Dsst;

	// Gen dR (mu)
	Float_t Gen_dR_Kp_mu;
	Float_t Gen_dR_Km_mu;
	Float_t Gen_dR_pi_mu;
	Float_t Gen_dR_Ds_mu;
	Float_t Gen_dR_gamma_mu;   // DsstW channel only
	Float_t Gen_dR_Dsst_mu;    // DsstW channel only

	// Gen MET
	Float_t Gen_met_phi;
	Float_t Gen_met_px;
	Float_t Gen_met_py;
	Float_t Gen_met_pt;
	Float_t Gen_met_sumEt;

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

	// muon count
	Int_t num_rs_mu;
	Int_t num_ws_mu;

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

	Bool_t mu_match;

	// ------------- association (packedPFCandidateToGenAssociation) -------------
	// Which gen particle a candidate is remains decided by the map plus the
	// pdgId + TwoBodyDecay ancestry chain, exactly as in H_phipimu. What is
	// dropped here is the dR < 0.05 / |dpt| < 0.3 gate: the map, not geometry,
	// identifies the particle, so the gate is post-hoc QA and is applied
	// offline instead. That keeps the working point scannable without a rerun.
	Bool_t Kp_pf_match;
	Bool_t Km_pf_match;
	Bool_t pi_pf_match;
	Bool_t gamma_pf_match;   // DsstW channel only; geometric already, left as-is

	// index into packedPFCandidates (null if no match)
	Int_t Kp_pf_idx;
	Int_t Km_pf_idx;
	Int_t pi_pf_idx;

	// genRef.key() into prunedGenParticles. Compare offline against
	// Gen_<dau>_idx: the Association is built on prunedGenParticles (see
	// PackedCandidateGenAssociationProducer.cc), the same collection the gen
	// filter indexes, so the two must agree. The cross-tag below rests on that.
	Int_t Kp_pf_key;
	Int_t Km_pf_key;
	Int_t pi_pf_key;

	// ------------- geometric (gen-driven, no association map) -----------------
	// For each gen daughter, keep every packedPFCandidate with |pdgId| == 211 and
	// that daughter's charge which passes BOTH dR < 0.05 and |dpt| < 0.3, then
	// take the smallest dR among the survivors. PF calls every charged hadron
	// +-211, so K/pi identity is unavailable; charge is.
	//
	// Both cuts sit INSIDE the selection, matching the gamma matching in
	// H_phipimu. This is not the same as taking the nearest PF and cutting
	// afterwards: a soft track that happens to be nearest must not block a good
	// candidate slightly further out. Consequently these branches are NOT
	// scannable offline -- only the chosen candidate is stored, so a looser
	// working point needs a rerun.
	//
	// The direction is not reversible. This matching runs gen-driven (daughter ->
	// nearest PF) against prunedGenParticles, the same collection the association
	// keys into, which is what makes the cross-tag an exact index compare.
	// Running it PF-driven (PF -> nearest gen) would have to read
	// packedGenParticles instead: prunedGenParticles applies "drop status == 1"
	// and then re-keeps selectively, so it is not a complete final-state list. A
	// PF-driven scan over it would miss the genuinely-nearest (pruned-away) gen
	// particle, misassign the PF to one of our daughters, and systematically
	// overstate the geometric efficiency.
	Bool_t Kp_geo_match;
	Bool_t Km_geo_match;
	Bool_t pi_geo_match;

	// index into packedPFCandidates (null if no match). K+ and pi carry the same
	// charge and can be won by the SAME PF; that collision is recorded as equal
	// indices and resolved offline, not here.
	Int_t Kp_geo_idx;
	Int_t Km_geo_idx;
	Int_t pi_geo_idx;

	// gen-match quality, both methods (null if no match). dR to the gen daughter and
	// the momentum residual (pt_reco - pt_gen) / pt_gen.
	// association is UNGATED, so its tails are visible and its working point is
	// scannable offline. geometric is bounded by its own matching criteria
	// (dR < 0.05, |dpt| < 0.3) because those cuts decide which candidate wins.
	Float_t Kp_pf_dR;
	Float_t Kp_pf_dpt;
	Float_t Km_pf_dR;
	Float_t Km_pf_dpt;
	Float_t pi_pf_dR;
	Float_t pi_pf_dpt;

	Float_t Kp_geo_dR;
	Float_t Kp_geo_dpt;
	Float_t Km_geo_dR;
	Float_t Km_geo_dpt;
	Float_t pi_geo_dR;
	Float_t pi_geo_dpt;

	// match surviving the charged-hadron pool cuts of H_phipimu
	// (hasTrackDetails, trackHighPurity, |pdgId|=211, pt>2 for K and pt>1 for pi)
	Bool_t Kp_pf_pool;
	Bool_t Km_pf_pool;
	Bool_t pi_pf_pool;
	Bool_t Kp_geo_pool;
	Bool_t Km_geo_pool;
	Bool_t pi_geo_pool;

	// kinematics of the geometric pick, so the offline macro can rebuild
	// M(KKpi) under K/K/pi mass hypotheses. That mass is the
	// association-independent closure test of whether the assoc-null bucket
	// below is genuine (a Ds peak) or fake (a broad shape).
	Int_t   Kp_geo_charge;
	Float_t Kp_geo_px;
	Float_t Kp_geo_py;
	Float_t Kp_geo_pz;
	Float_t Kp_geo_pt;
	Float_t Kp_geo_eta;
	Int_t   Km_geo_charge;
	Float_t Km_geo_px;
	Float_t Km_geo_py;
	Float_t Km_geo_pz;
	Float_t Km_geo_pt;
	Float_t Km_geo_eta;
	Int_t   pi_geo_charge;
	Float_t pi_geo_px;
	Float_t pi_geo_py;
	Float_t pi_geo_pz;
	Float_t pi_geo_pt;
	Float_t pi_geo_eta;

	// ---------------- the cross-tag ----------------
	// For the PF geometric picked, look up ITS OWN association and store
	// the raw key + pdgId. Classification happens offline, never here:
	//   key == Gen_<dau>_idx -> confirmed correct
	//   key == null (isNull) -> the geo pick has no association at all
	//   otherwise            -> the geo pick is truthfully a different particle
	//
	// The key==null bucket is where geometry could genuinely beat the map: a real
	// track the association loses to the associator's hit-sharing cut
	// (quickTrackAssociatorByHits Cut_RecoToSim=0.75 -- share >75% of hits with
	// the TrackingParticle or get no association at all). Since geometric now applies
	// |dpt| < 0.3 during selection, anything in this bucket already has sane
	// momentum, so it is a real recovery candidate rather than soft PU/UE junk.
	// (prunedTrackingParticles pruning is NOT the mechanism: that selector is
	// "keep++ (400<|pdgId|<600)" + "drop status!=1" + "drop charge==0", which
	// keeps the Ds/Ds* descendants -- our K/K/pi -- by design.)
	Int_t Kp_geo_assoc_key;
	Int_t Km_geo_assoc_key;
	Int_t pi_geo_assoc_key;
	Int_t Kp_geo_assoc_pdgId;
	Int_t Km_geo_assoc_pdgId;
	Int_t pi_geo_assoc_pdgId;

	// the gen daughters' own indices in prunedGenParticles, so the offline macro can
	// do the cross-tag compare without re-deriving them
	Int_t Gen_Kp_idx;
	Int_t Gen_Km_idx;
	Int_t Gen_pi_idx;
};

#endif
