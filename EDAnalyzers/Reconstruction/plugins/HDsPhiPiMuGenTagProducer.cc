#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Utilities/interface/EDMException.h"


#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"

#include "TMath.h"
#include "TLorentzVector.h"

#include "EDAnalyzers/Reconstruction/interface/H_phipimu_GenTag_Tree.h"

// Custom NanoAOD migration of the H_phipimu_GenTag tuplizer (background gen-Ds tagging).
// Same transformation as HDsPhiPiMuProducer: identical reconstruction, emitted as a
// pat::CompositeCandidateCollection; H_phipimu_GenTag_Tree kept as an in-memory scratch
// struct; fillCand() serialises every branch to a same-named userFloat/userInt.
class HDsPhiPiMuGenTagProducer : public edm::stream::EDProducer<>
{
public:
	explicit HDsPhiPiMuGenTagProducer(const edm::ParameterSet &);
	~HDsPhiPiMuGenTagProducer() override { delete ftree; };

	static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
	void produce(edm::Event &, const edm::EventSetup &) override;
	void fillCand(pat::CompositeCandidate &c) const;

	// Truth of one track in the K+K-pi candidate, resolved through
	// packedPFCandidateToGenAssociation. All fields are `null`/-1/false when the
	// association returns a null ref -- see the header comment on the *_gen_* branches
	// for why that is the normal case for an ordinary fragmentation track.
	struct PFTruth
	{
		Int_t pdgId    = null;   // true pdgId of the matched gen particle
		Int_t momPdgId = null;   // its immediate mother's pdgId (333 for a phi)
		Int_t momKey   = -1;     // identity of that mother within prunedGenParticles
		Int_t hfPdgId  = null;   // nearest charm/bottom hadron ancestor
		Int_t hfKey    = -1;     // identity of that ancestor
		Bool_t hfIsMom = false;  // that ancestor IS the immediate mother (direct decay)
		Bool_t valid   = false;  // the association returned a usable gen particle
	};

	static Bool_t isHeavyFlavour(const Int_t pdgId);
	// Walks up the mother chain and returns the nearest charm/bottom hadron ancestor, or
	// nullptr. Identity is carried as a GenParticleRef key, never as a pdgId: an event can
	// hold two Ds, and comparing pdgId would call tracks from different ones "same origin".
	const reco::GenParticle *nearestHeavyFlavour(const reco::GenParticle &gp, Int_t &key, Int_t depth = 0) const;
	PFTruth pfTruth(const pat::PackedCandidate &pf, const pat::PackedCandidateRef &pfRef,
	                const edm::Association<reco::GenParticleCollection> &assoc) const;
	Float_t bScore(const pat::Jet &jet) const;
	Measurement1D distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const;
	Measurement1D distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const;

	edm::EDGetTokenT<pat::PackedCandidateCollection> PFCandToken_;
	edm::EDGetTokenT<pat::MuonCollection> MuonToken_;
	edm::EDGetTokenT<pat::ElectronCollection> ElectronToken_;
	edm::EDGetTokenT<pat::METCollection> MetToken_;
	edm::EDGetTokenT<pat::JetCollection> JetToken_;
	edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PileupToken_;
	edm::EDGetTokenT<reco::VertexCollection> primvtxToken_;
	edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
	edm::EDGetTokenT<reco::GenParticleCollection> GenPartToken_;
	edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> associationToken_;
	edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

	// false (default): opposite-sign selection -- pi is required opposite the best muon,
	// which is the real H->Ds(*)W signal correlation (Ds and W have opposite charge).
	// true: same-sign control region -- pi required SAME as the best muon, which by charge
	// conservation excludes real signal and leaves wrong-charge combinatoric background.
	bool sameSignPi_;
	std::vector<std::string> bDiscNames_;
	double bScoreWPLoose_, bScoreWPMedium_, bScoreWPTight_;
	std::string eleIdCutLoose_, eleIdCutMedium_, eleIdCutTight_, eleIdMvaWP80_, eleIdMvaWP90_;

	const Float_t Mass_K    = 0.493677;
	const Float_t Mass_pi   = 0.13957039;
	const Float_t Mass_mu   = 0.1056583755;
	const Float_t Mass_phi   = 1.01946;
	const Float_t Mass_Ds   = 1.96835;
	const Float_t Mass_Dsst = 2.1066;

	Float_t Mass_K_sigma   = 0.000016;
	Float_t Mass_pi_sigma  = 0.000016;

	H_phipimu_GenTag_Tree *ftree;   // in-memory scratch struct only; no TTree, no CreateBranches()
};

HDsPhiPiMuGenTagProducer::HDsPhiPiMuGenTagProducer(const edm::ParameterSet &iConfig)
{
	edm::InputTag PFCandTag_("packedPFCandidates");
	PFCandToken_ = consumes<pat::PackedCandidateCollection>(PFCandTag_);

	edm::InputTag MuonTag_("slimmedMuons");
	MuonToken_ = consumes<pat::MuonCollection>(MuonTag_);

	edm::InputTag ElectronTag_("slimmedElectrons");
	ElectronToken_ = consumes<pat::ElectronCollection>(ElectronTag_);

	// Puppi MET + standard-NANO jet collection (JMEcorr retired; jet ID via userInt).
	edm::InputTag MetTag_("slimmedMETsPuppi");
	MetToken_ = consumes<pat::METCollection>(MetTag_);

	edm::InputTag JetTag_("updatedJetsPuppiWithUserData");
	JetToken_ = consumes<pat::JetCollection>(JetTag_);



	// Pileup. Same gating as the generator weight: declared unconditionally, read only for
	// simulation. "slimmedAddPileupInfo" is the MiniAOD label ("addPileupInfo" is AOD).
	edm::InputTag PileupTag_("slimmedAddPileupInfo");
	PileupToken_ = consumes<std::vector<PileupSummaryInfo>>(PileupTag_);

	edm::InputTag primvtxTag_("offlineSlimmedPrimaryVertices");
	primvtxToken_ = consumes<reco::VertexCollection>(primvtxTag_);

	edm::InputTag triggerTag_("TriggerResults", "", "HLT");
	triggerToken_ = consumes<edm::TriggerResults>(triggerTag_);

	edm::InputTag GenPartTag_("prunedGenParticles");
	GenPartToken_ = consumes<reco::GenParticleCollection>(GenPartTag_);

	edm::InputTag assotiationTag_("packedPFCandidateToGenAssociation");
	associationToken_ = consumes<edm::Association<reco::GenParticleCollection>>(assotiationTag_);

	ttBuilderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"));

	sameSignPi_ = iConfig.existsAs<bool>("sameSignPi") ? iConfig.getParameter<bool>("sameSignPi") : false;

	// b-tag discriminator name(s); summed per jet to form the b-score. 2022 uses
	// ParticleNet: the combined B-vs-all discriminator is a single stored value.
	// Override in the cmsRun config if the MiniAOD embeds a different label. A wrong
	// name throws at runtime (fail loud), so confirm via jet.getPairDiscri() on the sample.
	bDiscNames_ = iConfig.existsAs<std::vector<std::string>>("bDiscNames")
		? iConfig.getParameter<std::vector<std::string>>("bDiscNames")
		: std::vector<std::string>{"pfParticleNetFromMiniAODAK4PuppiCentralDiscriminatorsJetTags:BvsAll"};

	// Working points for counting n_b_jets_{loose,medium,tight} (jets with b-score above
	// each threshold): the ParticleNet L/M/T WPs of the era, from BTV btagging.json
	// (particleNet_wp_values). REQUIRED, no default -- they are era-dependent (2022:
	// 0.0470/0.2450/0.6734, 2022EE: 0.0499/0.2605/0.6915) and the config sets them per YEAR.
	// Until 2026-09-15 this read existsAs<double>(...) with the 2022 values as defaults, but
	// fillDescriptions never declared the parameters, so the config could not override them
	// at all ("Illegal parameter") and every 2022EE output silently counted n_b_jets with the
	// 2022 WPs. Declared required in fillDescriptions now, so a missing WP fails loud.
	bScoreWPLoose_  = iConfig.getParameter<double>("bScoreWPLoose");
	bScoreWPMedium_ = iConfig.getParameter<double>("bScoreWPMedium");
	bScoreWPTight_  = iConfig.getParameter<double>("bScoreWPTight");

	// Electron ID names for the ss/os electron counts (via electron.electronID(name) > 0.5).
	// Names are era-specific (2022 defaults below); a wrong name throws at runtime (fail
	// loud, lists the available IDs), so confirm against the sample or override in the config.
	auto eleId = [&](const char *key, const char *def) {
		return iConfig.existsAs<std::string>(key) ? iConfig.getParameter<std::string>(key) : std::string(def);
	};
	eleIdCutLoose_  = eleId("eleIdCutLoose",  "cutBasedElectronID-RunIIIWinter22-V1-loose");
	eleIdCutMedium_ = eleId("eleIdCutMedium", "cutBasedElectronID-RunIIIWinter22-V1-medium");
	eleIdCutTight_  = eleId("eleIdCutTight",  "cutBasedElectronID-RunIIIWinter22-V1-tight");
	eleIdMvaWP80_   = eleId("eleIdMvaWP80",   "mvaEleID-RunIIIWinter22-iso-V1-wp80");
	eleIdMvaWP90_   = eleId("eleIdMvaWP90",   "mvaEleID-RunIIIWinter22-iso-V1-wp90");

	ftree = new H_phipimu_GenTag_Tree(nullptr);   // scratch struct only

	produces<pat::CompositeCandidateCollection>();
}

void HDsPhiPiMuGenTagProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	auto out = std::make_unique<pat::CompositeCandidateCollection>();

	// Pileup: plain gen-level quantity kept per candidate. genWeight dropped (NanoAOD-native);
	// the inclusive nTrueInt_all histogram is dropped (pileup reweighting is a post-NanoAOD step).
	Float_t nTrueInt = null;
	Int_t nPUInt = null;
	if (!iEvent.isRealData())
	{
		edm::Handle<std::vector<PileupSummaryInfo>> PileupHandle;
		iEvent.getByToken(PileupToken_, PileupHandle);
		for (const PileupSummaryInfo &pu : *PileupHandle)
		{
			if (pu.getBunchCrossing() != 0) continue;   // in-time crossing only
			nTrueInt = pu.getTrueNumInteractions();
			nPUInt   = pu.getPU_NumInteractions();
			break;
		}
	}

	edm::ESHandle<TransientTrackBuilder> ttBuilderHandle = iSetup.getHandle(ttBuilderToken_);
	const TransientTrackBuilder *ttBuilder = ttBuilderHandle.product();

	ftree->Init();

	ftree->nTrueInt = nTrueInt;
	ftree->nPUInt = nPUInt;
	ftree->event_id = iEvent.id().event();
	ftree->run_id   = iEvent.id().run();
	ftree->lumi_id  = iEvent.luminosityBlock();
	Int_t cand_counter = 0;

	edm::Handle<pat::PackedCandidateCollection> PFCandHandle;
	iEvent.getByToken(PFCandToken_, PFCandHandle);

	edm::Handle<pat::MuonCollection> MuonHandle;
	iEvent.getByToken(MuonToken_, MuonHandle);

	//==================================================//
	// Event-level gen charm counts
	//==================================================//
	// How many true charm hadrons the event contained, independent of whether anything was
	// reconstructed. Together with the per-candidate tags this gives the reconstruction
	// efficiency directly. Last copies only, so a particle is not counted once per copy
	// record. MiniAOD pruning keeps every c/b hadron unconditionally (no pT, eta or origin
	// requirement), so these counts are complete.
	if (!iEvent.isRealData())
	{
		edm::Handle<reco::GenParticleCollection> GenPartHandle;
		iEvent.getByToken(GenPartToken_, GenPartHandle);
		for (const reco::GenParticle &gp : *GenPartHandle)
		{
			if (!gp.statusFlags().isLastCopy()) continue;
			const Int_t apid = std::abs(gp.pdgId());
			if (apid == 431)
			{
				ftree->n_gen_Ds++;
				for (size_t d = 0; d < gp.numberOfDaughters(); d++)
				{
					if (gp.daughter(d) && gp.daughter(d)->pdgId() == 333) { ftree->n_gen_Ds_phipi++; break; }
				}
			}
			else if (apid == 411 || apid == 421) ftree->n_gen_D++;
		}
	}

	//==================================================//
	// Triggers
	//==================================================//
	edm::Handle<edm::TriggerResults> triggerHandle;
	iEvent.getByToken(triggerToken_, triggerHandle);

	Bool_t trig_HLT_IsoMu24 = false;
	Bool_t trig_HLT_Mu50    = false;
	const edm::TriggerNames &names = iEvent.triggerNames(*triggerHandle);
	for (Size_t i = 0; i < names.size(); i++)
	{
		TString trigName = names.triggerName(i);
		Bool_t pass = triggerHandle->accept(i);
		if (trigName.Contains("HLT_IsoMu24_v")) trig_HLT_IsoMu24 = pass;
		else if (trigName.Contains("HLT_Mu50_v")) trig_HLT_Mu50 = pass;
	}
	if (!trig_HLT_IsoMu24 && !trig_HLT_Mu50)
		{ iEvent.put(std::move(out)); return; }

	ftree->trig_HLT_IsoMu24 = trig_HLT_IsoMu24;
	ftree->trig_HLT_Mu50    = trig_HLT_Mu50;

	//==================================================//
	// Primary vertex
	//==================================================//
	edm::Handle<reco::VertexCollection> primvtxHandle;
	iEvent.getByToken(primvtxToken_, primvtxHandle);
	if (primvtxHandle->empty()) { iEvent.put(std::move(out)); return; }
	const reco::Vertex &primvtx = primvtxHandle->front();
	if (!primvtx.isValid()) { iEvent.put(std::move(out)); return; }
	GlobalPoint primvtx_pos(primvtx.x(), primvtx.y(), primvtx.z());
	GlobalError primvtx_cov(primvtx.covariance());

	ftree->primvtx_isValid    = primvtx.isValid();
	ftree->primvtx_isFake     = primvtx.isFake();
	ftree->primvtx_chi2       = primvtx.chi2();
	ftree->primvtx_ndof       = primvtx.ndof();
	ftree->primvtx_chi2ndof   = primvtx.chi2() / primvtx.ndof();
	ftree->primvtx_vx         = primvtx.x();
	ftree->primvtx_vy         = primvtx.y();
	ftree->primvtx_vz         = primvtx.z();
	ftree->primvtx_vxerr      = primvtx.xError();
	ftree->primvtx_vyerr      = primvtx.yError();
	ftree->primvtx_vzerr      = primvtx.zError();

	//==================================================//
	// MET
	//==================================================//
	edm::Handle<pat::METCollection> MetHandle;
	iEvent.getByToken(MetToken_, MetHandle);
	if (MetHandle->empty()) { iEvent.put(std::move(out)); return; }
	const pat::MET &Met = MetHandle->front();

	ftree->met_sig   = Met.metSignificance();
	ftree->met_phi   = Met.phi();
	ftree->met_px    = Met.px();
	ftree->met_py    = Met.py();
	ftree->met_pt    = Met.pt();
	ftree->met_sumEt = Met.sumEt();

	//==================================================//
	// Jets
	//==================================================//
	edm::Handle<pat::JetCollection> JetHandle;
	iEvent.getByToken(JetToken_, JetHandle);

	Float_t jet_ht = 0, jet_ht_pt20 = 0, jet_leadpt = -1, jet_maxbscore = -999;
	for (Int_t idx_jet = 0; idx_jet < static_cast<Int_t>(JetHandle->size()); idx_jet++)
	{
		const pat::Jet &Jet = (*JetHandle)[idx_jet];
		if (Jet.userInt("tightId") == 0) continue;   // standard-NANO tight jet ID
		ftree->num_jets++;
		if (Jet.pt() > 20)  ftree->num_jets_pt20++;
		if (Jet.pt() > 30)  ftree->num_jets_pt30++;
		if (Jet.pt() > 40)  ftree->num_jets_pt40++;
		if (Jet.pt() > 50)  ftree->num_jets_pt50++;
		if (Jet.pt() > 60)  ftree->num_jets_pt60++;
		if (Jet.pt() > 70)  ftree->num_jets_pt70++;
		if (Jet.pt() > 80)  ftree->num_jets_pt80++;
		if (Jet.pt() > 90)  ftree->num_jets_pt90++;
		if (Jet.pt() > 100) ftree->num_jets_pt100++;

		jet_ht += Jet.pt();
		if (Jet.pt() > 20) jet_ht_pt20 += Jet.pt();
		if (Jet.pt() > jet_leadpt) jet_leadpt = Jet.pt();
		// b-score only for central jets: the ParticleNet Central tagger is defined for
		// |eta|<2.5 (forward jets lack the discriminator and would throw / can't be b-tagged).
		if (std::abs(Jet.eta()) < 2.5)
		{
			Float_t score = bScore(Jet);
			if (score > jet_maxbscore) jet_maxbscore = score;
		}
	}
	ftree->jet_ht = jet_ht;
	ftree->jet_ht_pt20 = jet_ht_pt20;
	if (ftree->num_jets > 0) ftree->jet_leadpt = jet_leadpt;
	if (jet_maxbscore > -999) ftree->jet_maxbscore = jet_maxbscore;

	//==================================================//
	// Muon
	//==================================================//
	Int_t best_mu_idx = -1;
	Float_t best_mu_pt = -1;

	for (Int_t idx_mu = 0; idx_mu < static_cast<Int_t>(MuonHandle->size()); idx_mu++)
	{
		const pat::Muon &muon = (*MuonHandle)[idx_mu];
		// Pre-cut on the UNCORRECTED pt, deliberately below the analysis cut. Was 26 until
		// 2026-09-14, then 20, now 24. The analysis requires CORRECTED pt > 26 offline (26 = MUO trigger-SF pt
		// floor, flow=error) -- every consumer of this nano MUST apply that cut, or events whose
		// best muon is below the trigger plateau leak in.
		// Why a margin is needed at all: the official MuonScaRe nominal correction leaves
		// pt < 26 untouched (pt_scale / pt_resol low_pt_threshold = 26), so nominal cannot move
		// a muon up across 26; but pt_scale_var (the scale systematic) has no threshold and
		// shifts any pt by up to 0.028 GeV at pt~26, so a 26 pre-cut would under-populate the
		// scale-up variation. 24 GeV (set 2026-09-15) leaves a 2 GeV margin and matches the
		// HLT_IsoMu24 threshold. If MUO ever lowers low_pt_threshold below 26, nominal corrections
		// would move sub-26 muons up too (up to ~+8.6% seen above 26) and this margin must be revisited.
		if (muon.pt() < 24) continue;
		if (!muon.isTightMuon(primvtx)) continue;
		if (muon.pt() > best_mu_pt) { best_mu_idx = idx_mu; best_mu_pt = muon.pt(); }
	}
	if (best_mu_idx == -1) { iEvent.put(std::move(out)); return; }
	const pat::Muon &mu = (*MuonHandle)[best_mu_idx];
	TLorentzVector mu_tlv;
	mu_tlv.SetXYZM(mu.px(), mu.py(), mu.pz(), Mass_mu);

	std::pair<Bool_t, Measurement1D> mu_ipresult   = IPTools::absoluteImpactParameter3D((*ttBuilder).build(mu.muonBestTrack()), primvtx);
	std::pair<Bool_t, Measurement1D> mu_ipxyresult = IPTools::absoluteTransverseImpactParameter((*ttBuilder).build(mu.muonBestTrack()), primvtx);

	Float_t mu_IsoR03_PFIso = mu.pfIsolationR03().sumChargedHadronPt + std::max(0.0f, mu.pfIsolationR03().sumNeutralHadronEt + mu.pfIsolationR03().sumPhotonEt - 0.5f * mu.pfIsolationR03().sumPUPt);
	Float_t mu_IsoR04_PFIso = mu.pfIsolationR04().sumChargedHadronPt + std::max(0.0f, mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5f * mu.pfIsolationR04().sumPUPt);

	// The charge the candidate loop below requires of the pi: opposite the best muon in the
	// default OS selection, same as it in the SS control region. Everything charge-related
	// keys off this one value, so os/ss stay defined relative to the pi in both modes.
	const Int_t pi_charge_ref = sameSignPi_ ? mu.charge() : -mu.charge();

	// ss/os are relative to the pi charge, so they mean the same thing in OS and SS mode
	// (in OS "same as pi" is opposite the best muon; in SS it is the same as the best muon).
	for (Int_t idx_mu = 0; idx_mu < static_cast<Int_t>(MuonHandle->size()); idx_mu++)
	{
		const pat::Muon &muon = (*MuonHandle)[idx_mu];
		Bool_t os = (muon.charge() == -pi_charge_ref);  // opposite pi
		Bool_t ss = (muon.charge() ==  pi_charge_ref);  // same as pi
		if (muon.isLooseMuon())        { if (os) ftree->num_os_mu_loose++;  else if (ss) ftree->num_ss_mu_loose++;  }
		if (muon.isMediumMuon())       { if (os) ftree->num_os_mu_medium++; else if (ss) ftree->num_ss_mu_medium++; }
		if (muon.isTightMuon(primvtx)) { if (os) ftree->num_os_mu_tight++;  else if (ss) ftree->num_ss_mu_tight++;  }
	}

	// Electron ss/os counts (same charge convention relative to pi as the muons above).
	// Counting only -- electrons do not enter the event selection.
	edm::Handle<pat::ElectronCollection> ElectronHandle;
	iEvent.getByToken(ElectronToken_, ElectronHandle);
	for (Int_t idx_e = 0; idx_e < static_cast<Int_t>(ElectronHandle->size()); idx_e++)
	{
		const pat::Electron &ele = (*ElectronHandle)[idx_e];
		Bool_t os = (ele.charge() == -pi_charge_ref);
		Bool_t ss = (ele.charge() ==  pi_charge_ref);
		if (ele.electronID(eleIdCutLoose_)  > 0.5f) { if (os) ftree->num_os_e_cutloose++;  else if (ss) ftree->num_ss_e_cutloose++;  }
		if (ele.electronID(eleIdCutMedium_) > 0.5f) { if (os) ftree->num_os_e_cutmedium++; else if (ss) ftree->num_ss_e_cutmedium++; }
		if (ele.electronID(eleIdCutTight_)  > 0.5f) { if (os) ftree->num_os_e_cuttight++;  else if (ss) ftree->num_ss_e_cuttight++;  }
		if (ele.electronID(eleIdMvaWP80_)   > 0.5f) { if (os) ftree->num_os_e_mvawp80++;   else if (ss) ftree->num_ss_e_mvawp80++;   }
		if (ele.electronID(eleIdMvaWP90_)   > 0.5f) { if (os) ftree->num_os_e_mvawp90++;   else if (ss) ftree->num_ss_e_mvawp90++;   }
	}

	ftree->mu_charge   = mu.charge();
	ftree->mu_eta      = mu.eta();   ftree->mu_phi  = mu.phi();
	ftree->mu_px       = mu.px();    ftree->mu_py   = mu.py();   ftree->mu_pz = mu.pz();
	ftree->mu_pt       = mu.pt();    ftree->mu_p    = mu.p();
	ftree->mu_isHighPt = mu.isHighPtMuon(primvtx);
	ftree->mu_isLoose  = mu.isLooseMuon();
	ftree->mu_isMedium = mu.isMediumMuon();
	ftree->mu_isSoft   = mu.isSoftMuon(primvtx);
	ftree->mu_isTight  = mu.isTightMuon(primvtx);
	ftree->mu_isPF     = mu.isPFMuon();
	ftree->mu_isTracker = mu.isTrackerMuon();
	ftree->mu_isGlobal  = mu.isGlobalMuon();
	if (mu_ipresult.first)
	{
		ftree->mu_ip    = mu_ipresult.second.value();
		ftree->mu_iperr = mu_ipresult.second.error();
		ftree->mu_ipsig = mu_ipresult.second.significance();
	}
	if (mu_ipxyresult.first)
	{
		ftree->mu_ipxy    = mu_ipxyresult.second.value();
		ftree->mu_ipxyerr = mu_ipxyresult.second.error();
		ftree->mu_ipxysig = mu_ipxyresult.second.significance();
	}
	ftree->mu_IsoR03_sumChargedHadronPt = mu.pfIsolationR03().sumChargedHadronPt;
	ftree->mu_IsoR03_sumNeutralHadronEt = mu.pfIsolationR03().sumNeutralHadronEt;
	ftree->mu_IsoR03_sumPhotonEt        = mu.pfIsolationR03().sumPhotonEt;
	ftree->mu_IsoR03_sumPUPt            = mu.pfIsolationR03().sumPUPt;
	ftree->mu_IsoR03_PFIso              = mu_IsoR03_PFIso;
	ftree->mu_IsoR03_RelIso             = mu_IsoR03_PFIso / mu.pt();
	ftree->mu_IsoR04_sumChargedHadronPt = mu.pfIsolationR04().sumChargedHadronPt;
	ftree->mu_IsoR04_sumNeutralHadronEt = mu.pfIsolationR04().sumNeutralHadronEt;
	ftree->mu_IsoR04_sumPhotonEt        = mu.pfIsolationR04().sumPhotonEt;
	ftree->mu_IsoR04_sumPUPt            = mu.pfIsolationR04().sumPUPt;
	ftree->mu_IsoR04_PFIso              = mu_IsoR04_PFIso;
	ftree->mu_IsoR04_RelIso             = mu_IsoR04_PFIso / mu.pt();
	ftree->W_tmass = std::sqrt(2.0f * mu.pt() * Met.pt() * (1.0f - std::cos(mu.phi() - Met.phi())));

	//==================================================//
	// PF hadron candidates
	//==================================================//
	std::vector<Int_t> idx_pf;
	for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[i];
		if (!pf.hasTrackDetails())  continue;
		if (!pf.trackHighPurity())  continue;
		if (pf.pt() < 1)            continue;
		if (abs(pf.pdgId()) != 211) continue;
		idx_pf.push_back(i);
	}
	if (idx_pf.size() < 3) { iEvent.put(std::move(out)); return; }

	// Truth of every track in the pool, resolved ONCE here and cached. This must not move
	// inside the K+/K-/pi loops below: those are O(n^3) over this pool, and an ancestry walk
	// per iteration would dominate the runtime. The loops only ever index this vector.
	// Left all-default on data, where neither the association nor prunedGenParticles exists.
	std::vector<PFTruth> pf_truth(idx_pf.size());
	if (!iEvent.isRealData())
	{
		edm::Handle<edm::Association<reco::GenParticleCollection>> associationHandle;
		iEvent.getByToken(associationToken_, associationHandle);
		for (size_t i = 0; i < idx_pf.size(); i++)
		{
			const pat::PackedCandidateRef pfRef(PFCandHandle, idx_pf[i]);
			pf_truth[i] = pfTruth((*PFCandHandle)[idx_pf[i]], pfRef, *associationHandle);
		}
	}

	// Photon candidates
	std::vector<Int_t> idx_photon;
	for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
	{
		if (abs((*PFCandHandle)[i].pdgId()) != 22) continue;
		if ((*PFCandHandle)[i].pt() < 0.5) continue;
		idx_photon.push_back(i);
	}

	KinematicParticleFactoryFromTransientTrack factory;
	KinematicParticleVertexFitter KPVFitter;
	float chi = 0, ndf = 0;

	//==================================================//
	// K+ loop
	//==================================================//
	for (Int_t iKp = 0; iKp < static_cast<Int_t>(idx_pf.size()); iKp++)
	{
		const Int_t idx_Kp = idx_pf[iKp];
		const pat::PackedCandidate &Kp_pf = (*PFCandHandle)[idx_Kp];
		if (Kp_pf.charge() == mu.charge()) continue;
		if (Kp_pf.pt() < 2) continue;

		TLorentzVector Kp_p4;
		Kp_p4.SetXYZM(Kp_pf.px(), Kp_pf.py(), Kp_pf.pz(), Mass_K);
		reco::TransientTrack Kp_tt = (*ttBuilder).build(Kp_pf.bestTrack());
		Int_t Kp_charge   = Kp_pf.charge();
		Int_t Kp_PVassoc  = Kp_pf.pvAssociationQuality();

		//==================================================//
		// K- loop
		//==================================================//
		for (Int_t iKm = 0; iKm < static_cast<Int_t>(idx_pf.size()); iKm++)
		{
			if (iKp == iKm) continue;
			const Int_t idx_Km = idx_pf[iKm];
			const pat::PackedCandidate &Km_pf = (*PFCandHandle)[idx_Km];
			if (Km_pf.charge() != mu.charge()) continue;
			if (Km_pf.pt() < 2) continue;

			TLorentzVector Km_p4;
			Km_p4.SetXYZM(Km_pf.px(), Km_pf.py(), Km_pf.pz(), Mass_K);
			reco::TransientTrack Km_tt = (*ttBuilder).build(Km_pf.bestTrack());
			Int_t Km_charge  = Km_pf.charge();
			Int_t Km_PVassoc = Km_pf.pvAssociationQuality();

			//==================================================//
			// pi loop
			//==================================================//
			for (Int_t ipi = 0; ipi < static_cast<Int_t>(idx_pf.size()); ipi++)
			{
				if ((iKp == ipi) || (iKm == ipi)) continue;
				const Int_t idx_pi = idx_pf[ipi];
				const pat::PackedCandidate &pi_pf = (*PFCandHandle)[idx_pi];
				// OS (default): pi opposite the best muon, the real signal correlation.
				// SS: pi same as the best muon, which excludes real signal by charge
				// conservation and leaves wrong-charge combinatoric background.
				// (The Kp/Km charge cuts above are deliberately NOT switched -- they only
				// decide which K gets labelled Kp vs Km, orthogonal to the pi choice.)
				if (pi_pf.charge() != pi_charge_ref) continue;

				TLorentzVector pi_p4;
				pi_p4.SetXYZM(pi_pf.px(), pi_pf.py(), pi_pf.pz(), Mass_pi);
				if ((Kp_p4 + Km_p4 + pi_p4).Pt() < 5) continue;
				// pre-fit KKpi window: loose pre-selection on the raw mass, kept wider than the
				// post-fit Ds window (below) so it never binds the SR/CR sideband.
				if (std::abs((Kp_p4 + Km_p4 + pi_p4).M() - Mass_Ds) > 0.5) continue;
				// pre-fit phi window: loose upper cut on the raw KK mass (wider than the post-fit
				// +50 MeV) to drop non-phi KK pairs before the expensive vertex fit. Pure speedup.
				if ((Kp_p4 + Km_p4).M() > Mass_phi + 0.08) continue;

				reco::TransientTrack pi_tt = (*ttBuilder).build(pi_pf.bestTrack());
				Int_t pi_charge  = pi_pf.charge();
				Int_t pi_PVassoc = pi_pf.pvAssociationQuality();

				//--------------------------------------------------//
				// Kinematic fit on Ds
				//--------------------------------------------------//
				std::vector<RefCountedKinematicParticle> DsDaughters;
				DsDaughters.push_back(factory.particle(Kp_tt, Mass_K,  chi, ndf, Mass_K_sigma));
				DsDaughters.push_back(factory.particle(Km_tt, Mass_K,  chi, ndf, Mass_K_sigma));
				DsDaughters.push_back(factory.particle(pi_tt, Mass_pi, chi, ndf, Mass_pi_sigma));
				RefCountedKinematicTree DsTree = KPVFitter.fit(DsDaughters);

				if (!DsTree->isValid()) continue;
				DsTree->movePointerToTheTop();
				RefCountedKinematicVertex Ds_vertex = DsTree->currentDecayVertex();
				if (!Ds_vertex->vertexIsValid()) continue;
				GlobalPoint Ds_vtx_pos = Ds_vertex->position();
				GlobalError Ds_vtx_cov = Ds_vertex->error();

				if (TMath::Prob(Ds_vertex->chiSquared(), Ds_vertex->degreesOfFreedom()) < 0.001) continue;

				RefCountedKinematicParticle Ds = DsTree->currentParticle();
				GlobalVector Ds_mtm = Ds->currentState().globalMomentum();
				const reco::TransientTrack Ds_rett = Ds->refittedTransientTrack();

				// post-fit Ds window: wide enough to hold the SR + Ds-mass sidebands (CR);
				// the SR/CR boundary is applied offline on Ds_mass.
				if (std::abs(Ds->currentState().mass() - Mass_Ds) > 0.3) continue;
				if (Ds_mtm.perp() < 7) continue;

				GlobalVector displacement = Ds_vtx_pos - primvtx_pos;
				Float_t dira = displacement.dot(Ds_mtm) / (displacement.mag() * Ds_mtm.mag());
				if (dira < 0.99) continue;

				TLorentzVector Ds_tlv;
				Ds_tlv.SetXYZM(Ds_mtm.x(), Ds_mtm.y(), Ds_mtm.z(), Ds->currentState().mass());

				TLorentzVector Dsmu_tlv = Ds_tlv + mu_tlv;
				if (Dsmu_tlv.M() > 125) continue;

				DsTree->movePointerToTheFirstChild();
				RefCountedKinematicParticle Kp = DsTree->currentParticle();
				GlobalVector Kp_mtm = Kp->currentState().globalMomentum();
				const reco::TransientTrack Kp_rett = Kp->refittedTransientTrack();

				DsTree->movePointerToTheNextChild();
				RefCountedKinematicParticle Km = DsTree->currentParticle();
				GlobalVector Km_mtm = Km->currentState().globalMomentum();
				const reco::TransientTrack Km_rett = Km->refittedTransientTrack();

				DsTree->movePointerToTheNextChild();
				RefCountedKinematicParticle pi = DsTree->currentParticle();
				GlobalVector pi_mtm = pi->currentState().globalMomentum();
				const reco::TransientTrack pi_rett = pi->refittedTransientTrack();

				TLorentzVector Kp_tlv, Km_tlv, pi_tlv, phi_tlv;
				Kp_tlv.SetXYZM(Kp_mtm.x(), Kp_mtm.y(), Kp_mtm.z(), Mass_K);
				Km_tlv.SetXYZM(Km_mtm.x(), Km_mtm.y(), Km_mtm.z(), Mass_K);
				pi_tlv.SetXYZM(pi_mtm.x(), pi_mtm.y(), pi_mtm.z(), Mass_pi);
				phi_tlv = Kp_tlv + Km_tlv;
				if (phi_tlv.Pt() < 3) continue;
				// phi window: upper cut only. The left side is physically bound by the KK
				// threshold (2*m_K = m_phi - 32 MeV), so no lower cut is needed. Upper +50 MeV.
				if (phi_tlv.M() > Mass_phi + 0.05) continue;

				//--------------------------------------------------//
				// IP
				//--------------------------------------------------//
				std::pair<Bool_t, Measurement1D> Kp_ipresult   = IPTools::absoluteImpactParameter3D(Kp_rett, primvtx);
				std::pair<Bool_t, Measurement1D> Km_ipresult   = IPTools::absoluteImpactParameter3D(Km_rett, primvtx);
				std::pair<Bool_t, Measurement1D> pi_ipresult   = IPTools::absoluteImpactParameter3D(pi_rett, primvtx);
				std::pair<Bool_t, Measurement1D> Ds_ipresult   = IPTools::absoluteImpactParameter3D(Ds_rett, primvtx);
				std::pair<Bool_t, Measurement1D> Kp_ipxyresult = IPTools::absoluteTransverseImpactParameter(Kp_rett, primvtx);
				std::pair<Bool_t, Measurement1D> Km_ipxyresult = IPTools::absoluteTransverseImpactParameter(Km_rett, primvtx);
				std::pair<Bool_t, Measurement1D> pi_ipxyresult = IPTools::absoluteTransverseImpactParameter(pi_rett, primvtx);
				std::pair<Bool_t, Measurement1D> Ds_ipxyresult = IPTools::absoluteTransverseImpactParameter(Ds_rett, primvtx);

				Measurement1D Ds_FD_XY = distanceXY(primvtx_pos, primvtx_cov, Ds_vtx_pos, Ds_vtx_cov);
				Measurement1D Ds_FD_3D = distance3D(primvtx_pos, primvtx_cov, Ds_vtx_pos, Ds_vtx_cov);

				//--------------------------------------------------//
				// Ds isolation
				//--------------------------------------------------//
				Float_t Ds_IsoR03_sumChHadPt = 0, Ds_IsoR03_sumNeuHadEt = 0, Ds_IsoR03_sumPhoEt = 0, Ds_IsoR03_sumPUPt = 0;
				Float_t Ds_IsoR04_sumChHadPt = 0, Ds_IsoR04_sumNeuHadEt = 0, Ds_IsoR04_sumPhoEt = 0, Ds_IsoR04_sumPUPt = 0;
				for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
				{
					if (i == idx_Kp || i == idx_Km || i == idx_pi) continue;
					const pat::PackedCandidate &pf = (*PFCandHandle)[i];
					Float_t iso_dR = reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), pf.eta(), pf.phi());
					if (iso_dR < 0.4)
					{
						if (std::abs(pf.pdgId()) == 211 && pf.fromPV() > 1) Ds_IsoR04_sumChHadPt  += pf.pt();
						else if (pf.pdgId() == 130)                          Ds_IsoR04_sumNeuHadEt += pf.et();
						else if (pf.pdgId() == 22)                           Ds_IsoR04_sumPhoEt    += pf.et();
						if (pf.fromPV() <= 1) Ds_IsoR04_sumPUPt += pf.pt();
						if (iso_dR < 0.3)
						{
							if (std::abs(pf.pdgId()) == 211 && pf.fromPV() > 1) Ds_IsoR03_sumChHadPt  += pf.pt();
							else if (pf.pdgId() == 130)                          Ds_IsoR03_sumNeuHadEt += pf.et();
							else if (pf.pdgId() == 22)                           Ds_IsoR03_sumPhoEt    += pf.et();
							if (pf.fromPV() <= 1) Ds_IsoR03_sumPUPt += pf.pt();
						}
					}
				}
				Float_t Ds_IsoR03_PFIso = Ds_IsoR03_sumChHadPt + std::max(0.0f, Ds_IsoR03_sumNeuHadEt + Ds_IsoR03_sumPhoEt - 0.5f * Ds_IsoR03_sumPUPt);
				Float_t Ds_IsoR03_RelIso = Ds_IsoR03_PFIso / Ds_mtm.perp();
				Float_t Ds_IsoR04_PFIso = Ds_IsoR04_sumChHadPt + std::max(0.0f, Ds_IsoR04_sumNeuHadEt + Ds_IsoR04_sumPhoEt - 0.5f * Ds_IsoR04_sumPUPt);
				Float_t Ds_IsoR04_RelIso = Ds_IsoR04_PFIso / Ds_mtm.perp();

				if (Ds_IsoR04_RelIso > 1.0f) continue;

				//--------------------------------------------------//
				// Unconstrained gamma summary (no iso/mass cut)
				//--------------------------------------------------//
				{
					const Float_t dM_nominal = Mass_Dsst - Mass_Ds;
					Int_t uc_best_idx = -1;
					Float_t uc_best_dist = 999.0f;
					Int_t n_window = 0;
					for (Int_t igamma = 0; igamma < static_cast<Int_t>(idx_photon.size()); igamma++)
					{
						const Int_t idx_g = idx_photon[igamma];
						const pat::PackedCandidate &gpf = (*PFCandHandle)[idx_g];
						TLorentzVector g_tlv;
						g_tlv.SetXYZM(gpf.px(), gpf.py(), gpf.pz(), 0.0f);
						Float_t dM = (Ds_tlv + g_tlv).M() - Ds_tlv.M();
						if (dM < 0.05f || dM > 0.25f) continue;
						n_window++;
						Float_t dist = std::abs(dM - dM_nominal);
						if (dist < uc_best_dist) { uc_best_idx = idx_g; uc_best_dist = dist; }
					}
					ftree->n_gamma_window = n_window;
					if (uc_best_idx >= 0)
					{
						const pat::PackedCandidate &gpf = (*PFCandHandle)[uc_best_idx];
						TLorentzVector g_tlv;
						g_tlv.SetXYZM(gpf.px(), gpf.py(), gpf.pz(), 0.0f);
						ftree->gamma_uc_dM     = (Ds_tlv + g_tlv).M() - Ds_tlv.M();
						ftree->gamma_uc_pt     = g_tlv.Pt();
						ftree->gamma_uc_dR_Ds  = reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), gpf.eta(), gpf.phi());
					}
				}

				//--------------------------------------------------//
				// Fill Ds-level branches (common to all gamma rows)
				//--------------------------------------------------//
				ftree->cand_id = cand_counter++;

				// Jet nearest the Ds (targets b->Ds: is the Ds inside a b jet?). Central jets
				// only (|eta|<2.5), where the ParticleNet b-score is defined; the Ds is central.
				{
					Float_t near_dR = 999, near_bscore = null;
					for (Int_t idx_jet = 0; idx_jet < static_cast<Int_t>(JetHandle->size()); idx_jet++)
					{
						const pat::Jet &Jet = (*JetHandle)[idx_jet];
						if (Jet.userInt("tightId") == 0) continue;
						if (std::abs(Jet.eta()) >= 2.5) continue;
						Float_t dR = reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), Jet.eta(), Jet.phi());
						if (dR < near_dR) { near_dR = dR; near_bscore = bScore(Jet); }
					}
					if (near_dR < 999) { ftree->Ds_nearestjet_dR = near_dR; ftree->Ds_nearestjet_bscore = near_bscore; }
				}

				// Central b-tagged jets excluding jets within dR03/dR04 of this Ds (avoids
				// counting the jet the Ds itself is clustered into as an extra b-jet).
				// Per-candidate: recomputed for each Ds using this row's Ds_mtm.
				{
					Int_t nl03 = 0, nl04 = 0, nm03 = 0, nm04 = 0, nt03 = 0, nt04 = 0;
					Int_t nl03pt = 0, nl04pt = 0, nm03pt = 0, nm04pt = 0, nt03pt = 0, nt04pt = 0;
					for (Int_t idx_jet = 0; idx_jet < static_cast<Int_t>(JetHandle->size()); idx_jet++)
					{
						const pat::Jet &Jet = (*JetHandle)[idx_jet];
						if (Jet.userInt("tightId") == 0) continue;
						if (std::abs(Jet.eta()) >= 2.5) continue;
						Float_t dR = reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), Jet.eta(), Jet.phi());
						Bool_t pastDR03 = dR > 0.3f, pastDR04 = dR > 0.4f;
						if (!pastDR03 && !pastDR04) continue;
						Float_t score = bScore(Jet);
						Bool_t isLoose = score > bScoreWPLoose_, isMedium = score > bScoreWPMedium_, isTight = score > bScoreWPTight_;
						Bool_t pt20 = Jet.pt() > 20;
						if (pastDR03)
						{
							if (isLoose)  { nl03++; if (pt20) nl03pt++; }
							if (isMedium) { nm03++; if (pt20) nm03pt++; }
							if (isTight)  { nt03++; if (pt20) nt03pt++; }
						}
						if (pastDR04)
						{
							if (isLoose)  { nl04++; if (pt20) nl04pt++; }
							if (isMedium) { nm04++; if (pt20) nm04pt++; }
							if (isTight)  { nt04++; if (pt20) nt04pt++; }
						}
					}
					ftree->n_b_jets_loose_dR03 = nl03;   ftree->n_b_jets_loose_dR04 = nl04;
					ftree->n_b_jets_medium_dR03 = nm03;  ftree->n_b_jets_medium_dR04 = nm04;
					ftree->n_b_jets_tight_dR03 = nt03;   ftree->n_b_jets_tight_dR04 = nt04;
					ftree->n_b_jets_loose_pt20_dR03 = nl03pt;   ftree->n_b_jets_loose_pt20_dR04 = nl04pt;
					ftree->n_b_jets_medium_pt20_dR03 = nm03pt;  ftree->n_b_jets_medium_pt20_dR04 = nm04pt;
					ftree->n_b_jets_tight_pt20_dR03 = nt03pt;   ftree->n_b_jets_tight_pt20_dR04 = nt04pt;
				}

				ftree->Kp_charge  = Kp_charge;
				ftree->Kp_eta     = Kp_mtm.eta();  ftree->Kp_phi  = Kp_mtm.phi();
				ftree->Kp_px      = Kp_mtm.x();    ftree->Kp_py   = Kp_mtm.y();   ftree->Kp_pz   = Kp_mtm.z();
				ftree->Kp_pt      = Kp_mtm.perp(); ftree->Kp_p    = Kp_mtm.mag(); ftree->Kp_mass = Kp->currentState().mass();
				ftree->Kp_PVassoc = Kp_PVassoc;

				ftree->Km_charge  = Km_charge;
				ftree->Km_eta     = Km_mtm.eta();  ftree->Km_phi  = Km_mtm.phi();
				ftree->Km_px      = Km_mtm.x();    ftree->Km_py   = Km_mtm.y();   ftree->Km_pz   = Km_mtm.z();
				ftree->Km_pt      = Km_mtm.perp(); ftree->Km_p    = Km_mtm.mag(); ftree->Km_mass = Km->currentState().mass();
				ftree->Km_PVassoc = Km_PVassoc;

				ftree->pi_charge  = pi_charge;
				ftree->pi_eta     = pi_mtm.eta();  ftree->pi_phi  = pi_mtm.phi();
				ftree->pi_px      = pi_mtm.x();    ftree->pi_py   = pi_mtm.y();   ftree->pi_pz   = pi_mtm.z();
				ftree->pi_pt      = pi_mtm.perp(); ftree->pi_p    = pi_mtm.mag(); ftree->pi_mass = pi->currentState().mass();
				ftree->pi_PVassoc = pi_PVassoc;

				ftree->phi_charge = Kp_charge + Km_charge;
				ftree->phi_eta    = phi_tlv.Eta();  ftree->phi_phi  = phi_tlv.Phi();
				ftree->phi_px     = phi_tlv.Px();   ftree->phi_py   = phi_tlv.Py();   ftree->phi_pz   = phi_tlv.Pz();
				ftree->phi_pt     = phi_tlv.Pt();   ftree->phi_p    = phi_tlv.P();    ftree->phi_mass = phi_tlv.M();

				ftree->Ds_charge = Kp_charge + Km_charge + pi_charge;
				ftree->Ds_eta    = Ds_mtm.eta();   ftree->Ds_phi  = Ds_mtm.phi();
				ftree->Ds_px     = Ds_mtm.x();     ftree->Ds_py   = Ds_mtm.y();     ftree->Ds_pz   = Ds_mtm.z();
				ftree->Ds_pt     = Ds_mtm.perp();  ftree->Ds_p    = Ds_mtm.mag();   ftree->Ds_mass = Ds->currentState().mass();

				ftree->Dsmu_charge = Kp_charge + Km_charge + pi_charge + mu.charge();
				ftree->Dsmu_eta    = Dsmu_tlv.Eta();  ftree->Dsmu_phi  = Dsmu_tlv.Phi();
				ftree->Dsmu_px     = Dsmu_tlv.Px();   ftree->Dsmu_py   = Dsmu_tlv.Py();   ftree->Dsmu_pz   = Dsmu_tlv.Pz();
				ftree->Dsmu_pt     = Dsmu_tlv.Pt();   ftree->Dsmu_p    = Dsmu_tlv.P();    ftree->Dsmu_mass = Dsmu_tlv.M();
				ftree->Dsmu_Ht     = Kp_mtm.perp() + Km_mtm.perp() + pi_mtm.perp() + mu_tlv.Pt();

				ftree->Ds_vx = Ds_vtx_pos.x(); ftree->Ds_vy = Ds_vtx_pos.y(); ftree->Ds_vz = Ds_vtx_pos.z();
				ftree->Ds_vxerr = Ds_vtx_cov.cxx(); ftree->Ds_vyerr = Ds_vtx_cov.cyy(); ftree->Ds_vzerr = Ds_vtx_cov.czz();
				ftree->Ds_chi2    = Ds_vertex->chiSquared();
				ftree->Ds_ndf     = Ds_vertex->degreesOfFreedom();
				ftree->Ds_chi2ndf = Ds_vertex->chiSquared() / Ds_vertex->degreesOfFreedom();
				ftree->Ds_prob    = TMath::Prob(Ds_vertex->chiSquared(), Ds_vertex->degreesOfFreedom());

				ftree->dR_Kp_Km = reco::deltaR(Kp_mtm, Km_mtm); ftree->dR_Kp_pi = reco::deltaR(Kp_mtm, pi_mtm);
				ftree->dR_Kp_Ds = reco::deltaR(Kp_mtm, Ds_mtm); ftree->dR_Km_pi = reco::deltaR(Km_mtm, pi_mtm);
				ftree->dR_Km_Ds = reco::deltaR(Km_mtm, Ds_mtm); ftree->dR_pi_Ds = reco::deltaR(pi_mtm, Ds_mtm);
				ftree->dR_phi_pi = reco::deltaR(phi_tlv.Eta(), phi_tlv.Phi(), pi_mtm.eta(), pi_mtm.phi());
				ftree->dR_phi_Ds = reco::deltaR(phi_tlv.Eta(), phi_tlv.Phi(), Ds_mtm.eta(), Ds_mtm.phi());
				ftree->dR_Kp_mu = reco::deltaR(Kp_mtm.eta(), Kp_mtm.phi(), mu.eta(), mu.phi());
				ftree->dR_Km_mu = reco::deltaR(Km_mtm.eta(), Km_mtm.phi(), mu.eta(), mu.phi());
				ftree->dR_pi_mu = reco::deltaR(pi_mtm.eta(), pi_mtm.phi(), mu.eta(), mu.phi());
				ftree->dR_phi_mu = reco::deltaR(phi_tlv.Eta(), phi_tlv.Phi(), mu.eta(), mu.phi());
				ftree->dR_Ds_mu = reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), mu.eta(), mu.phi());

				if (Kp_ipresult.first)   { ftree->Kp_ip    = Kp_ipresult.second.value();   ftree->Kp_iperr    = Kp_ipresult.second.error();   ftree->Kp_ipsig    = Kp_ipresult.second.significance();   }
				if (Kp_ipxyresult.first) { ftree->Kp_ipxy  = Kp_ipxyresult.second.value(); ftree->Kp_ipxyerr  = Kp_ipxyresult.second.error(); ftree->Kp_ipxysig  = Kp_ipxyresult.second.significance(); }
				if (Km_ipresult.first)   { ftree->Km_ip    = Km_ipresult.second.value();   ftree->Km_iperr    = Km_ipresult.second.error();   ftree->Km_ipsig    = Km_ipresult.second.significance();   }
				if (Km_ipxyresult.first) { ftree->Km_ipxy  = Km_ipxyresult.second.value(); ftree->Km_ipxyerr  = Km_ipxyresult.second.error(); ftree->Km_ipxysig  = Km_ipxyresult.second.significance(); }
				if (pi_ipresult.first)   { ftree->pi_ip    = pi_ipresult.second.value();   ftree->pi_iperr    = pi_ipresult.second.error();   ftree->pi_ipsig    = pi_ipresult.second.significance();   }
				if (pi_ipxyresult.first) { ftree->pi_ipxy  = pi_ipxyresult.second.value(); ftree->pi_ipxyerr  = pi_ipxyresult.second.error(); ftree->pi_ipxysig  = pi_ipxyresult.second.significance(); }
				if (Ds_ipresult.first)   { ftree->Ds_ip    = Ds_ipresult.second.value();   ftree->Ds_iperr    = Ds_ipresult.second.error();   ftree->Ds_ipsig    = Ds_ipresult.second.significance();   }
				if (Ds_ipxyresult.first) { ftree->Ds_ipxy  = Ds_ipxyresult.second.value(); ftree->Ds_ipxyerr  = Ds_ipxyresult.second.error(); ftree->Ds_ipxysig  = Ds_ipxyresult.second.significance(); }

				ftree->dira      = dira;
				ftree->Ds_FDxy   = Ds_FD_XY.value(); ftree->Ds_FDxyerr = Ds_FD_XY.error(); ftree->Ds_FDxysig = Ds_FD_XY.significance();
				ftree->Ds_FD     = Ds_FD_3D.value(); ftree->Ds_FDerr   = Ds_FD_3D.error(); ftree->Ds_FDsig   = Ds_FD_3D.significance();

				ftree->Ds_IsoR03_sumChargedHadronPt = Ds_IsoR03_sumChHadPt;
				ftree->Ds_IsoR03_sumNeutralHadronEt = Ds_IsoR03_sumNeuHadEt;
				ftree->Ds_IsoR03_sumPhotonEt        = Ds_IsoR03_sumPhoEt;
				ftree->Ds_IsoR03_sumPUPt            = Ds_IsoR03_sumPUPt;
				ftree->Ds_IsoR03_PFIso              = Ds_IsoR03_PFIso;
				ftree->Ds_IsoR03_RelIso             = Ds_IsoR03_RelIso;
				ftree->Ds_IsoR04_sumChargedHadronPt = Ds_IsoR04_sumChHadPt;
				ftree->Ds_IsoR04_sumNeutralHadronEt = Ds_IsoR04_sumNeuHadEt;
				ftree->Ds_IsoR04_sumPhotonEt        = Ds_IsoR04_sumPhoEt;
				ftree->Ds_IsoR04_sumPUPt            = Ds_IsoR04_sumPUPt;
				ftree->Ds_IsoR04_PFIso              = Ds_IsoR04_PFIso;
				ftree->Ds_IsoR04_RelIso             = Ds_IsoR04_RelIso;

				//--------------------------------------------------//
				// Truth tagging of this (Kp, Km, pi) combination
				//--------------------------------------------------//
				const PFTruth &tKp = pf_truth[iKp];
				const PFTruth &tKm = pf_truth[iKm];
				const PFTruth &tpi = pf_truth[ipi];

				ftree->Kp_gen_pdgId = tKp.pdgId; ftree->Kp_gen_mom_pdgId = tKp.momPdgId;
				ftree->Kp_gen_mom_key = tKp.momKey; ftree->Kp_gen_hf_pdgId = tKp.hfPdgId;
				ftree->Kp_gen_hf_key = tKp.hfKey; ftree->Kp_gen_hf_isMom = tKp.hfIsMom;
				ftree->Km_gen_pdgId = tKm.pdgId; ftree->Km_gen_mom_pdgId = tKm.momPdgId;
				ftree->Km_gen_mom_key = tKm.momKey; ftree->Km_gen_hf_pdgId = tKm.hfPdgId;
				ftree->Km_gen_hf_key = tKm.hfKey; ftree->Km_gen_hf_isMom = tKm.hfIsMom;
				ftree->pi_gen_pdgId = tpi.pdgId; ftree->pi_gen_mom_pdgId = tpi.momPdgId;
				ftree->pi_gen_mom_key = tpi.momKey; ftree->pi_gen_hf_pdgId = tpi.hfPdgId;
				ftree->pi_gen_hf_key = tpi.hfKey; ftree->pi_gen_hf_isMom = tpi.hfIsMom;

				ftree->cand_gen_nmatched = (tKp.valid ? 1 : 0) + (tKm.valid ? 1 : 0) + (tpi.valid ? 1 : 0);

				// Distinct heavy-flavour ancestors, counted by KEY. Using pdgId here would
				// merge two different Ds in the same event into one "common origin".
				const Int_t hfKeys[3] = {tKp.hfKey, tKm.hfKey, tpi.hfKey};
				Int_t n_hf = 0, n_distinct = 0;
				for (Int_t a = 0; a < 3; a++)
				{
					if (hfKeys[a] < 0) continue;
					n_hf++;
					Bool_t seen = false;
					for (Int_t b = 0; b < a; b++) if (hfKeys[b] == hfKeys[a]) seen = true;
					if (!seen) n_distinct++;
				}
				ftree->cand_gen_n_distinct_hf = n_distinct;

				// The two K slots sharing one immediate mother is what makes a real phi.
				ftree->cand_gen_KK_same_mom = (tKp.momKey >= 0) && (tKp.momKey == tKm.momKey);
				ftree->cand_gen_KK_mom_pdgId = ftree->cand_gen_KK_same_mom ? tKp.momPdgId : null;

				const Bool_t one_origin = (n_hf == 3) && (n_distinct == 1);
				ftree->cand_gen_common_hf_pdgId = one_origin ? tKp.hfPdgId : null;
				ftree->cand_gen_all_direct = one_origin && tKp.hfIsMom && tKm.hfIsMom && tpi.hfIsMom;

				// Convenience enum; every branch of it is re-derivable offline from the
				// raw fields above, which is why those are the primary product.
				const Int_t absHf = one_origin ? std::abs(tKp.hfPdgId) : 0;
				const Bool_t isDs = (absHf == 431);
				const Bool_t isD  = (absHf == 411) || (absHf == 421);
				// phi topology: both K from one phi, and the pi a direct daughter of the charm
				const Bool_t viaPhi = one_origin && ftree->cand_gen_KK_same_mom
					&& (ftree->cand_gen_KK_mom_pdgId == 333) && tpi.hfIsMom;
				// The candidate loop assigns the K,K,pi mass hypotheses blind, so these two
				// counters are the misID bookkeeping: 2K+1pi is the intended assignment,
				// 1K+2pi is one track taken for a kaon, 0K+3pi is two.
				const Int_t nTrueK = (std::abs(tKp.pdgId) == 321 ? 1 : 0)
					+ (std::abs(tKm.pdgId) == 321 ? 1 : 0) + (std::abs(tpi.pdgId) == 321 ? 1 : 0);
				const Int_t nTruePi = (std::abs(tKp.pdgId) == 211 ? 1 : 0)
					+ (std::abs(tKm.pdgId) == 211 ? 1 : 0) + (std::abs(tpi.pdgId) == 211 ? 1 : 0);
				ftree->cand_gen_ntrue_K = nTrueK;
				ftree->cand_gen_ntrue_pi = nTruePi;

				const Bool_t direct = ftree->cand_gen_all_direct;
				if (n_hf == 0)                       ftree->cand_gen_mode = 0;
				else if (n_distinct >= 2)            ftree->cand_gen_mode = 7;
				else if (!one_origin)                ftree->cand_gen_mode = 8;
				else if (isDs && viaPhi)             ftree->cand_gen_mode = 1;
				else if (isDs && direct && nTrueK == 2 && nTruePi == 1) ftree->cand_gen_mode = 2;
				else if (isD  && viaPhi)             ftree->cand_gen_mode = 3;
				else if (isD  && direct && nTrueK == 2 && nTruePi == 1) ftree->cand_gen_mode = 4;
				else if (isD  && direct && nTrueK == 1 && nTruePi == 2) ftree->cand_gen_mode = 5;
				// Ds-side twin of mode 5: Ds -> K pi pi with one pi taken for a kaon.
				else if (isDs && direct && nTrueK == 1 && nTruePi == 2) ftree->cand_gen_mode = 9;
				// Both K slots misID'd. Measured to dominate what used to fall into mode 6.
				else if ((isDs || isD) && direct && nTrueK == 0 && nTruePi == 3)
				                                     ftree->cand_gen_mode = 10;
				else                                 ftree->cand_gen_mode = 6;

				//--------------------------------------------------//
				// Gamma loop: one row per (Ds, gamma) pair
				// Cuts: dM in [50,250] MeV, Dsstmu_mass<125, Dsst_IsoR04<1
				//--------------------------------------------------//
				bool any_gamma = false;
				for (Int_t igamma = 0; igamma < static_cast<Int_t>(idx_photon.size()); igamma++)
				{
					const Int_t idx_g = idx_photon[igamma];
					const pat::PackedCandidate &gamma_pf = (*PFCandHandle)[idx_g];
					TLorentzVector gamma_tlv;
					gamma_tlv.SetXYZM(gamma_pf.px(), gamma_pf.py(), gamma_pf.pz(), 0.0);
					TLorentzVector Dsst_tlv   = Ds_tlv + gamma_tlv;
					if (std::abs(Dsst_tlv.M() - Mass_Dsst) > 0.15f) continue;
					Float_t dM = Dsst_tlv.M() - Ds_tlv.M();
					if (dM < 0.05f || dM > 0.25f) continue;

					TLorentzVector Dsstmu_tlv = Dsst_tlv + mu_tlv;
					if (Dsstmu_tlv.M() > 125.0f) continue;

					// Dsst isolation (exclude Kp, Km, pi, gamma)
					Float_t Dsst_IsoR03_sumChHadPt = 0, Dsst_IsoR03_sumNeuHadEt = 0, Dsst_IsoR03_sumPhoEt = 0, Dsst_IsoR03_sumPUPt = 0;
					Float_t Dsst_IsoR04_sumChHadPt = 0, Dsst_IsoR04_sumNeuHadEt = 0, Dsst_IsoR04_sumPhoEt = 0, Dsst_IsoR04_sumPUPt = 0;
					for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
					{
						if (i == idx_Kp || i == idx_Km || i == idx_pi || i == idx_g) continue;
						const pat::PackedCandidate &pf = (*PFCandHandle)[i];
						Float_t iso_dR = reco::deltaR(Dsst_tlv.Eta(), Dsst_tlv.Phi(), pf.eta(), pf.phi());
						if (iso_dR < 0.4)
						{
							if (std::abs(pf.pdgId()) == 211 && pf.fromPV() > 1) Dsst_IsoR04_sumChHadPt  += pf.pt();
							else if (pf.pdgId() == 130)                          Dsst_IsoR04_sumNeuHadEt += pf.et();
							else if (pf.pdgId() == 22)                           Dsst_IsoR04_sumPhoEt    += pf.et();
							if (pf.fromPV() <= 1) Dsst_IsoR04_sumPUPt += pf.pt();
							if (iso_dR < 0.3)
							{
								if (std::abs(pf.pdgId()) == 211 && pf.fromPV() > 1) Dsst_IsoR03_sumChHadPt  += pf.pt();
								else if (pf.pdgId() == 130)                          Dsst_IsoR03_sumNeuHadEt += pf.et();
								else if (pf.pdgId() == 22)                           Dsst_IsoR03_sumPhoEt    += pf.et();
								if (pf.fromPV() <= 1) Dsst_IsoR03_sumPUPt += pf.pt();
							}
						}
					}
					Float_t Dsst_IsoR03_PFIso = Dsst_IsoR03_sumChHadPt + std::max(0.0f, Dsst_IsoR03_sumNeuHadEt + Dsst_IsoR03_sumPhoEt - 0.5f * Dsst_IsoR03_sumPUPt);
					Float_t Dsst_IsoR03_RelIso = Dsst_IsoR03_PFIso / Dsst_tlv.Pt();
					Float_t Dsst_IsoR04_PFIso = Dsst_IsoR04_sumChHadPt + std::max(0.0f, Dsst_IsoR04_sumNeuHadEt + Dsst_IsoR04_sumPhoEt - 0.5f * Dsst_IsoR04_sumPUPt);
					Float_t Dsst_IsoR04_RelIso = Dsst_IsoR04_PFIso / Dsst_tlv.Pt();
					if (Dsst_IsoR04_RelIso > 1.0f) continue;

					// Fill gamma/Dsst/Dsstmu branches
					ftree->has_gamma      = true;
					ftree->gamma_charge   = gamma_pf.charge();
					ftree->gamma_eta      = gamma_tlv.Eta();  ftree->gamma_phi  = gamma_tlv.Phi();
					ftree->gamma_px       = gamma_tlv.Px();   ftree->gamma_py   = gamma_tlv.Py();   ftree->gamma_pz   = gamma_tlv.Pz();
					ftree->gamma_pt       = gamma_tlv.Pt();   ftree->gamma_p    = gamma_tlv.P();    ftree->gamma_mass = gamma_tlv.M();
					ftree->gamma_PVassoc  = gamma_pf.pvAssociationQuality();

					ftree->Dsst_charge = Kp_charge + Km_charge + pi_charge;
					ftree->Dsst_eta    = Dsst_tlv.Eta();  ftree->Dsst_phi  = Dsst_tlv.Phi();
					ftree->Dsst_px     = Dsst_tlv.Px();   ftree->Dsst_py   = Dsst_tlv.Py();   ftree->Dsst_pz   = Dsst_tlv.Pz();
					ftree->Dsst_pt     = Dsst_tlv.Pt();   ftree->Dsst_p    = Dsst_tlv.P();    ftree->Dsst_mass = Dsst_tlv.M();

					ftree->Dsstmu_charge = Kp_charge + Km_charge + pi_charge + mu.charge();
					ftree->Dsstmu_eta    = Dsstmu_tlv.Eta();  ftree->Dsstmu_phi  = Dsstmu_tlv.Phi();
					ftree->Dsstmu_px     = Dsstmu_tlv.Px();   ftree->Dsstmu_py   = Dsstmu_tlv.Py();   ftree->Dsstmu_pz   = Dsstmu_tlv.Pz();
					ftree->Dsstmu_pt     = Dsstmu_tlv.Pt();   ftree->Dsstmu_p    = Dsstmu_tlv.P();    ftree->Dsstmu_mass = Dsstmu_tlv.M();
					ftree->Dsstmu_Ht     = Kp_mtm.perp() + Km_mtm.perp() + pi_mtm.perp() + gamma_tlv.Pt() + mu_tlv.Pt();

					ftree->dR_Kp_gamma   = reco::deltaR(Kp_mtm.eta(), Kp_mtm.phi(), gamma_pf.eta(), gamma_pf.phi());
					ftree->dR_Kp_Dsst    = reco::deltaR(Kp_mtm.eta(), Kp_mtm.phi(), Dsst_tlv.Eta(), Dsst_tlv.Phi());
					ftree->dR_Km_gamma   = reco::deltaR(Km_mtm.eta(), Km_mtm.phi(), gamma_pf.eta(), gamma_pf.phi());
					ftree->dR_Km_Dsst    = reco::deltaR(Km_mtm.eta(), Km_mtm.phi(), Dsst_tlv.Eta(), Dsst_tlv.Phi());
					ftree->dR_pi_gamma   = reco::deltaR(pi_mtm.eta(), pi_mtm.phi(), gamma_pf.eta(), gamma_pf.phi());
					ftree->dR_pi_Dsst    = reco::deltaR(pi_mtm.eta(), pi_mtm.phi(), Dsst_tlv.Eta(), Dsst_tlv.Phi());
					ftree->dR_phi_gamma  = reco::deltaR(phi_tlv.Eta(), phi_tlv.Phi(), gamma_pf.eta(), gamma_pf.phi());
					ftree->dR_phi_Dsst   = reco::deltaR(phi_tlv.Eta(), phi_tlv.Phi(), Dsst_tlv.Eta(), Dsst_tlv.Phi());
					ftree->dR_Ds_gamma   = reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), gamma_pf.eta(), gamma_pf.phi());
					ftree->dR_Ds_Dsst    = reco::deltaR(Ds_mtm.eta(), Ds_mtm.phi(), Dsst_tlv.Eta(), Dsst_tlv.Phi());
					ftree->dR_gamma_Dsst = reco::deltaR(gamma_pf.eta(), gamma_pf.phi(), Dsst_tlv.Eta(), Dsst_tlv.Phi());
					ftree->dR_gamma_mu   = reco::deltaR(gamma_pf.eta(), gamma_pf.phi(), mu.eta(), mu.phi());
					ftree->dR_Dsst_mu    = reco::deltaR(Dsst_tlv.Eta(), Dsst_tlv.Phi(), mu.eta(), mu.phi());

					ftree->Dsst_IsoR03_sumChargedHadronPt = Dsst_IsoR03_sumChHadPt;
					ftree->Dsst_IsoR03_sumNeutralHadronEt = Dsst_IsoR03_sumNeuHadEt;
					ftree->Dsst_IsoR03_sumPhotonEt        = Dsst_IsoR03_sumPhoEt;
					ftree->Dsst_IsoR03_sumPUPt            = Dsst_IsoR03_sumPUPt;
					ftree->Dsst_IsoR03_PFIso              = Dsst_IsoR03_PFIso;
					ftree->Dsst_IsoR03_RelIso             = Dsst_IsoR03_RelIso;
					ftree->Dsst_IsoR04_sumChargedHadronPt = Dsst_IsoR04_sumChHadPt;
					ftree->Dsst_IsoR04_sumNeutralHadronEt = Dsst_IsoR04_sumNeuHadEt;
					ftree->Dsst_IsoR04_sumPhotonEt        = Dsst_IsoR04_sumPhoEt;
					ftree->Dsst_IsoR04_sumPUPt            = Dsst_IsoR04_sumPUPt;
					ftree->Dsst_IsoR04_PFIso              = Dsst_IsoR04_PFIso;
					ftree->Dsst_IsoR04_RelIso             = Dsst_IsoR04_RelIso;

					{ pat::CompositeCandidate c; fillCand(c); out->push_back(c); }
					ftree->InitGamma();
					any_gamma = true;
				}

				// No passing gamma: one row with has_gamma=false
				if (!any_gamma)
				{
					{ pat::CompositeCandidate c; fillCand(c); out->push_back(c); }
					ftree->InitGamma();
				}
			} // pi loop
		} // Km loop
	} // Kp loop

	const int nCand = static_cast<int>(out->size());
	for (auto &c : *out) c.addUserInt("nCandidates", nCand);

	iEvent.put(std::move(out));
}

void HDsPhiPiMuGenTagProducer::fillCand(pat::CompositeCandidate &c) const
{
	c.addUserInt("cand_id", ftree->cand_id);
	c.addUserInt("n_gen_Ds", ftree->n_gen_Ds);
	c.addUserInt("n_gen_D", ftree->n_gen_D);
	c.addUserInt("n_gen_Ds_phipi", ftree->n_gen_Ds_phipi);
	c.addUserInt("trig_HLT_IsoMu24", ftree->trig_HLT_IsoMu24);
	c.addUserInt("trig_HLT_Mu50", ftree->trig_HLT_Mu50);
	c.addUserInt("primvtx_isValid", ftree->primvtx_isValid);
	c.addUserInt("primvtx_isFake", ftree->primvtx_isFake);
	c.addUserFloat("primvtx_chi2", ftree->primvtx_chi2);
	c.addUserInt("primvtx_ndof", ftree->primvtx_ndof);
	c.addUserFloat("primvtx_chi2ndof", ftree->primvtx_chi2ndof);
	c.addUserFloat("primvtx_vx", ftree->primvtx_vx);
	c.addUserFloat("primvtx_vy", ftree->primvtx_vy);
	c.addUserFloat("primvtx_vz", ftree->primvtx_vz);
	c.addUserFloat("primvtx_vxerr", ftree->primvtx_vxerr);
	c.addUserFloat("primvtx_vyerr", ftree->primvtx_vyerr);
	c.addUserFloat("primvtx_vzerr", ftree->primvtx_vzerr);
	c.addUserFloat("met_sig", ftree->met_sig);
	c.addUserFloat("met_phi", ftree->met_phi);
	c.addUserFloat("met_px", ftree->met_px);
	c.addUserFloat("met_py", ftree->met_py);
	c.addUserFloat("met_pt", ftree->met_pt);
	c.addUserFloat("met_sumEt", ftree->met_sumEt);
	c.addUserFloat("W_tmass", ftree->W_tmass);
	c.addUserInt("num_jets", ftree->num_jets);
	c.addUserInt("num_jets_pt20", ftree->num_jets_pt20);
	c.addUserInt("num_jets_pt30", ftree->num_jets_pt30);
	c.addUserInt("num_jets_pt40", ftree->num_jets_pt40);
	c.addUserInt("num_jets_pt50", ftree->num_jets_pt50);
	c.addUserInt("num_jets_pt60", ftree->num_jets_pt60);
	c.addUserInt("num_jets_pt70", ftree->num_jets_pt70);
	c.addUserInt("num_jets_pt80", ftree->num_jets_pt80);
	c.addUserInt("num_jets_pt90", ftree->num_jets_pt90);
	c.addUserInt("num_jets_pt100", ftree->num_jets_pt100);
	c.addUserFloat("jet_ht", ftree->jet_ht);
	c.addUserFloat("jet_ht_pt20", ftree->jet_ht_pt20);
	c.addUserFloat("jet_leadpt", ftree->jet_leadpt);
	c.addUserFloat("jet_maxbscore", ftree->jet_maxbscore);
	c.addUserFloat("nTrueInt", ftree->nTrueInt);
	c.addUserInt("nPUInt", ftree->nPUInt);
	c.addUserInt("num_os_mu_loose", ftree->num_os_mu_loose);
	c.addUserInt("num_ss_mu_loose", ftree->num_ss_mu_loose);
	c.addUserInt("num_os_mu_medium", ftree->num_os_mu_medium);
	c.addUserInt("num_ss_mu_medium", ftree->num_ss_mu_medium);
	c.addUserInt("num_os_mu_tight", ftree->num_os_mu_tight);
	c.addUserInt("num_ss_mu_tight", ftree->num_ss_mu_tight);
	c.addUserInt("num_os_e_cutloose", ftree->num_os_e_cutloose);
	c.addUserInt("num_ss_e_cutloose", ftree->num_ss_e_cutloose);
	c.addUserInt("num_os_e_cutmedium", ftree->num_os_e_cutmedium);
	c.addUserInt("num_ss_e_cutmedium", ftree->num_ss_e_cutmedium);
	c.addUserInt("num_os_e_cuttight", ftree->num_os_e_cuttight);
	c.addUserInt("num_ss_e_cuttight", ftree->num_ss_e_cuttight);
	c.addUserInt("num_os_e_mvawp80", ftree->num_os_e_mvawp80);
	c.addUserInt("num_ss_e_mvawp80", ftree->num_ss_e_mvawp80);
	c.addUserInt("num_os_e_mvawp90", ftree->num_os_e_mvawp90);
	c.addUserInt("num_ss_e_mvawp90", ftree->num_ss_e_mvawp90);
	c.addUserInt("mu_charge", ftree->mu_charge);
	c.addUserFloat("mu_eta", ftree->mu_eta);
	c.addUserFloat("mu_phi", ftree->mu_phi);
	c.addUserFloat("mu_px", ftree->mu_px);
	c.addUserFloat("mu_py", ftree->mu_py);
	c.addUserFloat("mu_pz", ftree->mu_pz);
	c.addUserFloat("mu_pt", ftree->mu_pt);
	c.addUserFloat("mu_p", ftree->mu_p);
	c.addUserInt("mu_isHighPt", ftree->mu_isHighPt);
	c.addUserInt("mu_isLoose", ftree->mu_isLoose);
	c.addUserInt("mu_isMedium", ftree->mu_isMedium);
	c.addUserInt("mu_isSoft", ftree->mu_isSoft);
	c.addUserInt("mu_isTight", ftree->mu_isTight);
	c.addUserInt("mu_isPF", ftree->mu_isPF);
	c.addUserInt("mu_isTracker", ftree->mu_isTracker);
	c.addUserInt("mu_isGlobal", ftree->mu_isGlobal);
	c.addUserFloat("mu_ip", ftree->mu_ip);
	c.addUserFloat("mu_iperr", ftree->mu_iperr);
	c.addUserFloat("mu_ipsig", ftree->mu_ipsig);
	c.addUserFloat("mu_ipxy", ftree->mu_ipxy);
	c.addUserFloat("mu_ipxyerr", ftree->mu_ipxyerr);
	c.addUserFloat("mu_ipxysig", ftree->mu_ipxysig);
	c.addUserFloat("mu_IsoR03_sumChargedHadronPt", ftree->mu_IsoR03_sumChargedHadronPt);
	c.addUserFloat("mu_IsoR03_sumNeutralHadronEt", ftree->mu_IsoR03_sumNeutralHadronEt);
	c.addUserFloat("mu_IsoR03_sumPhotonEt", ftree->mu_IsoR03_sumPhotonEt);
	c.addUserFloat("mu_IsoR03_sumPUPt", ftree->mu_IsoR03_sumPUPt);
	c.addUserFloat("mu_IsoR03_PFIso", ftree->mu_IsoR03_PFIso);
	c.addUserFloat("mu_IsoR03_RelIso", ftree->mu_IsoR03_RelIso);
	c.addUserFloat("mu_IsoR04_sumChargedHadronPt", ftree->mu_IsoR04_sumChargedHadronPt);
	c.addUserFloat("mu_IsoR04_sumNeutralHadronEt", ftree->mu_IsoR04_sumNeutralHadronEt);
	c.addUserFloat("mu_IsoR04_sumPhotonEt", ftree->mu_IsoR04_sumPhotonEt);
	c.addUserFloat("mu_IsoR04_sumPUPt", ftree->mu_IsoR04_sumPUPt);
	c.addUserFloat("mu_IsoR04_PFIso", ftree->mu_IsoR04_PFIso);
	c.addUserFloat("mu_IsoR04_RelIso", ftree->mu_IsoR04_RelIso);
	c.addUserInt("Kp_charge", ftree->Kp_charge);
	c.addUserFloat("Kp_eta", ftree->Kp_eta);
	c.addUserFloat("Kp_phi", ftree->Kp_phi);
	c.addUserFloat("Kp_px", ftree->Kp_px);
	c.addUserFloat("Kp_py", ftree->Kp_py);
	c.addUserFloat("Kp_pz", ftree->Kp_pz);
	c.addUserFloat("Kp_pt", ftree->Kp_pt);
	c.addUserFloat("Kp_p", ftree->Kp_p);
	c.addUserFloat("Kp_mass", ftree->Kp_mass);
	c.addUserInt("Kp_PVassoc", ftree->Kp_PVassoc);
	c.addUserInt("Km_charge", ftree->Km_charge);
	c.addUserFloat("Km_eta", ftree->Km_eta);
	c.addUserFloat("Km_phi", ftree->Km_phi);
	c.addUserFloat("Km_px", ftree->Km_px);
	c.addUserFloat("Km_py", ftree->Km_py);
	c.addUserFloat("Km_pz", ftree->Km_pz);
	c.addUserFloat("Km_pt", ftree->Km_pt);
	c.addUserFloat("Km_p", ftree->Km_p);
	c.addUserFloat("Km_mass", ftree->Km_mass);
	c.addUserInt("Km_PVassoc", ftree->Km_PVassoc);
	c.addUserInt("pi_charge", ftree->pi_charge);
	c.addUserFloat("pi_eta", ftree->pi_eta);
	c.addUserFloat("pi_phi", ftree->pi_phi);
	c.addUserFloat("pi_px", ftree->pi_px);
	c.addUserFloat("pi_py", ftree->pi_py);
	c.addUserFloat("pi_pz", ftree->pi_pz);
	c.addUserFloat("pi_pt", ftree->pi_pt);
	c.addUserFloat("pi_p", ftree->pi_p);
	c.addUserFloat("pi_mass", ftree->pi_mass);
	c.addUserInt("pi_PVassoc", ftree->pi_PVassoc);
	c.addUserInt("phi_charge", ftree->phi_charge);
	c.addUserFloat("phi_eta", ftree->phi_eta);
	c.addUserFloat("phi_phi", ftree->phi_phi);
	c.addUserFloat("phi_px", ftree->phi_px);
	c.addUserFloat("phi_py", ftree->phi_py);
	c.addUserFloat("phi_pz", ftree->phi_pz);
	c.addUserFloat("phi_pt", ftree->phi_pt);
	c.addUserFloat("phi_p", ftree->phi_p);
	c.addUserFloat("phi_mass", ftree->phi_mass);
	c.addUserInt("Ds_charge", ftree->Ds_charge);
	c.addUserFloat("Ds_eta", ftree->Ds_eta);
	c.addUserFloat("Ds_phi", ftree->Ds_phi);
	c.addUserFloat("Ds_px", ftree->Ds_px);
	c.addUserFloat("Ds_py", ftree->Ds_py);
	c.addUserFloat("Ds_pz", ftree->Ds_pz);
	c.addUserFloat("Ds_pt", ftree->Ds_pt);
	c.addUserFloat("Ds_p", ftree->Ds_p);
	c.addUserFloat("Ds_mass", ftree->Ds_mass);
	c.addUserInt("Dsmu_charge", ftree->Dsmu_charge);
	c.addUserFloat("Dsmu_eta", ftree->Dsmu_eta);
	c.addUserFloat("Dsmu_phi", ftree->Dsmu_phi);
	c.addUserFloat("Dsmu_px", ftree->Dsmu_px);
	c.addUserFloat("Dsmu_py", ftree->Dsmu_py);
	c.addUserFloat("Dsmu_pz", ftree->Dsmu_pz);
	c.addUserFloat("Dsmu_pt", ftree->Dsmu_pt);
	c.addUserFloat("Dsmu_p", ftree->Dsmu_p);
	c.addUserFloat("Dsmu_mass", ftree->Dsmu_mass);
	c.addUserFloat("Dsmu_Ht", ftree->Dsmu_Ht);
	c.addUserInt("n_gamma_window", ftree->n_gamma_window);
	c.addUserFloat("gamma_uc_dM", ftree->gamma_uc_dM);
	c.addUserFloat("gamma_uc_pt", ftree->gamma_uc_pt);
	c.addUserFloat("gamma_uc_dR_Ds", ftree->gamma_uc_dR_Ds);
	c.addUserInt("has_gamma", ftree->has_gamma);
	c.addUserInt("gamma_charge", ftree->gamma_charge);
	c.addUserFloat("gamma_eta", ftree->gamma_eta);
	c.addUserFloat("gamma_phi", ftree->gamma_phi);
	c.addUserFloat("gamma_px", ftree->gamma_px);
	c.addUserFloat("gamma_py", ftree->gamma_py);
	c.addUserFloat("gamma_pz", ftree->gamma_pz);
	c.addUserFloat("gamma_pt", ftree->gamma_pt);
	c.addUserFloat("gamma_p", ftree->gamma_p);
	c.addUserFloat("gamma_mass", ftree->gamma_mass);
	c.addUserInt("gamma_PVassoc", ftree->gamma_PVassoc);
	c.addUserInt("Dsst_charge", ftree->Dsst_charge);
	c.addUserFloat("Dsst_eta", ftree->Dsst_eta);
	c.addUserFloat("Dsst_phi", ftree->Dsst_phi);
	c.addUserFloat("Dsst_px", ftree->Dsst_px);
	c.addUserFloat("Dsst_py", ftree->Dsst_py);
	c.addUserFloat("Dsst_pz", ftree->Dsst_pz);
	c.addUserFloat("Dsst_pt", ftree->Dsst_pt);
	c.addUserFloat("Dsst_p", ftree->Dsst_p);
	c.addUserFloat("Dsst_mass", ftree->Dsst_mass);
	c.addUserInt("Dsstmu_charge", ftree->Dsstmu_charge);
	c.addUserFloat("Dsstmu_eta", ftree->Dsstmu_eta);
	c.addUserFloat("Dsstmu_phi", ftree->Dsstmu_phi);
	c.addUserFloat("Dsstmu_px", ftree->Dsstmu_px);
	c.addUserFloat("Dsstmu_py", ftree->Dsstmu_py);
	c.addUserFloat("Dsstmu_pz", ftree->Dsstmu_pz);
	c.addUserFloat("Dsstmu_pt", ftree->Dsstmu_pt);
	c.addUserFloat("Dsstmu_p", ftree->Dsstmu_p);
	c.addUserFloat("Dsstmu_mass", ftree->Dsstmu_mass);
	c.addUserFloat("Dsstmu_Ht", ftree->Dsstmu_Ht);
	c.addUserFloat("Ds_vx", ftree->Ds_vx);
	c.addUserFloat("Ds_vy", ftree->Ds_vy);
	c.addUserFloat("Ds_vz", ftree->Ds_vz);
	c.addUserFloat("Ds_vxerr", ftree->Ds_vxerr);
	c.addUserFloat("Ds_vyerr", ftree->Ds_vyerr);
	c.addUserFloat("Ds_vzerr", ftree->Ds_vzerr);
	c.addUserFloat("Ds_chi2", ftree->Ds_chi2);
	c.addUserInt("Ds_ndf", ftree->Ds_ndf);
	c.addUserFloat("Ds_chi2ndf", ftree->Ds_chi2ndf);
	c.addUserFloat("Ds_prob", ftree->Ds_prob);
	c.addUserFloat("dR_Kp_Km", ftree->dR_Kp_Km);
	c.addUserFloat("dR_Kp_pi", ftree->dR_Kp_pi);
	c.addUserFloat("dR_Kp_Ds", ftree->dR_Kp_Ds);
	c.addUserFloat("dR_Km_pi", ftree->dR_Km_pi);
	c.addUserFloat("dR_Km_Ds", ftree->dR_Km_Ds);
	c.addUserFloat("dR_pi_Ds", ftree->dR_pi_Ds);
	c.addUserFloat("dR_phi_pi", ftree->dR_phi_pi);
	c.addUserFloat("dR_phi_Ds", ftree->dR_phi_Ds);
	c.addUserFloat("dR_Kp_mu", ftree->dR_Kp_mu);
	c.addUserFloat("dR_Km_mu", ftree->dR_Km_mu);
	c.addUserFloat("dR_pi_mu", ftree->dR_pi_mu);
	c.addUserFloat("dR_phi_mu", ftree->dR_phi_mu);
	c.addUserFloat("dR_Ds_mu", ftree->dR_Ds_mu);
	c.addUserFloat("dR_Kp_gamma", ftree->dR_Kp_gamma);
	c.addUserFloat("dR_Kp_Dsst", ftree->dR_Kp_Dsst);
	c.addUserFloat("dR_Km_gamma", ftree->dR_Km_gamma);
	c.addUserFloat("dR_Km_Dsst", ftree->dR_Km_Dsst);
	c.addUserFloat("dR_pi_gamma", ftree->dR_pi_gamma);
	c.addUserFloat("dR_pi_Dsst", ftree->dR_pi_Dsst);
	c.addUserFloat("dR_phi_gamma", ftree->dR_phi_gamma);
	c.addUserFloat("dR_phi_Dsst", ftree->dR_phi_Dsst);
	c.addUserFloat("dR_Ds_gamma", ftree->dR_Ds_gamma);
	c.addUserFloat("dR_Ds_Dsst", ftree->dR_Ds_Dsst);
	c.addUserFloat("dR_gamma_Dsst", ftree->dR_gamma_Dsst);
	c.addUserFloat("dR_gamma_mu", ftree->dR_gamma_mu);
	c.addUserFloat("dR_Dsst_mu", ftree->dR_Dsst_mu);
	c.addUserFloat("Kp_ip", ftree->Kp_ip);
	c.addUserFloat("Kp_iperr", ftree->Kp_iperr);
	c.addUserFloat("Kp_ipsig", ftree->Kp_ipsig);
	c.addUserFloat("Kp_ipxy", ftree->Kp_ipxy);
	c.addUserFloat("Kp_ipxyerr", ftree->Kp_ipxyerr);
	c.addUserFloat("Kp_ipxysig", ftree->Kp_ipxysig);
	c.addUserFloat("Km_ip", ftree->Km_ip);
	c.addUserFloat("Km_iperr", ftree->Km_iperr);
	c.addUserFloat("Km_ipsig", ftree->Km_ipsig);
	c.addUserFloat("Km_ipxy", ftree->Km_ipxy);
	c.addUserFloat("Km_ipxyerr", ftree->Km_ipxyerr);
	c.addUserFloat("Km_ipxysig", ftree->Km_ipxysig);
	c.addUserFloat("pi_ip", ftree->pi_ip);
	c.addUserFloat("pi_iperr", ftree->pi_iperr);
	c.addUserFloat("pi_ipsig", ftree->pi_ipsig);
	c.addUserFloat("pi_ipxy", ftree->pi_ipxy);
	c.addUserFloat("pi_ipxyerr", ftree->pi_ipxyerr);
	c.addUserFloat("pi_ipxysig", ftree->pi_ipxysig);
	c.addUserFloat("Ds_ip", ftree->Ds_ip);
	c.addUserFloat("Ds_iperr", ftree->Ds_iperr);
	c.addUserFloat("Ds_ipsig", ftree->Ds_ipsig);
	c.addUserFloat("Ds_ipxy", ftree->Ds_ipxy);
	c.addUserFloat("Ds_ipxyerr", ftree->Ds_ipxyerr);
	c.addUserFloat("Ds_ipxysig", ftree->Ds_ipxysig);
	c.addUserFloat("dira", ftree->dira);
	c.addUserFloat("Ds_FDxy", ftree->Ds_FDxy);
	c.addUserFloat("Ds_FDxyerr", ftree->Ds_FDxyerr);
	c.addUserFloat("Ds_FDxysig", ftree->Ds_FDxysig);
	c.addUserFloat("Ds_FD", ftree->Ds_FD);
	c.addUserFloat("Ds_FDerr", ftree->Ds_FDerr);
	c.addUserFloat("Ds_FDsig", ftree->Ds_FDsig);
	c.addUserFloat("Ds_IsoR03_sumChargedHadronPt", ftree->Ds_IsoR03_sumChargedHadronPt);
	c.addUserFloat("Ds_IsoR03_sumNeutralHadronEt", ftree->Ds_IsoR03_sumNeutralHadronEt);
	c.addUserFloat("Ds_IsoR03_sumPhotonEt", ftree->Ds_IsoR03_sumPhotonEt);
	c.addUserFloat("Ds_IsoR03_sumPUPt", ftree->Ds_IsoR03_sumPUPt);
	c.addUserFloat("Ds_IsoR03_PFIso", ftree->Ds_IsoR03_PFIso);
	c.addUserFloat("Ds_IsoR03_RelIso", ftree->Ds_IsoR03_RelIso);
	c.addUserFloat("Ds_IsoR04_sumChargedHadronPt", ftree->Ds_IsoR04_sumChargedHadronPt);
	c.addUserFloat("Ds_IsoR04_sumNeutralHadronEt", ftree->Ds_IsoR04_sumNeutralHadronEt);
	c.addUserFloat("Ds_IsoR04_sumPhotonEt", ftree->Ds_IsoR04_sumPhotonEt);
	c.addUserFloat("Ds_IsoR04_sumPUPt", ftree->Ds_IsoR04_sumPUPt);
	c.addUserFloat("Ds_IsoR04_PFIso", ftree->Ds_IsoR04_PFIso);
	c.addUserFloat("Ds_IsoR04_RelIso", ftree->Ds_IsoR04_RelIso);
	c.addUserFloat("Ds_nearestjet_dR", ftree->Ds_nearestjet_dR);
	c.addUserFloat("Ds_nearestjet_bscore", ftree->Ds_nearestjet_bscore);
	c.addUserInt("n_b_jets_loose_dR03", ftree->n_b_jets_loose_dR03);
	c.addUserInt("n_b_jets_loose_dR04", ftree->n_b_jets_loose_dR04);
	c.addUserInt("n_b_jets_medium_dR03", ftree->n_b_jets_medium_dR03);
	c.addUserInt("n_b_jets_medium_dR04", ftree->n_b_jets_medium_dR04);
	c.addUserInt("n_b_jets_tight_dR03", ftree->n_b_jets_tight_dR03);
	c.addUserInt("n_b_jets_tight_dR04", ftree->n_b_jets_tight_dR04);
	c.addUserInt("n_b_jets_loose_pt20_dR03", ftree->n_b_jets_loose_pt20_dR03);
	c.addUserInt("n_b_jets_loose_pt20_dR04", ftree->n_b_jets_loose_pt20_dR04);
	c.addUserInt("n_b_jets_medium_pt20_dR03", ftree->n_b_jets_medium_pt20_dR03);
	c.addUserInt("n_b_jets_medium_pt20_dR04", ftree->n_b_jets_medium_pt20_dR04);
	c.addUserInt("n_b_jets_tight_pt20_dR03", ftree->n_b_jets_tight_pt20_dR03);
	c.addUserInt("n_b_jets_tight_pt20_dR04", ftree->n_b_jets_tight_pt20_dR04);
	c.addUserFloat("Dsst_IsoR03_sumChargedHadronPt", ftree->Dsst_IsoR03_sumChargedHadronPt);
	c.addUserFloat("Dsst_IsoR03_sumNeutralHadronEt", ftree->Dsst_IsoR03_sumNeutralHadronEt);
	c.addUserFloat("Dsst_IsoR03_sumPhotonEt", ftree->Dsst_IsoR03_sumPhotonEt);
	c.addUserFloat("Dsst_IsoR03_sumPUPt", ftree->Dsst_IsoR03_sumPUPt);
	c.addUserFloat("Dsst_IsoR03_PFIso", ftree->Dsst_IsoR03_PFIso);
	c.addUserFloat("Dsst_IsoR03_RelIso", ftree->Dsst_IsoR03_RelIso);
	c.addUserFloat("Dsst_IsoR04_sumChargedHadronPt", ftree->Dsst_IsoR04_sumChargedHadronPt);
	c.addUserFloat("Dsst_IsoR04_sumNeutralHadronEt", ftree->Dsst_IsoR04_sumNeutralHadronEt);
	c.addUserFloat("Dsst_IsoR04_sumPhotonEt", ftree->Dsst_IsoR04_sumPhotonEt);
	c.addUserFloat("Dsst_IsoR04_sumPUPt", ftree->Dsst_IsoR04_sumPUPt);
	c.addUserFloat("Dsst_IsoR04_PFIso", ftree->Dsst_IsoR04_PFIso);
	c.addUserFloat("Dsst_IsoR04_RelIso", ftree->Dsst_IsoR04_RelIso);
	c.addUserInt("Kp_gen_pdgId", ftree->Kp_gen_pdgId);
	c.addUserInt("Km_gen_pdgId", ftree->Km_gen_pdgId);
	c.addUserInt("pi_gen_pdgId", ftree->pi_gen_pdgId);
	c.addUserInt("Kp_gen_mom_pdgId", ftree->Kp_gen_mom_pdgId);
	c.addUserInt("Km_gen_mom_pdgId", ftree->Km_gen_mom_pdgId);
	c.addUserInt("pi_gen_mom_pdgId", ftree->pi_gen_mom_pdgId);
	c.addUserInt("Kp_gen_mom_key", ftree->Kp_gen_mom_key);
	c.addUserInt("Km_gen_mom_key", ftree->Km_gen_mom_key);
	c.addUserInt("pi_gen_mom_key", ftree->pi_gen_mom_key);
	c.addUserInt("Kp_gen_hf_pdgId", ftree->Kp_gen_hf_pdgId);
	c.addUserInt("Km_gen_hf_pdgId", ftree->Km_gen_hf_pdgId);
	c.addUserInt("pi_gen_hf_pdgId", ftree->pi_gen_hf_pdgId);
	c.addUserInt("Kp_gen_hf_key", ftree->Kp_gen_hf_key);
	c.addUserInt("Km_gen_hf_key", ftree->Km_gen_hf_key);
	c.addUserInt("pi_gen_hf_key", ftree->pi_gen_hf_key);
	c.addUserInt("Kp_gen_hf_isMom", ftree->Kp_gen_hf_isMom);
	c.addUserInt("Km_gen_hf_isMom", ftree->Km_gen_hf_isMom);
	c.addUserInt("pi_gen_hf_isMom", ftree->pi_gen_hf_isMom);
	c.addUserInt("cand_gen_nmatched", ftree->cand_gen_nmatched);
	c.addUserInt("cand_gen_n_distinct_hf", ftree->cand_gen_n_distinct_hf);
	c.addUserInt("cand_gen_KK_same_mom", ftree->cand_gen_KK_same_mom);
	c.addUserInt("cand_gen_KK_mom_pdgId", ftree->cand_gen_KK_mom_pdgId);
	c.addUserInt("cand_gen_common_hf_pdgId", ftree->cand_gen_common_hf_pdgId);
	c.addUserInt("cand_gen_all_direct", ftree->cand_gen_all_direct);
	c.addUserInt("cand_gen_ntrue_K", ftree->cand_gen_ntrue_K);
	c.addUserInt("cand_gen_ntrue_pi", ftree->cand_gen_ntrue_pi);
	c.addUserInt("cand_gen_mode", ftree->cand_gen_mode);

	c.setP4(math::PtEtaPhiMLorentzVector(ftree->Ds_pt, ftree->Ds_eta, ftree->Ds_phi, ftree->Ds_mass));
	c.setCharge(ftree->Ds_charge);
}

Bool_t HDsPhiPiMuGenTagProducer::isHeavyFlavour(const Int_t pdgId)
{
	const Int_t a = std::abs(pdgId);
	// Charm and bottom mesons (4xx, 5xx) and baryons (4xxx, 5xxx). Same ranges the MiniAOD
	// pruning uses to keep these hadrons unconditionally, which is what makes their
	// first-generation daughters resolvable through the association at all.
	return ((a > 400 && a < 600) || (a > 4000 && a < 6000));
}

const reco::GenParticle *HDsPhiPiMuGenTagProducer::nearestHeavyFlavour(const reco::GenParticle &gp, Int_t &key, Int_t depth) const
{
	// Depth guard: gen mother chains can loop through copy records, and this runs once per
	// pool track per event.
	if (depth > 20) return nullptr;
	for (size_t i = 0; i < gp.numberOfMothers(); i++)
	{
		const reco::GenParticleRef momRef = gp.motherRef(i);
		if (momRef.isNull()) continue;
		if (isHeavyFlavour(momRef->pdgId()))
		{
			key = static_cast<Int_t>(momRef.key());
			return &(*momRef);
		}
		if (const reco::GenParticle *up = nearestHeavyFlavour(*momRef, key, depth + 1)) return up;
	}
	return nullptr;
}

HDsPhiPiMuGenTagProducer::PFTruth HDsPhiPiMuGenTagProducer::pfTruth(const pat::PackedCandidate &pf,
                                                   const pat::PackedCandidateRef &pfRef,
                                                   const edm::Association<reco::GenParticleCollection> &assoc) const
{
	PFTruth t;
	const reco::GenParticleRef genRef = assoc[pfRef];
	// A null ref is the normal outcome for an ordinary fragmentation track: the association
	// re-keys into prunedGenParticles, whose "drop status == 1" removes every stable particle
	// that is not re-kept as a c/b-hadron or phi/K_S/Lambda daughter (or a lepton etc.).
	if (genRef.isNull()) return t;

	const reco::GenParticle &gp = *genRef;
	t.valid  = true;
	t.pdgId  = gp.pdgId();

	if (gp.numberOfMothers() > 0)
	{
		const reco::GenParticleRef momRef = gp.motherRef(0);
		if (momRef.isNonnull())
		{
			t.momPdgId = momRef->pdgId();
			t.momKey   = static_cast<Int_t>(momRef.key());
		}
	}

	Int_t hfKey = -1;
	if (const reco::GenParticle *hf = nearestHeavyFlavour(gp, hfKey))
	{
		t.hfPdgId = hf->pdgId();
		t.hfKey   = hfKey;
		// "Direct" means the charm/bottom hadron is the immediate mother, i.e. a genuine
		// three-body decay rather than one proceeding through a phi (or any other resonance).
		t.hfIsMom = (t.momKey >= 0) && (t.momKey == hfKey);
	}
	return t;
}

Float_t HDsPhiPiMuGenTagProducer::bScore(const pat::Jet &jet) const
{
	Float_t score = 0;
	for (const std::string &name : bDiscNames_) score += jet.bDiscriminator(name);
	return score;
}

Measurement1D HDsPhiPiMuGenTagProducer::distanceXY(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
{
	AlgebraicSymMatrix33 error = vtx1PositionError.matrix() + vtx2PositionError.matrix();
	GlobalVector diff = vtx1Position - vtx2Position;
	AlgebraicVector3 vDiff;
	vDiff[0] = diff.x(); vDiff[1] = diff.y(); vDiff[2] = 0.;
	Double_t dist = sqrt(pow(diff.x(), 2) + pow(diff.y(), 2));
	Double_t err2 = ROOT::Math::Similarity(error, vDiff);
	Double_t err = (dist != 0) ? sqrt(err2) / dist : 0;
	return Measurement1D(dist, err);
}

Measurement1D HDsPhiPiMuGenTagProducer::distance3D(const GlobalPoint &vtx1Position, const GlobalError &vtx1PositionError, const GlobalPoint &vtx2Position, const GlobalError &vtx2PositionError) const
{
	AlgebraicSymMatrix33 error = vtx1PositionError.matrix() + vtx2PositionError.matrix();
	GlobalVector diff = vtx1Position - vtx2Position;
	AlgebraicVector3 vDiff;
	vDiff[0] = diff.x(); vDiff[1] = diff.y(); vDiff[2] = diff.z();
	Double_t dist = diff.mag();
	Double_t err2 = ROOT::Math::Similarity(error, vDiff);
	Double_t err = (dist != 0) ? sqrt(err2) / dist : 0;
	return Measurement1D(dist, err);
}

void HDsPhiPiMuGenTagProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.add<bool>("sameSignPi", false);
	desc.add<double>("bScoreWPLoose");   // required, era-dependent (no default on purpose)
	desc.add<double>("bScoreWPMedium");
	desc.add<double>("bScoreWPTight");
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(HDsPhiPiMuGenTagProducer);
