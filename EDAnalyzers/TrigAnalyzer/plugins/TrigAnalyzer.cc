#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "EDAnalyzers/TrigAnalyzer/interface/RecoTree.h"

#include "TLorentzVector.h"

class TrigAnalyzer : public edm::one::EDAnalyzer<>
{
public:
	explicit TrigAnalyzer(const edm::ParameterSet &);
	~TrigAnalyzer() override {};

	static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
	void beginJob() override;
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	bool hasAncestor(const reco::GenParticle &gp, const int motherid, const int otherparticleid) const;
	void endJob() override;

	edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
	edm::EDGetTokenT<reco::VertexCollection> primvtxToken_;
	edm::EDGetTokenT<reco::GenParticleCollection> GenPartToken_;
	edm::EDGetTokenT<pat::PackedCandidateCollection> PFCandToken_;
	edm::EDGetTokenT<pat::MuonCollection> MuonToken_;
	edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> associationToken_;
	edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

	const Double_t Mass_Ds = 1.96835;
	const Double_t Mass_phi = 1.019460;
	const Double_t Mass_K = 0.493677;
	const Double_t Mass_pi = 0.13957039;
	const Double_t Mass_mu = 0.1056583755;
	unsigned long nEventsProcessed_;
	unsigned long nEventsValidKp_;
	unsigned long nEventsValidKm_;
	unsigned long nEventsValidpi_;
	unsigned long nEventsValidmu_;

	std::map<TString, int> triggerCounts_;
	std::map<TString, int> triggerTotal_;

	const edm::Service<TFileService> fs;
	RecoTree *ftree;

	KalmanVertexFitter KVFitter;
	// VertexReProducer *revertex;
};

TrigAnalyzer::TrigAnalyzer(const edm::ParameterSet &iConfig)
	: KVFitter(true)
{
	edm::InputTag TriggerBitsTag_("TriggerResults", "", "HLT");
	triggerToken_ = consumes<edm::TriggerResults>(TriggerBitsTag_);

	edm::InputTag primvtxTag_("offlineSlimmedPrimaryVertices");
	primvtxToken_ = consumes<reco::VertexCollection>(primvtxTag_);

	edm::InputTag GenPartTag_("prunedGenParticles");
	GenPartToken_ = consumes<reco::GenParticleCollection>(GenPartTag_);

	edm::InputTag PFCandTag_("packedPFCandidates");
	PFCandToken_ = consumes<pat::PackedCandidateCollection>(PFCandTag_);

	edm::InputTag MuonTag_("slimmedMuons");
	MuonToken_ = consumes<pat::MuonCollection>(MuonTag_);

	edm::InputTag assotiationTag_("packedPFCandidateToGenAssociation");
	associationToken_ = consumes<edm::Association<reco::GenParticleCollection>>(assotiationTag_);

	edm::ESInputTag ttBuilderTag_("", "TransientTrackBuilder");
	ttBuilderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(ttBuilderTag_);
}

void TrigAnalyzer::beginJob()
{
	nEventsProcessed_ = 0;
	nEventsValidKp_ = 0;
	nEventsValidKm_ = 0;
	nEventsValidpi_ = 0;
	nEventsValidmu_ = 0;
	ftree = new RecoTree(fs->make<TTree>("Events", "Events"));
	ftree->CreateBranches();
}

void TrigAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	edm::Handle<edm::TriggerResults> triggerHandle;
	iEvent.getByToken(triggerToken_, triggerHandle);

	edm::Handle<reco::VertexCollection> primvtxHandle;
	iEvent.getByToken(primvtxToken_, primvtxHandle);

	edm::Handle<reco::GenParticleCollection> GenPartHandle;
	iEvent.getByToken(GenPartToken_, GenPartHandle);

	edm::Handle<pat::PackedCandidateCollection> PFCandHandle;
	iEvent.getByToken(PFCandToken_, PFCandHandle);

	edm::Handle<pat::MuonCollection> MuonHandle;
	iEvent.getByToken(MuonToken_, MuonHandle);

	edm::Handle<edm::Association<reco::GenParticleCollection>> associationHandle;
	iEvent.getByToken(associationToken_, associationHandle);

	edm::ESHandle<TransientTrackBuilder> ttBuilderHandle = iSetup.getHandle(ttBuilderToken_);
	const TransientTrackBuilder *ttBuilder = ttBuilderHandle.product();

	// const edm::TriggerNames &names = iEvent.triggerNames(*triggerHandle);
	// for (UInt_t i = 0; i < names.size(); i++)
	// {
	// 	TString trigName = TString(names.triggerName(i));
	// 	Bool_t pass = (triggerHandle->accept(i) ? true : false);
	// 	// if (trigName.Contains("HLT") && (trigName.Contains("Mu") || trigName.Contains("mu")) && pass)
	// 	edm::LogPrint("TrigAnalyzer") << trigName << ": " << pass;
	// }

	ftree->Init();

	const reco::Vertex &primvtx = primvtxHandle->front();
	if (!primvtx.isValid())
		return;
	nEventsProcessed_++;

	ftree->primvtx_isValid = primvtx.isValid();
	ftree->primvtx_isFake = primvtx.isFake();
	ftree->primvtx_chi2 = primvtx.chi2();
	ftree->primvtx_ndof = primvtx.ndof();
	ftree->primvtx_chi2ndof = primvtx.chi2() / primvtx.ndof();
	ftree->primvtx_vx = primvtx.x();
	ftree->primvtx_vy = primvtx.y();
	ftree->primvtx_vz = primvtx.z();
	ftree->primvtx_vxerr = primvtx.xError();
	ftree->primvtx_vyerr = primvtx.yError();
	ftree->primvtx_vzerr = primvtx.zError();
	GlobalPoint primvtx_pos(primvtx.x(), primvtx.y(), primvtx.z());
	GlobalError primvtx_cov(primvtx.covariance());

	size_t idx_Kp = -1;
	size_t idx_Km = -1;
	size_t idx_pi = -1;
	size_t idx_mu = -1;

	Int_t num_Kp = 0;
	Int_t num_Km = 0;
	Int_t num_pi = 0;
	Int_t num_mu = 0;

	for (size_t i = 0; i < PFCandHandle->size(); i++)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[i];
		if (abs(pf.pdgId()) != 211)
			continue;
		if (!pf.hasTrackDetails())
			continue;

		pat::PackedCandidateRef packedPFCandRef(PFCandHandle, i);
		reco::GenParticleRef genRef = (*associationHandle)[packedPFCandRef];

		if (genRef.isNull())
			continue;

		const reco::GenParticle &gp = *genRef;

		if (abs(gp.pdgId()) == 211)
		{
			if (!hasAncestor(gp, 431, 333))
				continue;
			if (!hasAncestor(gp, 25, 24))
				continue;

			num_pi++;
			idx_pi = i;
		}
		else if (abs(gp.pdgId()) == 321)
		{
			if (!hasAncestor(gp, 333, 321))
				continue;
			if (!hasAncestor(gp, 431, 211))
				continue;
			if (!hasAncestor(gp, 25, 24))
				continue;

			if (pf.charge() == 1)
			{
				num_Kp++;
				idx_Kp = i;
			}
			else if (pf.charge() == -1)
			{
				num_Km++;
				idx_Km = i;
			}
		}
	}

	if (num_Kp != 1)
		return;
	nEventsValidKp_++;
	if (num_Km != 1)
		return;
	nEventsValidKm_++;
	if (num_pi != 1)
		return;
	nEventsValidpi_++;

	size_t idx_K1 = -1;
	size_t idx_K2 = -1;
	if ((*PFCandHandle)[idx_Kp].charge() == (*PFCandHandle)[idx_pi].charge())
	{
		idx_K1 = idx_Kp;
		idx_K2 = idx_Km;
	}
	else if ((*PFCandHandle)[idx_Km].charge() == (*PFCandHandle)[idx_pi].charge())
	{
		idx_K1 = idx_Km;
		idx_K2 = idx_Kp;
	}

	const pat::PackedCandidate &K1_pf = (*PFCandHandle)[idx_K1];
	const pat::PackedCandidate &K2_pf = (*PFCandHandle)[idx_K2];
	const pat::PackedCandidate &pi_pf = (*PFCandHandle)[idx_pi];

	// std::cout << "K1 fromPV: " << K1_pf.fromPV() << ", K1 pvAssociationQuality " << K1_pf.pvAssociationQuality() << std::endl;
	// std::cout << "K2 fromPV: " << K2_pf.fromPV() << ", K2 pvAssociationQuality " << K2_pf.pvAssociationQuality() << std::endl;
	// std::cout << "pi fromPV: " << pi_pf.fromPV() << ", pi pvAssociationQuality " << pi_pf.pvAssociationQuality() << std::endl;

	for (size_t i = 0; i < MuonHandle->size(); i++)
	{
		const pat::Muon &muon = (*MuonHandle)[i];
		if (muon.charge() == pi_pf.charge())
			continue;

		if (!muon.genParticle())
			continue;

		const reco::GenParticle &gp = *(muon.genParticle());

		if (abs(gp.pdgId()) != 13)
			continue;
		if (!hasAncestor(gp, 24, 14))
			continue;
		if (!hasAncestor(gp, 25, 431))
			continue;

		num_mu++;
		idx_mu = i;
	}

	if (num_mu != 1)
		return;
	nEventsValidmu_++;
	const pat::Muon &muon = (*MuonHandle)[idx_mu];

	std::vector<reco::TransientTrack> phi_ttrk = {
		(*ttBuilder).build(K1_pf.pseudoTrack()),
		(*ttBuilder).build(K2_pf.pseudoTrack())};

	TransientVertex phi_tvtx = KVFitter.vertex(phi_ttrk);
	if (!phi_tvtx.isValid())
		return;
	if (!phi_tvtx.hasRefittedTracks())
		return;
	reco::Vertex phi_vtx = reco::Vertex(phi_tvtx);

	std::vector<reco::TransientTrack> phiFit_ttrk = phi_tvtx.refittedTracks();

	std::vector<reco::TransientTrack> Ds_ttrk = {
		phiFit_ttrk[0],
		phiFit_ttrk[1],
		(*ttBuilder).build(pi_pf.pseudoTrack())};

	TransientVertex Ds_tvtx = KVFitter.vertex(Ds_ttrk);

	if (!Ds_tvtx.isValid())
		return;
	if (!Ds_tvtx.hasRefittedTracks())
		return;
	reco::Vertex Ds_vtx = reco::Vertex(Ds_tvtx);

	std::vector<reco::TransientTrack> DsFit_ttrk = Ds_tvtx.refittedTracks();

	TLorentzVector K1_tlv, K2_tlv, pi_tlv, phi_tlv, Ds_tlv;
	K1_tlv.SetXYZM(DsFit_ttrk[0].track().px(), DsFit_ttrk[0].track().py(), DsFit_ttrk[0].track().pz(), Mass_K);
	K2_tlv.SetXYZM(DsFit_ttrk[1].track().px(), DsFit_ttrk[1].track().py(), DsFit_ttrk[1].track().pz(), Mass_K);
	pi_tlv.SetXYZM(DsFit_ttrk[2].track().px(), DsFit_ttrk[2].track().py(), DsFit_ttrk[2].track().pz(), Mass_pi);
	phi_tlv = K1_tlv + K2_tlv;
	Ds_tlv = K1_tlv + K2_tlv + pi_tlv;

	ftree->K1_charge = K1_pf.charge();
	ftree->K1_eta = K1_tlv.Eta();
	ftree->K1_phi = K1_tlv.Phi();
	ftree->K1_px = K1_tlv.Px();
	ftree->K1_py = K1_tlv.Py();
	ftree->K1_pz = K1_tlv.Pz();
	ftree->K1_pt = K1_tlv.Pt();
	ftree->K1_p = K1_tlv.P();

	ftree->K2_charge = K2_pf.charge();
	ftree->K2_eta = K2_tlv.Eta();
	ftree->K2_phi = K2_tlv.Phi();
	ftree->K2_px = K2_tlv.Px();
	ftree->K2_py = K2_tlv.Py();
	ftree->K2_pz = K2_tlv.Pz();
	ftree->K2_pt = K2_tlv.Pt();
	ftree->K2_p = K2_tlv.P();

	ftree->pi_charge = pi_pf.charge();
	ftree->pi_eta = pi_tlv.Eta();
	ftree->pi_phi = pi_tlv.Phi();
	ftree->pi_px = pi_tlv.Px();
	ftree->pi_py = pi_tlv.Py();
	ftree->pi_pz = pi_tlv.Pz();
	ftree->pi_pt = pi_tlv.Pt();
	ftree->pi_p = pi_tlv.P();

	ftree->phi_charge = K1_pf.charge() + K2_pf.charge();
	ftree->phi_eta = phi_tlv.Eta();
	ftree->phi_phi = phi_tlv.Phi();
	ftree->phi_px = phi_tlv.Px();
	ftree->phi_py = phi_tlv.Py();
	ftree->phi_pz = phi_tlv.Pz();
	ftree->phi_pt = phi_tlv.Pt();
	ftree->phi_p = phi_tlv.P();
	ftree->phi_m = phi_tlv.M();

	ftree->Ds_charge = K1_pf.charge() + K2_pf.charge() + pi_pf.charge();
	ftree->Ds_eta = Ds_tlv.Eta();
	ftree->Ds_phi = Ds_tlv.Phi();
	ftree->Ds_px = Ds_tlv.Px();
	ftree->Ds_py = Ds_tlv.Py();
	ftree->Ds_pz = Ds_tlv.Pz();
	ftree->Ds_pt = Ds_tlv.Pt();
	ftree->Ds_p = Ds_tlv.P();
	ftree->Ds_m = Ds_tlv.M();

	ftree->phiFit_isValid = phi_vtx.isValid();
	ftree->phiFit_isFake = phi_vtx.isFake();
	ftree->phiFit_chi2 = phi_vtx.chi2();
	ftree->phiFit_ndof = phi_vtx.ndof();
	ftree->phiFit_chi2ndof = phi_vtx.chi2() / phi_vtx.ndof();
	ftree->phiFit_vx = phi_vtx.x();
	ftree->phiFit_vy = phi_vtx.y();
	ftree->phiFit_vz = phi_vtx.z();
	ftree->phiFit_vxerr = phi_vtx.xError();
	ftree->phiFit_vyerr = phi_vtx.yError();
	ftree->phiFit_vzerr = phi_vtx.zError();

	ftree->DsFit_isValid = Ds_vtx.isValid();
	ftree->DsFit_isFake = Ds_vtx.isFake();
	ftree->DsFit_chi2 = Ds_vtx.chi2();
	ftree->DsFit_ndof = Ds_vtx.ndof();
	ftree->DsFit_chi2ndof = Ds_vtx.chi2() / Ds_vtx.ndof();
	ftree->DsFit_vx = Ds_vtx.x();
	ftree->DsFit_vy = Ds_vtx.y();
	ftree->DsFit_vz = Ds_vtx.z();
	ftree->DsFit_vxerr = Ds_vtx.xError();
	ftree->DsFit_vyerr = Ds_vtx.yError();
	ftree->DsFit_vzerr = Ds_vtx.zError();

	ftree->dR_K1_K2 = reco::deltaR(K1_tlv.Eta(), K1_tlv.Phi(), K2_tlv.Eta(), K2_tlv.Phi());
	ftree->dR_K1_pi = reco::deltaR(K1_tlv.Eta(), K1_tlv.Phi(), pi_tlv.Eta(), pi_tlv.Phi());
	ftree->dR_K1_phi = reco::deltaR(K1_tlv.Eta(), K1_tlv.Phi(), phi_tlv.Eta(), phi_tlv.Phi());
	ftree->dR_K1_Ds = reco::deltaR(K1_tlv.Eta(), K1_tlv.Phi(), Ds_tlv.Eta(), Ds_tlv.Phi());
	ftree->dR_K2_pi = reco::deltaR(K2_tlv.Eta(), K2_tlv.Phi(), pi_tlv.Eta(), pi_tlv.Phi());
	ftree->dR_K2_phi = reco::deltaR(K2_tlv.Eta(), K2_tlv.Phi(), phi_tlv.Eta(), phi_tlv.Phi());
	ftree->dR_K2_Ds = reco::deltaR(K2_tlv.Eta(), K2_tlv.Phi(), Ds_tlv.Eta(), Ds_tlv.Phi());
	ftree->dR_pi_phi = reco::deltaR(pi_tlv.Eta(), pi_tlv.Phi(), phi_tlv.Eta(), phi_tlv.Phi());
	ftree->dR_pi_Ds = reco::deltaR(pi_tlv.Eta(), pi_tlv.Phi(), Ds_tlv.Eta(), Ds_tlv.Phi());
	ftree->dR_phi_Ds = reco::deltaR(phi_tlv.Eta(), phi_tlv.Phi(), Ds_tlv.Eta(), Ds_tlv.Phi());

	ftree->dxy_phi_Ds = std::hypot(phi_vtx.x() - Ds_vtx.x(), phi_vtx.y() - Ds_vtx.y());
	ftree->dz_phi_Ds = abs(phi_vtx.z() - Ds_vtx.z());

	std::pair<bool, Measurement1D> K1_ipresult = IPTools::absoluteImpactParameter3D(DsFit_ttrk[0], primvtx);
	if (K1_ipresult.first)
	{
		ftree->K1_ip = K1_ipresult.second.value();
		ftree->K1_iperr = K1_ipresult.second.error();
		ftree->K1_ipchi2 = K1_ipresult.second.significance();
	}
	std::pair<bool, Measurement1D> K2_ipresult = IPTools::absoluteImpactParameter3D(DsFit_ttrk[1], primvtx);
	if (K2_ipresult.first)
	{
		ftree->K2_ip = K2_ipresult.second.value();
		ftree->K2_iperr = K2_ipresult.second.error();
		ftree->K2_ipchi2 = K2_ipresult.second.significance();
	}
	std::pair<bool, Measurement1D> pi_ipresult = IPTools::absoluteImpactParameter3D(DsFit_ttrk[2], primvtx);
	if (pi_ipresult.first)
	{
		ftree->pi_ip = pi_ipresult.second.value();
		ftree->pi_iperr = pi_ipresult.second.error();
		ftree->pi_ipchi2 = pi_ipresult.second.significance();
	}

	GlobalPoint Ds_pos(Ds_vtx.x(), Ds_vtx.y(), Ds_vtx.z());
	GlobalError Ds_cov(Ds_vtx.covariance());

	GlobalVector displacement = Ds_pos - primvtx_pos;
	GlobalVector PDirec(Ds_tlv.Px(), Ds_tlv.Py(), Ds_tlv.Pz());
	ftree->dira = displacement.dot(PDirec) / (displacement.mag() * PDirec.mag());
	ftree->dira_angle = std::acos(ftree->dira);

	// // VertexDistance vtxdist;
	VertexDistanceXY vtxdistXY;
	Measurement1D Ds_FDXY = vtxdistXY.signedDistance(primvtx, Ds_vtx, PDirec);
	ftree->Ds_FDxy = Ds_FDXY.value();
	ftree->Ds_FDxyerr = Ds_FDXY.error();
	ftree->Ds_FDxysig = Ds_FDXY.significance();

	VertexDistance3D vtxdist3D;
	Measurement1D Ds_FD_3D = vtxdist3D.signedDistance(primvtx, Ds_vtx, PDirec);
	ftree->Ds_FD = Ds_FD_3D.value();
	ftree->Ds_FDerr = Ds_FD_3D.error();
	ftree->Ds_FDsig = Ds_FD_3D.significance();
	// // Measurement1D Ds_FD = vtxdist.distance(primvtx, Ds_vtx);

	// // Measurement1D Ds_FD_3D = vtxdist3D.signedDistance(primvtx, Ds_vtx, PDirec);
	// // Measurement1D Ds_FD_3D = vtxdist3D.distance(primvtx_pos, primvtx_cov, Ds_pos, Ds_cov);
	// // ftree->Ds_FD_xy = vtxdist.distance(primvtx, Ds_vtx);

	Float_t IsoR03_sumChargedHadronPt = 0;
	Float_t IsoR03_sumNeutralHadronEt = 0;
	Float_t IsoR03_sumPhotonEt = 0;
	Float_t IsoR03_sumPUPt = 0;
	Float_t IsoR04_sumChargedHadronPt = 0;
	Float_t IsoR04_sumNeutralHadronEt = 0;
	Float_t IsoR04_sumPhotonEt = 0;
	Float_t IsoR04_sumPUPt = 0;

	for (size_t i = 0; i < PFCandHandle->size(); i++)
	{
		if (i == idx_K1)
			continue;
		if (i == idx_K2)
			continue;
		if (i == idx_pi)
			continue;
		const pat::PackedCandidate &pf = (*PFCandHandle)[i];

		if (reco::deltaR(Ds_tlv.Eta(), Ds_tlv.Phi(), pf.eta(), pf.phi()) < 0.3)
		{
			if (abs(pf.pdgId()) == 211 && (pf.fromPV() > 1))
				IsoR03_sumChargedHadronPt += pf.pt();
			else if (pf.pdgId() == 130)
				IsoR03_sumNeutralHadronEt += pf.et();
			else if (pf.pdgId() == 22)
				IsoR03_sumPhotonEt += pf.et();
			if (pf.fromPV() <= 1)
				IsoR03_sumPUPt += pf.pt();
		}
		if (reco::deltaR(Ds_tlv.Eta(), Ds_tlv.Phi(), pf.eta(), pf.phi()) < 0.4)
		{
			if (abs(pf.pdgId()) == 211 && (pf.fromPV() > 1))
				IsoR04_sumChargedHadronPt += pf.pt();
			else if (pf.pdgId() == 130)
				IsoR04_sumNeutralHadronEt += pf.et();
			else if (pf.pdgId() == 22)
				IsoR04_sumPhotonEt += pf.et();
			if (pf.fromPV() <= 1)
				IsoR04_sumPUPt += pf.pt();
		}
	}

	ftree->Ds_IsoR03_sumChargedHadronPt = IsoR03_sumChargedHadronPt;
	ftree->Ds_IsoR03_sumNeutralHadronEt = IsoR03_sumNeutralHadronEt;
	ftree->Ds_IsoR03_sumPhotonEt = IsoR03_sumPhotonEt;
	ftree->Ds_IsoR03_sumPUPt = IsoR03_sumPUPt;
	ftree->Ds_IsoR03_PFIso = IsoR03_sumChargedHadronPt + std::max(0.0, IsoR03_sumChargedHadronPt + IsoR03_sumPhotonEt - 0.5 * IsoR03_sumPUPt);
	ftree->Ds_IsoR03_RelIso = (IsoR03_sumChargedHadronPt + std::max(0.0, IsoR03_sumChargedHadronPt + IsoR03_sumPhotonEt - 0.5 * IsoR03_sumPUPt)) / Ds_tlv.Pt();

	ftree->Ds_IsoR04_sumChargedHadronPt = IsoR04_sumChargedHadronPt;
	ftree->Ds_IsoR04_sumNeutralHadronEt = IsoR04_sumNeutralHadronEt;
	ftree->Ds_IsoR04_sumPhotonEt = IsoR04_sumPhotonEt;
	ftree->Ds_IsoR04_sumPUPt = IsoR04_sumPUPt;
	ftree->Ds_IsoR04_PFIso = IsoR04_sumChargedHadronPt + std::max(0.0, IsoR04_sumChargedHadronPt + IsoR04_sumPhotonEt - 0.5 * IsoR04_sumPUPt);
	ftree->Ds_IsoR04_RelIso = (IsoR04_sumChargedHadronPt + std::max(0.0, IsoR04_sumChargedHadronPt + IsoR04_sumPhotonEt - 0.5 * IsoR04_sumPUPt)) / Ds_tlv.Pt();

	ftree->mu_charge = muon.charge();
	ftree->mu_eta = muon.eta();
	ftree->mu_phi = muon.phi();
	ftree->mu_px = muon.px();
	ftree->mu_py = muon.py();
	ftree->mu_pz = muon.pz();
	ftree->mu_pt = muon.pt();
	ftree->mu_p = muon.p();
	ftree->mu_isHighPt = muon.isHighPtMuon(primvtx);
	ftree->mu_isLoose = muon.isLooseMuon();
	ftree->mu_isMedium = muon.isMediumMuon();
	ftree->mu_isSoft = muon.isSoftMuon(primvtx);
	ftree->mu_isTight = muon.isTightMuon(primvtx);
	ftree->mu_isPF = muon.isPFMuon();
	ftree->mu_isTracker = muon.isTrackerMuon();
	ftree->mu_isGlobal = muon.isGlobalMuon();

	ftree->dR_K1_mu = reco::deltaR(K1_tlv.Eta(), K1_tlv.Phi(), muon.eta(), muon.phi());
	ftree->dR_K2_mu = reco::deltaR(K2_tlv.Eta(), K2_tlv.Phi(), muon.eta(), muon.phi());
	ftree->dR_pi_mu = reco::deltaR(pi_tlv.Eta(), pi_tlv.Phi(), muon.eta(), muon.phi());
	ftree->dR_phi_mu = reco::deltaR(phi_tlv.Eta(), phi_tlv.Phi(), muon.eta(), muon.phi());
	ftree->dR_Ds_mu = reco::deltaR(Ds_tlv.Eta(), Ds_tlv.Phi(), muon.eta(), muon.phi());

	std::pair<bool, Measurement1D> mu_ipresult = IPTools::absoluteImpactParameter3D((*ttBuilder).build(muon.muonBestTrack()), primvtx);
	if (mu_ipresult.first)
	{
		ftree->mu_ip = mu_ipresult.second.value();
		ftree->mu_iperr = mu_ipresult.second.error();
		ftree->mu_ipchi2 = mu_ipresult.second.significance();
	}

	ftree->mu_IsoR03_sumChargedHadronPt = muon.pfIsolationR03().sumChargedHadronPt;
	ftree->mu_IsoR03_sumNeutralHadronEt = muon.pfIsolationR03().sumNeutralHadronEt;
	ftree->mu_IsoR03_sumPhotonEt = muon.pfIsolationR03().sumPhotonEt;
	ftree->mu_IsoR03_sumPUPt = muon.pfIsolationR03().sumPUPt;
	ftree->mu_IsoR03_PFIso = muon.pfIsolationR03().sumChargedHadronPt + std::max(0.0, muon.pfIsolationR03().sumNeutralHadronEt + muon.pfIsolationR03().sumPhotonEt - 0.5 * muon.pfIsolationR03().sumPUPt);
	ftree->mu_IsoR03_RelIso = (muon.pfIsolationR03().sumChargedHadronPt + std::max(0.0, muon.pfIsolationR03().sumNeutralHadronEt + muon.pfIsolationR03().sumPhotonEt - 0.5 * muon.pfIsolationR03().sumPUPt)) / muon.pt();
	ftree->mu_IsoR04_sumChargedHadronPt = muon.pfIsolationR04().sumChargedHadronPt;
	ftree->mu_IsoR04_sumNeutralHadronEt = muon.pfIsolationR04().sumNeutralHadronEt;
	ftree->mu_IsoR04_sumPhotonEt = muon.pfIsolationR04().sumPhotonEt;
	ftree->mu_IsoR04_sumPUPt = muon.pfIsolationR04().sumPUPt;
	ftree->mu_IsoR04_PFIso = muon.pfIsolationR04().sumChargedHadronPt + std::max(0.0, muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt);
	ftree->mu_IsoR04_RelIso = (muon.pfIsolationR04().sumChargedHadronPt + std::max(0.0, muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt();

	ftree->tree->Fill();

	const edm::TriggerNames &names = iEvent.triggerNames(*triggerHandle);

	for (Size_t i = 0; i < names.size(); i++)
	{
		TString trigName = names.triggerName(i);
		if (!trigName.Contains("HLT"))
			continue;
		if (!trigName.Contains("Mu") && !trigName.Contains("mu") && !trigName.Contains("MET"))
			continue;
		Bool_t pass = triggerHandle->accept(i);
		triggerTotal_[trigName]++;
		if (pass)
			triggerCounts_[trigName]++;
	}
}

bool TrigAnalyzer::hasAncestor(const reco::GenParticle &gp, const int motherid, const int otherparticleid) const
{
	bool hasancestor = false;

	if (gp.numberOfMothers() == 0)
		return false;

	if (gp.numberOfMothers() == 1 && abs(gp.mother(0)->pdgId()) == motherid)
	{
		if (gp.mother(0)->numberOfDaughters() == 2)
		{
			if (abs(gp.mother(0)->daughter(0)->pdgId()) == otherparticleid || abs(gp.mother(0)->daughter(1)->pdgId()) == otherparticleid)
				hasancestor = true;
		}
	}

	if (hasancestor)
		return true;
	else
		return hasAncestor(*gp.motherRef(0), motherid, otherparticleid);
}

void TrigAnalyzer::endJob()
{
	// std::cout << nEventsProcessed_ << ", Kp: " << nEventsValidKp_ << ", Km " << nEventsValidKm_ << ", pi" << nEventsValidpi_ << ", mu" << nEventsValidmu_ << std::endl;

	edm::LogPrint("TrigAnalyzer") << "\n========== Trigger Statistics ==========";
	edm::LogPrint("TrigAnalyzer") << std::setw(60) << "Trigger Name"
								  << std::setw(15) << "Fired"
								  << std::setw(15) << "Total"
								  << std::setw(15) << "Efficiency";
	edm::LogPrint("TrigAnalyzer") << std::string(105, '-');

	std::vector<std::pair<TString, double>> triggerEfficiencies;

	for (const auto &entry : triggerCounts_)
	{
		const TString trigName = entry.first;
		int fired = entry.second;
		int total = triggerTotal_[trigName];
		double efficiency = (total > 0) ? (double)fired / total * 100.0 : 0.0;

		triggerEfficiencies.push_back({trigName, efficiency});
	}

	std::sort(triggerEfficiencies.begin(), triggerEfficiencies.end(),
			  [](const auto &a, const auto &b)
			  { return a.second > b.second; });

	for (const auto &item : triggerEfficiencies)
	{
		const TString &trigName = item.first;
		int fired = triggerCounts_[trigName];
		int total = triggerTotal_[trigName];
		double efficiency = item.second;

		edm::LogPrint("TrigAnalyzer") << std::setw(80) << trigName
									  << std::setw(15) << fired
									  << std::setw(15) << total
									  << std::setw(14) << std::fixed << std::setprecision(2)
									  << efficiency << "%";
	}

	edm::LogPrint("TrigAnalyzer") << "========================================\n";
}

void TrigAnalyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(TrigAnalyzer);
