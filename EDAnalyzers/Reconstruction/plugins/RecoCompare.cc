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

#include "EDAnalyzers/Reconstruction/interface/AllRecoTree.h"

#include "TLorentzVector.h"
#include "TTree.h"

class RecoCompare : public edm::one::EDAnalyzer<>
{
public:
	explicit RecoCompare(const edm::ParameterSet &);
	~RecoCompare() override {};

	static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
	void beginJob() override;
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	void endJob() override;

	Bool_t hasAncestor(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const;

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

	// std::map<TString, Int_t> triggerCounts_;
	// std::map<TString, Int_t> triggerTotal_;

	const edm::Service<TFileService> fs;
	AllRecoTree *ftree;

	KalmanVertexFitter KVFitter;
};

RecoCompare::RecoCompare(const edm::ParameterSet &iConfig)
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

void RecoCompare::beginJob()
{
	ftree = new AllRecoTree(fs->make<TTree>("Events", "Events"));
	ftree->CreateBranches();
}

void RecoCompare::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
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

	ftree->Init();

	const reco::Vertex &primvtx = primvtxHandle->front();
	if (!primvtx.isValid())
		return;

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

	Int_t Kp_match_idx = null;
	Int_t Km_match_idx = null;
	Int_t pi_match_idx = null;
	Int_t mu_match_idx = null;

	Float_t Kp_deltaR = 1.0;
	Float_t Km_deltaR = 1.0;
	Float_t pi_deltaR = 1.0;
	Float_t mu_deltaR = 1.0;

	for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[i];
		const pat::PackedCandidateRef pfRef(PFCandHandle, i);
		const reco::GenParticleRef genRef = (*associationHandle)[pfRef];
		if (genRef.isNull())
			continue;
		const reco::GenParticle &gp = *genRef;
		Float_t deltaR = reco::deltaR(gp.eta(), gp.phi(), pf.eta(), pf.phi());
		if (deltaR > 0.05)
			continue;
		Float_t deltapt = abs((gp.pt() - pf.pt()) / gp.pt());
		if (deltapt > 0.3)
			continue;

		if ((gp.pdgId() == 321))
		{
			if (hasAncestor(gp, 333, 321, -321) && hasAncestor(gp, 431, 333, 211) && hasAncestor(gp, 25, 431, -24))
			{
				if (deltaR < Kp_deltaR)
				{
					Kp_match_idx = i;
					Kp_deltaR = deltaR;
				}
			}
		}
		else if ((gp.pdgId() == -321))
		{
			if (hasAncestor(gp, 333, 321, -321) && hasAncestor(gp, 431, 333, 211) && hasAncestor(gp, 25, 431, -24))
			{
				if (deltaR < Km_deltaR)
				{
					Km_match_idx = i;
					Km_deltaR = deltaR;
				}
			}
		}
		else if ((gp.pdgId() == 211))
		{
			if (hasAncestor(gp, 431, 333, 211) && hasAncestor(gp, 25, 431, -24))
			{
				if (deltaR < pi_deltaR)
				{
					pi_match_idx = i;
					pi_deltaR = deltaR;
				}
			}
		}
	}

	for (Int_t i = 0; i < static_cast<Int_t>(MuonHandle->size()); i++)
	{
		const pat::Muon &muon = (*MuonHandle)[i];
		if (!muon.genParticle())
			continue;
		const reco::GenParticle &gp = *(muon.genParticle());
		Float_t deltaR = reco::deltaR(gp.eta(), gp.phi(), muon.eta(), muon.phi());
		if (deltaR > 0.05)
			continue;
		Float_t deltapt = abs((gp.pt() - muon.pt()) / gp.pt());
		if (deltapt > 0.3)
			continue;

		if ((gp.pdgId() == 13) && hasAncestor(gp, -24, 13, -14) && hasAncestor(gp, 25, 431, -24))
		{
			if (deltaR < mu_deltaR)
			{
				mu_match_idx = i;
				mu_deltaR = deltaR;
			}
		}
	}

	std::vector<Int_t> idx_pf;
	idx_pf.clear();

	for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
	{
		const pat::PackedCandidate &pf = (*PFCandHandle)[i];
		if (!pf.hasTrackDetails())
			continue;
		if (!pf.trackHighPurity())
			continue;
		if (pf.pt() < 0.5)
			continue;
		if (pf.p() < 1)
			continue;
		if (abs(pf.pdgId()) != 211)
			continue;

		idx_pf.push_back(i);
	}

	if (idx_pf.size() < 3)
		return;

	for (Int_t iKp = 0; iKp < static_cast<Int_t>(idx_pf.size()); iKp++)
	{
		const Int_t idx_Kp = idx_pf[iKp];

		const pat::PackedCandidate &Kp_pf = (*PFCandHandle)[idx_Kp];
		if (Kp_pf.charge() != 1)
			continue;

		for (Int_t iKm = 0; iKm < static_cast<Int_t>(idx_pf.size()); iKm++)
		{
			if (iKp == iKm)
				continue;
			const Int_t idx_Km = idx_pf[iKm];

			const pat::PackedCandidate &Km_pf = (*PFCandHandle)[idx_Km];
			if (Km_pf.charge() != -1)
				continue;

			if (reco::deltaR(Kp_pf.eta(), Kp_pf.phi(), Km_pf.eta(), Km_pf.phi()) > 0.15)
				continue;

			std::vector<reco::TransientTrack> phi_ttrk = {
				(*ttBuilder).build(Kp_pf.pseudoTrack()),
				(*ttBuilder).build(Km_pf.pseudoTrack())};

			TransientVertex phi_tvtx = KVFitter.vertex(phi_ttrk);
			if (!phi_tvtx.isValid())
				continue;
			if (!phi_tvtx.hasRefittedTracks())
				continue;
			if (phi_tvtx.normalisedChiSquared() < 0)
				continue;
			if (phi_tvtx.normalisedChiSquared() > 10)
				continue;

			reco::Vertex phi_vtx = reco::Vertex(phi_tvtx);
			std::vector<reco::TransientTrack> phiFit_ttrk = phi_tvtx.refittedTracks();
			TLorentzVector phiFit_Kp_tlv, phiFit_Km_tlv, phiFit_phi_tlv;
			phiFit_Kp_tlv.SetXYZM(phiFit_ttrk[0].track().px(), phiFit_ttrk[0].track().py(), phiFit_ttrk[0].track().pz(), Mass_K);
			phiFit_Km_tlv.SetXYZM(phiFit_ttrk[1].track().px(), phiFit_ttrk[1].track().py(), phiFit_ttrk[1].track().pz(), Mass_K);
			phiFit_phi_tlv = phiFit_Kp_tlv + phiFit_Km_tlv;
			if (reco::deltaR(phiFit_Kp_tlv.Eta(), phiFit_Kp_tlv.Phi(), phiFit_phi_tlv.Eta(), phiFit_phi_tlv.Phi()) > 0.1)
				continue;
			if (reco::deltaR(phiFit_Km_tlv.Eta(), phiFit_Km_tlv.Phi(), phiFit_phi_tlv.Eta(), phiFit_phi_tlv.Phi()) > 0.1)
				continue;
			if ((phiFit_phi_tlv.M() - Mass_phi) > 0.1)
				continue;

			for (Int_t ipi = 0; ipi < static_cast<Int_t>(idx_pf.size()); ipi++)
			{
				if ((iKp == ipi) || (iKm == ipi))
					continue;
				const Int_t idx_pi = idx_pf[ipi];
				const pat::PackedCandidate &pi_pf = (*PFCandHandle)[idx_pi];

				if (reco::deltaR(Kp_pf.eta(), Kp_pf.phi(), pi_pf.eta(), pi_pf.phi()) > 0.6)
					continue;
				if (reco::deltaR(Km_pf.eta(), Km_pf.phi(), pi_pf.eta(), pi_pf.phi()) > 0.6)
					continue;
				if (reco::deltaR(phiFit_phi_tlv.Eta(), phiFit_phi_tlv.Phi(), pi_pf.eta(), pi_pf.phi()) > 0.6)
					continue;

				std::vector<reco::TransientTrack> Ds_ttrk = {
					phiFit_ttrk[0],
					phiFit_ttrk[1],
					(*ttBuilder).build(pi_pf.pseudoTrack())};

				TransientVertex Ds_tvtx = KVFitter.vertex(Ds_ttrk);
				if (!Ds_tvtx.isValid())
					continue;
				if (!Ds_tvtx.hasRefittedTracks())
					continue;
				if (Ds_tvtx.normalisedChiSquared() < 0)
					continue;
				if (Ds_tvtx.normalisedChiSquared() > 10)
					continue;
				reco::Vertex Ds_vtx = reco::Vertex(Ds_tvtx);

				std::vector<reco::TransientTrack> DsFit_ttrk = Ds_tvtx.refittedTracks();
				TLorentzVector Kp_tlv, Km_tlv, pi_tlv, phi_tlv, Ds_tlv;
				Kp_tlv.SetXYZM(DsFit_ttrk[0].track().px(), DsFit_ttrk[0].track().py(), DsFit_ttrk[0].track().pz(), Mass_K);
				Km_tlv.SetXYZM(DsFit_ttrk[1].track().px(), DsFit_ttrk[1].track().py(), DsFit_ttrk[1].track().pz(), Mass_K);
				pi_tlv.SetXYZM(DsFit_ttrk[2].track().px(), DsFit_ttrk[2].track().py(), DsFit_ttrk[2].track().pz(), Mass_pi);
				phi_tlv = Kp_tlv + Km_tlv;
				Ds_tlv = Kp_tlv + Km_tlv + pi_tlv;
				TLorentzVector phi_tlv_phiM, Ds_tlv_phiM;
				phi_tlv_phiM.SetXYZM(phi_tlv.Px(), phi_tlv.Py(), phi_tlv.Pz(), Mass_phi);
				Ds_tlv_phiM = phi_tlv_phiM + pi_tlv;

				if (reco::deltaR(Kp_tlv.Eta(), Kp_tlv.Phi(), Ds_tlv.Eta(), Ds_tlv.Phi()) > 0.4)
					continue;
				if (reco::deltaR(Km_tlv.Eta(), Km_tlv.Phi(), Ds_tlv.Eta(), Ds_tlv.Phi()) > 0.4)
					continue;
				if (reco::deltaR(pi_tlv.Eta(), pi_tlv.Phi(), Ds_tlv.Eta(), Ds_tlv.Phi()) > 0.6)
					continue;
				if (reco::deltaR(phi_tlv.Eta(), phi_tlv.Phi(), Ds_tlv.Eta(), Ds_tlv.Phi()) > 0.4)
					continue;

				if ((phi_tlv.M() - Mass_phi) > 0.05)
					continue;
				if ((Ds_tlv.M() - Mass_Ds) > 0.15)
					continue;

				Bool_t Kp_match = (idx_Kp == Kp_match_idx);
				Bool_t Km_match = (idx_Km == Km_match_idx);
				Bool_t pi_match = (idx_pi == pi_match_idx);

				std::pair<Bool_t, Measurement1D> Kp_ipresult = IPTools::absoluteImpactParameter3D(DsFit_ttrk[0], primvtx);
				std::pair<Bool_t, Measurement1D> Km_ipresult = IPTools::absoluteImpactParameter3D(DsFit_ttrk[1], primvtx);
				std::pair<Bool_t, Measurement1D> pi_ipresult = IPTools::absoluteImpactParameter3D(DsFit_ttrk[2], primvtx);

				GlobalPoint Ds_pos(Ds_vtx.x(), Ds_vtx.y(), Ds_vtx.z());
				GlobalError Ds_cov(Ds_vtx.covariance());

				GlobalVector displacement = Ds_pos - primvtx_pos;
				GlobalVector PDirec(Ds_tlv.Px(), Ds_tlv.Py(), Ds_tlv.Pz());

				VertexDistanceXY vtxdistXY;
				Measurement1D Ds_FDXY = vtxdistXY.signedDistance(primvtx, Ds_vtx, PDirec);
				VertexDistance3D vtxdist3D;
				Measurement1D Ds_FD_3D = vtxdist3D.signedDistance(primvtx, Ds_vtx, PDirec);

				Float_t direction_angle = displacement.dot(PDirec) / (displacement.mag() * PDirec.mag());

				Float_t IsoR03_sumChargedHadronPt = 0;
				Float_t IsoR03_sumNeutralHadronEt = 0;
				Float_t IsoR03_sumPhotonEt = 0;
				Float_t IsoR03_sumPUPt = 0;
				Float_t IsoR04_sumChargedHadronPt = 0;
				Float_t IsoR04_sumNeutralHadronEt = 0;
				Float_t IsoR04_sumPhotonEt = 0;
				Float_t IsoR04_sumPUPt = 0;

				for (Int_t i = 0; i < static_cast<Int_t>(PFCandHandle->size()); i++)
				{
					if (i == idx_Kp)
						continue;
					if (i == idx_Km)
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
				Float_t IsoR03_PFIso = IsoR03_sumChargedHadronPt + std::max(0.0, IsoR03_sumNeutralHadronEt + IsoR03_sumPhotonEt - 0.5 * IsoR03_sumPUPt);
				Float_t IsoR03_RelIso = (IsoR03_sumChargedHadronPt + std::max(0.0, IsoR03_sumNeutralHadronEt + IsoR03_sumPhotonEt - 0.5 * IsoR03_sumPUPt)) / Ds_tlv.Pt();
				Float_t IsoR04_PFIso = IsoR04_sumChargedHadronPt + std::max(0.0, IsoR04_sumNeutralHadronEt + IsoR04_sumPhotonEt - 0.5 * IsoR04_sumPUPt);
				Float_t IsoR04_RelIso = (IsoR04_sumChargedHadronPt + std::max(0.0, IsoR04_sumNeutralHadronEt + IsoR04_sumPhotonEt - 0.5 * IsoR04_sumPUPt)) / Ds_tlv.Pt();

				ftree->Kp_charge.push_back(Kp_pf.charge());
				ftree->Kp_eta.push_back(Kp_tlv.Eta());
				ftree->Kp_phi.push_back(Kp_tlv.Phi());
				ftree->Kp_px.push_back(Kp_tlv.Px());
				ftree->Kp_py.push_back(Kp_tlv.Py());
				ftree->Kp_pz.push_back(Kp_tlv.Pz());
				ftree->Kp_pt.push_back(Kp_tlv.Pt());
				ftree->Kp_p.push_back(Kp_tlv.P());

				ftree->Km_charge.push_back(Km_pf.charge());
				ftree->Km_eta.push_back(Km_tlv.Eta());
				ftree->Km_phi.push_back(Km_tlv.Phi());
				ftree->Km_px.push_back(Km_tlv.Px());
				ftree->Km_py.push_back(Km_tlv.Py());
				ftree->Km_pz.push_back(Km_tlv.Pz());
				ftree->Km_pt.push_back(Km_tlv.Pt());
				ftree->Km_p.push_back(Km_tlv.P());

				ftree->pi_charge.push_back(pi_pf.charge());
				ftree->pi_eta.push_back(pi_tlv.Eta());
				ftree->pi_phi.push_back(pi_tlv.Phi());
				ftree->pi_px.push_back(pi_tlv.Px());
				ftree->pi_py.push_back(pi_tlv.Py());
				ftree->pi_pz.push_back(pi_tlv.Pz());
				ftree->pi_pt.push_back(pi_tlv.Pt());
				ftree->pi_p.push_back(pi_tlv.P());

				ftree->phi_charge.push_back(Kp_pf.charge() + Km_pf.charge());
				ftree->phi_eta.push_back(phi_tlv.Eta());
				ftree->phi_phi.push_back(phi_tlv.Phi());
				ftree->phi_px.push_back(phi_tlv.Px());
				ftree->phi_py.push_back(phi_tlv.Py());
				ftree->phi_pz.push_back(phi_tlv.Pz());
				ftree->phi_pt.push_back(phi_tlv.Pt());
				ftree->phi_p.push_back(phi_tlv.P());
				ftree->phi_m.push_back(phi_tlv.M());

				ftree->Ds_charge.push_back(Kp_pf.charge() + Km_pf.charge() + pi_pf.charge());
				ftree->Ds_eta.push_back(Ds_tlv.Eta());
				ftree->Ds_phi.push_back(Ds_tlv.Phi());
				ftree->Ds_px.push_back(Ds_tlv.Px());
				ftree->Ds_py.push_back(Ds_tlv.Py());
				ftree->Ds_pz.push_back(Ds_tlv.Pz());
				ftree->Ds_pt.push_back(Ds_tlv.Pt());
				ftree->Ds_p.push_back(Ds_tlv.P());
				ftree->Ds_m.push_back(Ds_tlv.M());
				ftree->Ds_m_phiM.push_back(Ds_tlv_phiM.M());

				ftree->phiFit_isValid.push_back(phi_vtx.isValid());
				ftree->phiFit_isFake.push_back(phi_vtx.isFake());
				ftree->phiFit_chi2.push_back(phi_vtx.chi2());
				ftree->phiFit_ndof.push_back(phi_vtx.ndof());
				ftree->phiFit_chi2ndof.push_back(phi_vtx.chi2() / phi_vtx.ndof());
				ftree->phiFit_vx.push_back(phi_vtx.x());
				ftree->phiFit_vy.push_back(phi_vtx.y());
				ftree->phiFit_vz.push_back(phi_vtx.z());
				ftree->phiFit_vxerr.push_back(phi_vtx.xError());
				ftree->phiFit_vyerr.push_back(phi_vtx.yError());
				ftree->phiFit_vzerr.push_back(phi_vtx.zError());

				ftree->DsFit_isValid.push_back(Ds_vtx.isValid());
				ftree->DsFit_isFake.push_back(Ds_vtx.isFake());
				ftree->DsFit_chi2.push_back(Ds_vtx.chi2());
				ftree->DsFit_ndof.push_back(Ds_vtx.ndof());
				ftree->DsFit_chi2ndof.push_back(Ds_vtx.chi2() / Ds_vtx.ndof());
				ftree->DsFit_vx.push_back(Ds_vtx.x());
				ftree->DsFit_vy.push_back(Ds_vtx.y());
				ftree->DsFit_vz.push_back(Ds_vtx.z());
				ftree->DsFit_vxerr.push_back(Ds_vtx.xError());
				ftree->DsFit_vyerr.push_back(Ds_vtx.yError());
				ftree->DsFit_vzerr.push_back(Ds_vtx.zError());

				ftree->orig_dR_Kp_Km.push_back(reco::deltaR(Kp_pf.eta(), Kp_pf.phi(), Km_pf.eta(), Km_pf.phi()));
				ftree->orig_dR_Kp_pi.push_back(reco::deltaR(Kp_pf.eta(), Kp_pf.phi(), pi_pf.eta(), pi_pf.phi()));
				ftree->orig_dR_Km_pi.push_back(reco::deltaR(Km_pf.eta(), Km_pf.phi(), pi_pf.eta(), pi_pf.phi()));

				ftree->phiFit_dR_Kp_Km.push_back(reco::deltaR(phiFit_Kp_tlv.Eta(), phiFit_Kp_tlv.Phi(), phiFit_Km_tlv.Eta(), phiFit_Km_tlv.Phi()));
				ftree->phiFit_dR_Kp_pi.push_back(reco::deltaR(phiFit_Kp_tlv.Eta(), phiFit_Kp_tlv.Phi(), pi_pf.eta(), pi_pf.phi()));
				ftree->phiFit_dR_Kp_phi.push_back(reco::deltaR(phiFit_Kp_tlv.Eta(), phiFit_Kp_tlv.Phi(), phiFit_phi_tlv.Eta(), phiFit_phi_tlv.Phi()));
				ftree->phiFit_dR_Km_pi.push_back(reco::deltaR(phiFit_Km_tlv.Eta(), phiFit_Km_tlv.Phi(), pi_pf.eta(), pi_pf.phi()));
				ftree->phiFit_dR_Km_phi.push_back(reco::deltaR(phiFit_Km_tlv.Eta(), phiFit_Km_tlv.Phi(), phiFit_phi_tlv.Eta(), phiFit_phi_tlv.Phi()));
				ftree->phiFit_dR_pi_phi.push_back(reco::deltaR(pi_pf.eta(), pi_pf.phi(), phiFit_phi_tlv.Eta(), phiFit_phi_tlv.Phi()));

				ftree->dR_Kp_Km.push_back(reco::deltaR(Kp_tlv.Eta(), Kp_tlv.Phi(), Km_tlv.Eta(), Km_tlv.Phi()));
				ftree->dR_Kp_pi.push_back(reco::deltaR(Kp_tlv.Eta(), Kp_tlv.Phi(), pi_tlv.Eta(), pi_tlv.Phi()));
				ftree->dR_Kp_phi.push_back(reco::deltaR(Kp_tlv.Eta(), Kp_tlv.Phi(), phi_tlv.Eta(), phi_tlv.Phi()));
				ftree->dR_Kp_Ds.push_back(reco::deltaR(Kp_tlv.Eta(), Kp_tlv.Phi(), Ds_tlv.Eta(), Ds_tlv.Phi()));
				ftree->dR_Km_pi.push_back(reco::deltaR(Km_tlv.Eta(), Km_tlv.Phi(), pi_tlv.Eta(), pi_tlv.Phi()));
				ftree->dR_Km_phi.push_back(reco::deltaR(Km_tlv.Eta(), Km_tlv.Phi(), phi_tlv.Eta(), phi_tlv.Phi()));
				ftree->dR_Km_Ds.push_back(reco::deltaR(Km_tlv.Eta(), Km_tlv.Phi(), Ds_tlv.Eta(), Ds_tlv.Phi()));
				ftree->dR_pi_phi.push_back(reco::deltaR(pi_tlv.Eta(), pi_tlv.Phi(), phi_tlv.Eta(), phi_tlv.Phi()));
				ftree->dR_pi_Ds.push_back(reco::deltaR(pi_tlv.Eta(), pi_tlv.Phi(), Ds_tlv.Eta(), Ds_tlv.Phi()));
				ftree->dR_phi_Ds.push_back(reco::deltaR(phi_tlv.Eta(), phi_tlv.Phi(), Ds_tlv.Eta(), Ds_tlv.Phi()));

				ftree->dxy_phi_Ds.push_back(std::hypot(phi_vtx.x() - Ds_vtx.x(), phi_vtx.y() - Ds_vtx.y()));
				ftree->dz_phi_Ds.push_back(abs(phi_vtx.z() - Ds_vtx.z()));

				if (Kp_ipresult.first)
				{
					ftree->Kp_ip.push_back(Kp_ipresult.second.value());
					ftree->Kp_iperr.push_back(Kp_ipresult.second.error());
					ftree->Kp_ipchi2.push_back(Kp_ipresult.second.significance());
				}
				else
				{
					ftree->Kp_ip.push_back(null);
					ftree->Kp_iperr.push_back(null);
					ftree->Kp_ipchi2.push_back(null);
				}
				if (Km_ipresult.first)
				{
					ftree->Km_ip.push_back(Km_ipresult.second.value());
					ftree->Km_iperr.push_back(Km_ipresult.second.error());
					ftree->Km_ipchi2.push_back(Km_ipresult.second.significance());
				}
				else
				{
					ftree->Km_ip.push_back(null);
					ftree->Km_iperr.push_back(null);
					ftree->Km_ipchi2.push_back(null);
				}
				if (pi_ipresult.first)
				{
					ftree->pi_ip.push_back(pi_ipresult.second.value());
					ftree->pi_iperr.push_back(pi_ipresult.second.error());
					ftree->pi_ipchi2.push_back(pi_ipresult.second.significance());
				}
				else
				{
					ftree->pi_ip.push_back(null);
					ftree->pi_iperr.push_back(null);
					ftree->pi_ipchi2.push_back(null);
				}

				ftree->dira.push_back(direction_angle);
				ftree->dira_angle.push_back(std::acos(direction_angle));

				ftree->Ds_FDxy.push_back(Ds_FDXY.value());
				ftree->Ds_FDxyerr.push_back(Ds_FDXY.error());
				ftree->Ds_FDxysig.push_back(Ds_FDXY.significance());

				ftree->Ds_FD.push_back(Ds_FD_3D.value());
				ftree->Ds_FDerr.push_back(Ds_FD_3D.error());
				ftree->Ds_FDsig.push_back(Ds_FD_3D.significance());

				ftree->Ds_IsoR03_sumChargedHadronPt.push_back(IsoR03_sumChargedHadronPt);
				ftree->Ds_IsoR03_sumNeutralHadronEt.push_back(IsoR03_sumNeutralHadronEt);
				ftree->Ds_IsoR03_sumPhotonEt.push_back(IsoR03_sumPhotonEt);
				ftree->Ds_IsoR03_sumPUPt.push_back(IsoR03_sumPUPt);
				ftree->Ds_IsoR03_PFIso.push_back(IsoR03_PFIso);
				ftree->Ds_IsoR03_RelIso.push_back(IsoR03_RelIso);

				ftree->Ds_IsoR04_sumChargedHadronPt.push_back(IsoR04_sumChargedHadronPt);
				ftree->Ds_IsoR04_sumNeutralHadronEt.push_back(IsoR04_sumNeutralHadronEt);
				ftree->Ds_IsoR04_sumPhotonEt.push_back(IsoR04_sumPhotonEt);
				ftree->Ds_IsoR04_sumPUPt.push_back(IsoR04_sumPUPt);
				ftree->Ds_IsoR04_PFIso.push_back(IsoR04_PFIso);
				ftree->Ds_IsoR04_RelIso.push_back(IsoR04_RelIso);

				ftree->Kp_match.push_back(Kp_match);
				ftree->Km_match.push_back(Km_match);
				ftree->pi_match.push_back(pi_match);
			}
		}
	}

	for (Int_t idx_mu = 0; idx_mu < static_cast<Int_t>(MuonHandle->size()); idx_mu++)
	{
		const pat::Muon &muon = (*MuonHandle)[idx_mu];
		if (muon.charge() != -1)
			continue;
		if (muon.pt() < 0.5)
			continue;
		if (muon.p() < 1)
			continue;
		if (!muon.isTightMuon(primvtx))
			continue;

		Bool_t mu_match = (idx_mu == mu_match_idx);

		std::pair<Bool_t, Measurement1D> mu_ipresult = IPTools::absoluteImpactParameter3D((*ttBuilder).build(muon.muonBestTrack()), primvtx);

		ftree->mu_charge.push_back(muon.charge());
		ftree->mu_eta.push_back(muon.eta());
		ftree->mu_phi.push_back(muon.phi());
		ftree->mu_px.push_back(muon.px());
		ftree->mu_py.push_back(muon.py());
		ftree->mu_pz.push_back(muon.pz());
		ftree->mu_pt.push_back(muon.pt());
		ftree->mu_p.push_back(muon.p());
		ftree->mu_isHighPt.push_back(muon.isHighPtMuon(primvtx));
		ftree->mu_isLoose.push_back(muon.isLooseMuon());
		ftree->mu_isMedium.push_back(muon.isMediumMuon());
		ftree->mu_isSoft.push_back(muon.isSoftMuon(primvtx));
		ftree->mu_isTight.push_back(muon.isTightMuon(primvtx));
		ftree->mu_isPF.push_back(muon.isPFMuon());
		ftree->mu_isTracker.push_back(muon.isTrackerMuon());
		ftree->mu_isGlobal.push_back(muon.isGlobalMuon());

		if (mu_ipresult.first)
		{
			ftree->mu_ip.push_back(mu_ipresult.second.value());
			ftree->mu_iperr.push_back(mu_ipresult.second.error());
			ftree->mu_ipchi2.push_back(mu_ipresult.second.significance());
		}
		else
		{
			ftree->mu_ip.push_back(null);
			ftree->mu_iperr.push_back(null);
			ftree->mu_ipchi2.push_back(null);
		}

		ftree->mu_IsoR03_sumChargedHadronPt.push_back(muon.pfIsolationR03().sumChargedHadronPt);
		ftree->mu_IsoR03_sumNeutralHadronEt.push_back(muon.pfIsolationR03().sumNeutralHadronEt);
		ftree->mu_IsoR03_sumPhotonEt.push_back(muon.pfIsolationR03().sumPhotonEt);
		ftree->mu_IsoR03_sumPUPt.push_back(muon.pfIsolationR03().sumPUPt);
		ftree->mu_IsoR03_PFIso.push_back(muon.pfIsolationR03().sumChargedHadronPt + std::max(0.0, muon.pfIsolationR03().sumNeutralHadronEt + muon.pfIsolationR03().sumPhotonEt - 0.5 * muon.pfIsolationR03().sumPUPt));
		ftree->mu_IsoR03_RelIso.push_back((muon.pfIsolationR03().sumChargedHadronPt + std::max(0.0, muon.pfIsolationR03().sumNeutralHadronEt + muon.pfIsolationR03().sumPhotonEt - 0.5 * muon.pfIsolationR03().sumPUPt)) / muon.pt());

		ftree->mu_IsoR04_sumChargedHadronPt.push_back(muon.pfIsolationR04().sumChargedHadronPt);
		ftree->mu_IsoR04_sumNeutralHadronEt.push_back(muon.pfIsolationR04().sumNeutralHadronEt);
		ftree->mu_IsoR04_sumPhotonEt.push_back(muon.pfIsolationR04().sumPhotonEt);
		ftree->mu_IsoR04_sumPUPt.push_back(muon.pfIsolationR04().sumPUPt);
		ftree->mu_IsoR04_PFIso.push_back(muon.pfIsolationR04().sumChargedHadronPt + std::max(0.0, muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt));
		ftree->mu_IsoR04_RelIso.push_back((muon.pfIsolationR04().sumChargedHadronPt + std::max(0.0, muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt());

		ftree->mu_match.push_back(mu_match);
	}

	ftree->tree->Fill();

	// const edm::TriggerNames &names = iEvent.triggerNames(*triggerHandle);

	// for (Size_t i = 0; i < names.size(); i++)
	// {
	// 	TString trigName = names.triggerName(i);
	// 	if (!trigName.Contains("HLT"))
	// 		continue;
	// 	if (!trigName.Contains("Mu") && !trigName.Contains("mu") && !trigName.Contains("MET"))
	// 		continue;
	// 	Bool_t pass = triggerHandle->accept(i);
	// 	triggerTotal_[trigName]++;
	// 	if (pass)
	// 		triggerCounts_[trigName]++;
	// }
}

Bool_t RecoCompare::hasAncestor(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const
{
	if (gp.numberOfMothers() == 0)
		return false;

	for (size_t i = 0; i < gp.numberOfMothers(); i++)
	{
		const reco::Candidate *mother = gp.mother(i);
		if (!mother)
			continue;

		if (mother->pdgId() == motherid && mother->numberOfDaughters() >= 2)
		{
			Bool_t daughter1_flag = false;
			Bool_t daughter2_flag = false;
			for (size_t j = 0; j < mother->numberOfDaughters(); j++)
			{
				if (mother->daughter(j)->pdgId() == daughter1id)
					daughter1_flag = true;
				if (mother->daughter(j)->pdgId() == daughter2id)
					daughter2_flag = true;
			}

			if (daughter1_flag && daughter2_flag)
				return true;
		}

		if (hasAncestor(*(gp.motherRef(i)), motherid, daughter1id, daughter2id))
			return true;
	}

	return false;
}

void RecoCompare::endJob()
{
	// std::cout << "Total: " << nEventsProcessed_ << std::endl;
	// std::cout << "Kp match: " << nEventsKpMatch_ << ", nonmatch: " << nEventsKpNonMatch_ << std::endl;
	// std::cout << "Km match: " << nEventsKmMatch_ << ", nonmatch: " << nEventsKmNonMatch_ << std::endl;
	// std::cout << "pi match: " << nEventspiMatch_ << ", nonmatch: " << nEventspiNonMatch_ << std::endl;
	// std::cout << "phi match: " << nEventsphiMatch_ << ", nonmatch: " << nEventsphiNonMatch_ << std::endl;
	// std::cout << "Ds match: " << nEventsDsMatch_ << ", nonmatch: " << nEventsDsNonMatch_ << std::endl;
	// std::cout << "Mu match: " << nEventsMuMatch_ << ", nonmatch: " << nEventsMuNonMatch_ << std::endl;
	// std::cout << "H match: " << nEventsMatch_ << ", nonmatch: " << nEventsNonMatch_ << std::endl;

	// edm::LogPrint("RecoCompare") << "\n========== Trigger Statistics ==========";
	// edm::LogPrint("RecoCompare") << std::setw(60) << "Trigger Name"
	// 							  << std::setw(15) << "Fired"
	// 							  << std::setw(15) << "Total"
	// 							  << std::setw(15) << "Efficiency";
	// edm::LogPrint("RecoCompare") << std::string(105, '-');

	// std::vector<std::pair<TString, double>> triggerEfficiencies;

	// for (const auto &entry : triggerCounts_)
	// {
	// 	const TString trigName = entry.first;
	// 	Int_t fired = entry.second;
	// 	Int_t total = triggerTotal_[trigName];
	// 	double efficiency = (total > 0) ? (double)fired / total * 100.0 : 0.0;

	// 	triggerEfficiencies.push_back({trigName, efficiency});
	// }

	// std::sort(triggerEfficiencies.begin(), triggerEfficiencies.end(),
	// 		  [](const auto &a, const auto &b)
	// 		  { return a.second > b.second; });

	// for (const auto &item : triggerEfficiencies)
	// {
	// 	const TString &trigName = item.first;
	// 	Int_t fired = triggerCounts_[trigName];
	// 	Int_t total = triggerTotal_[trigName];
	// 	double efficiency = item.second;

	// 	edm::LogPrint("RecoCompare") << std::setw(80) << trigName
	// 								  << std::setw(15) << fired
	// 								  << std::setw(15) << total
	// 								  << std::setw(14) << std::fixed << std::setprecision(2)
	// 								  << efficiency << "%";
	// }

	// edm::LogPrint("RecoCompare") << "========================================\n";
}

void RecoCompare::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RecoCompare);
