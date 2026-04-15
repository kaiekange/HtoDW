#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>
#include <string>
#include <optional>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
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

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

class GenReco : public edm::one::EDAnalyzer<>
{
public:
	explicit GenReco(const edm::ParameterSet &);
	~GenReco() override {};

	static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
	void beginJob() override;
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	void endJob() override;

	edm::EDGetTokenT<reco::GenParticleCollection> GenPartToken_;
};

GenReco::GenReco(const edm::ParameterSet &iConfig)
{
	edm::InputTag GenPartTag_("prunedGenParticles");
	GenPartToken_ = consumes<reco::GenParticleCollection>(GenPartTag_);
}

void GenReco::beginJob()
{}

void GenReco::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	edm::Handle<reco::GenParticleCollection> GenPartHandle;
	iEvent.getByToken(GenPartToken_, GenPartHandle);

	std::cout << "hard process" << std::endl;
	for (Int_t i = 0; i < static_cast<Int_t>(GenPartHandle->size()); i++)
	{
		const reco::GenParticle &gp = (*GenPartHandle)[i];

		if (!gp.isHardProcess())
			continue;
		if ((abs(gp.pdgId()) != 431))
			continue;

		const reco::Candidate* current = &gp;
		while (true)
		{
			bool foundSame = false;
			const reco::Candidate* d = current->daughter(0);
			if (d && d->pdgId() == current->pdgId())
			{
				current = d;
				foundSame = true;
			}
			if (!foundSame) break;
		}

		std::cout << "pt: " << current->pt() << "; daughters: ";
		for (size_t i = 0; i < current->numberOfDaughters(); ++i)
			std::cout << current->daughter(i)->pdgId() << ", ";
		std::cout << "end" << std::endl;
	}

	std::cout << "random process" << std::endl;
	for (Int_t i = 0; i < static_cast<Int_t>(GenPartHandle->size()); i++)
	{
		const reco::GenParticle &gp = (*GenPartHandle)[i];

		if (gp.isHardProcess())
			continue;
		if ((abs(gp.pdgId()) != 431))
			continue;

		const reco::Candidate* current = &gp;
		while (true)
		{
			bool foundSame = false;
			const reco::Candidate* d = current->daughter(0);
			if (d && d->pdgId() == current->pdgId())
			{
				current = d;
				foundSame = true;
			}
			if (!foundSame) break;
		}

		std::cout << "pt: " << current->pt() << "; daughters: ";
		for (size_t i = 0; i < current->numberOfDaughters(); ++i)
			std::cout << current->daughter(i)->pdgId() << ", ";
		std::cout << "end" << std::endl;
	}
}

void GenReco::endJob()
{}

void GenReco::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(GenReco);
