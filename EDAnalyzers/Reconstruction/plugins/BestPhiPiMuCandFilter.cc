// Keep one Ds(*) candidate per event -- the highest-Ds_pt one -- and drop events with none.
//
// HDsPhiPiMuProducer emits every (Ds, gamma) combination that passes all cuts as a
// pat::CompositeCandidate. The old tuplizer wrote one TTree row per combination and then
// highest_pt_dedup.py::dedup_best_indices collapsed to one row per event: highest Ds_pt,
// ties broken on highest gamma_pt, remaining ties on loop order (np.argmax -> first). This
// filter reproduces that exactly -- the multiple-gamma tie-break matters when the winning Ds
// has more than one passing photon. Events with an empty collection return false (never
// reached tree->Fill() in the old analyzer either) so NanoAODOutputModule drops them.
//
// Ds_pt is bit-identical across every candidate of the same Ds (same fit, copied from the
// same ftree->Ds_pt), so the == tie test below is exact, matching the dedup's pt == pt.max().
//
// The single-element output is what the SimpleCandidateFlatTableProducer (singleton=True,
// which does a bare assert(size()==1)) is pointed at.

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

class BestPhiPiMuCandFilter : public edm::stream::EDFilter<>
{
public:
	explicit BestPhiPiMuCandFilter(const edm::ParameterSet &iConfig)
		: token_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("src")))
	{
		produces<pat::CompositeCandidateCollection>();
	}

	static void fillDescriptions(edm::ConfigurationDescriptions &descriptions)
	{
		edm::ParameterSetDescription desc;
		desc.add<edm::InputTag>("src", edm::InputTag("phipimuCand"));
		descriptions.addDefault(desc);
	}

private:
	bool filter(edm::Event &iEvent, const edm::EventSetup &) override
	{
		edm::Handle<pat::CompositeCandidateCollection> in;
		iEvent.getByToken(token_, in);

		auto out = std::make_unique<pat::CompositeCandidateCollection>();
		if (in->empty())
		{
			iEvent.put(std::move(out));
			return false;
		}

		size_t best = 0;
		for (size_t i = 1; i < in->size(); i++)
		{
			const float dp_i = (*in)[i].userFloat("Ds_pt");
			const float dp_b = (*in)[best].userFloat("Ds_pt");
			if (dp_i > dp_b) { best = i; continue; }
			// same Ds (bit-identical Ds_pt): keep the highest-gamma_pt row, matching
			// dedup_best_indices' np.argmax(gamma_pt[top]). Strict '>' leaves loop order
			// as the final tie-break (first row wins), matching np.argmax.
			if (dp_i == dp_b && (*in)[i].userFloat("gamma_pt") > (*in)[best].userFloat("gamma_pt")) best = i;
		}

		out->push_back((*in)[best]);
		iEvent.put(std::move(out));
		return true;
	}

	const edm::EDGetTokenT<pat::CompositeCandidateCollection> token_;
};

DEFINE_FWK_MODULE(BestPhiPiMuCandFilter);
