// Pileup reweighting DENOMINATOR for the custom NanoAOD: the genWeight-weighted nTrueInt
// distribution over ALL processed MC events, written to the Runs tree as a
// nanoaod::MergeableCounterTable (same mechanism as genWeightsTable's genEventSumw).
//
// Why it exists: the MC profile that pileup weights divide by must be taken BEFORE any
// event selection -- selection shifts pileup-correlated quantities (measured -19% to
// +4.5% on a PU proxy, see memory feedback_pu_reweight_denominator). The old tuplizer
// filled such a TH1 (nTrueInt_all) as the first line of analyze(); the NanoAOD migration
// dropped it, and nano's Events tree only holds events passing phipimu_step, so the
// denominator could no longer be recovered. This module restores it.
//
// Scheduling: must run on every event. It is placed in the unfiltered, scheduled sumw_step
// next to genWeightsTable (nano_phipimu_cff.customise), NOT behind phipimu_step's filter.
// Closure: sum(nTrueIntSumw) + nTrueIntSumwOverflow + nTrueIntSumwNoInfo == genEventSumw
// for the same file, since both are sum(genWeight) over all processed events.
//
// Binning: 100 unit bins, bin i = [i, i+1), i = 0..99, matching the old nTrueInt_all TH1
// (0-100, one bin per unit) and the unit-width edges of LUM puWeights.json (which covers
// [0,99] with flow=clamp). Values >= 100 go to a separate overflow column rather than being
// clamped here, and events without an in-time (bunch crossing 0) PileupSummaryInfo go to a
// separate no-info column, so nothing is folded in silently; clamp downstream if needed.
// The output label must end in "Table" to match NANOAODSIMEventContent's keep rule
// "keep nanoaodMergeableCounterTable_*Table_*_*".

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/NanoAOD/interface/MergeableCounterTable.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include <cmath>
#include <memory>
#include <vector>

namespace {
  constexpr unsigned kNBins = 100;

  struct NTrueIntCounter {
    std::vector<double> sumw = std::vector<double>(kNBins, 0.0);
    double sumwOverflow = 0.0;
    double sumwNoInfo = 0.0;
    long long count = 0;

    void clear() {
      std::fill(sumw.begin(), sumw.end(), 0.0);
      sumwOverflow = sumwNoInfo = 0.0;
      count = 0;
    }
    void merge(const NTrueIntCounter& other) {
      for (unsigned i = 0; i < kNBins; ++i)
        sumw[i] += other.sumw[i];
      sumwOverflow += other.sumwOverflow;
      sumwNoInfo += other.sumwNoInfo;
      count += other.count;
    }
  };
}  // namespace

class NTrueIntCounterTableProducer : public edm::global::EDProducer<edm::StreamCache<NTrueIntCounter>,
                                                                    edm::RunSummaryCache<NTrueIntCounter>,
                                                                    edm::EndRunProducer> {
public:
  explicit NTrueIntCounterTableProducer(const edm::ParameterSet& params)
      : genToken_(consumes<GenEventInfoProduct>(params.getParameter<edm::InputTag>("genEvent"))),
        puToken_(consumes<std::vector<PileupSummaryInfo>>(params.getParameter<edm::InputTag>("pileupInfo"))) {
    produces<nanoaod::MergeableCounterTable, edm::Transition::EndRun>();
  }

  void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup&) const override {
    NTrueIntCounter* counter = streamCache(id);
    const double w = iEvent.get(genToken_).weight();
    ++counter->count;

    for (const PileupSummaryInfo& pu : iEvent.get(puToken_)) {
      if (pu.getBunchCrossing() != 0)
        continue;  // in-time crossing only, same as the producers' per-event nTrueInt
      const float nTrueInt = pu.getTrueNumInteractions();
      if (nTrueInt >= 0 && nTrueInt < kNBins)
        counter->sumw[static_cast<unsigned>(std::floor(nTrueInt))] += w;
      else
        counter->sumwOverflow += w;  // >= 100 (or a negative value, which should not occur)
      return;
    }
    counter->sumwNoInfo += w;
  }

  std::unique_ptr<NTrueIntCounter> beginStream(edm::StreamID) const override {
    return std::make_unique<NTrueIntCounter>();
  }
  void streamBeginRun(edm::StreamID id, const edm::Run&, const edm::EventSetup&) const override {
    streamCache(id)->clear();
  }
  std::shared_ptr<NTrueIntCounter> globalBeginRunSummary(const edm::Run&, const edm::EventSetup&) const override {
    return std::make_shared<NTrueIntCounter>();
  }
  void streamEndRunSummary(edm::StreamID id,
                           const edm::Run&,
                           const edm::EventSetup&,
                           NTrueIntCounter* runCounter) const override {
    runCounter->merge(*streamCache(id));
  }
  void globalEndRunSummary(const edm::Run&, const edm::EventSetup&, NTrueIntCounter*) const override {}

  void globalEndRunProduce(edm::Run& iRun, const edm::EventSetup&, const NTrueIntCounter* runCounter) const override {
    auto out = std::make_unique<nanoaod::MergeableCounterTable>();
    out->addVFloat("nTrueIntSumw",
                   "sum of genWeight per unit bin of Pileup_nTrueInt, bin i = [i,i+1), i=0..99, over ALL "
                   "processed events (pre-selection pileup reweighting denominator)",
                   runCounter->sumw);
    out->addFloat("nTrueIntSumwOverflow", "sum of genWeight for Pileup_nTrueInt >= 100", runCounter->sumwOverflow);
    out->addFloat("nTrueIntSumwNoInfo",
                  "sum of genWeight for events without an in-time PileupSummaryInfo",
                  runCounter->sumwNoInfo);
    out->addInt("nTrueIntEventCount", "number of events entering the nTrueInt counter", runCounter->count);
    iRun.put(std::move(out));
  }

private:
  const edm::EDGetTokenT<GenEventInfoProduct> genToken_;
  const edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;
};

DEFINE_FWK_MODULE(NTrueIntCounterTableProducer);
