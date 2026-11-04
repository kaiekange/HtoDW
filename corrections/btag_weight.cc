// Fixed-WP b-tag event weight, C++ port of compute_btag_event_weight_multi_wp from
// gitlab.cern.ch/cms-analysis-corrections/BTV/examples lib/btag_sf.py
// (commit 9b42445d9fb2a14eb36486f5c32fbe5d6c8c1ec4, 2026-04-17). BTV provides no C++ version of
// this function; the Python one needs Python >= 3.10 (CMSSW_15_0_5 has 3.9).
// Line-by-line: same N+1 jet categories, same heavy/light split (hadronFlavour > 3 / < 4),
// same checks (p_MC, p_DATA >= 0 and weight >= 0 unless ignore_warnings), same replacement of a
// non-finite event weight by 1.
// Compiled through ROOT (corrections/btv.py) after correctionlib.register_pyroot_binding().
#include "correction.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "ROOT/RVec.hxx"

class BtagWeight {
 public:
  // algo_key: correctionlib key of the tagger ("particleNet" -> particleNet_comb/_light/_wp_values)
  // eff_key : correction name in the efficiency-map file ("particleNet_efficiency")
  // wps     : WP names; re-sorted by threshold, loosest first (as in btag_sf.py)
  BtagWeight(const std::string& sf_json, const std::string& eff_json, const std::string& algo_key,
             const std::string& eff_key, const std::vector<std::string>& wps)
      : sf_set_(correction::CorrectionSet::from_file(sf_json)),
        eff_set_(correction::CorrectionSet::from_file(eff_json)) {
    sf_heavy_ = sf_set_->at(algo_key + "_comb");
    sf_light_ = sf_set_->at(algo_key + "_light");
    eff_ = eff_set_->at(eff_key);
    auto wp_values = sf_set_->at(algo_key + "_wp_values");
    std::vector<std::pair<double, std::string>> v;
    for (const auto& wp : wps) v.emplace_back(wp_values->evaluate({wp}), wp);
    std::stable_sort(v.begin(), v.end(),
                     [](const auto& a, const auto& b) { return a.first < b.first; });
    for (const auto& x : v) {
      thr_.push_back(x.first);
      wps_.push_back(x.second);
    }
  }

  // Flat jet arrays of all events, plus the number of jets per event. Returns one weight per event.
  ROOT::RVec<double> weights(const ROOT::RVec<double>& pt, const ROOT::RVec<double>& eta,
                             const ROOT::RVec<int>& flav, const ROOT::RVec<double>& disc,
                             const ROOT::RVec<int>& counts, const std::string& variation_heavy,
                             const std::string& variation_light, bool ignore_warnings) const {
    const size_t nwp = wps_.size();
    ROOT::RVec<double> out(counts.size());
    size_t j = 0;
    bool neg_mc = false, neg_data = false;
    for (size_t e = 0; e < counts.size(); ++e) {
      double p_mc = 1.0, p_data = 1.0;
      for (int c = 0; c < counts[e]; ++c, ++j) {
        const double abseta = std::abs(eta[j]);
        std::vector<double> eff_mc(nwp), eff_data(nwp);
        for (size_t k = 0; k < nwp; ++k) {
          eff_mc[k] = eff_->evaluate({wps_[k], flav[j], abseta, pt[j]});
          const double sf = flav[j] > 3
              ? sf_heavy_->evaluate({variation_heavy, wps_[k], flav[j], abseta, pt[j]})
              : sf_light_->evaluate({variation_light, wps_[k], flav[j], abseta, pt[j]});
          eff_data[k] = sf * eff_mc[k];
        }
        // category: index of the tightest WP passed, -1 if the loosest fails
        int cat = -1;
        for (size_t k = 0; k < nwp; ++k)
          if (disc[j] > thr_[k]) cat = static_cast<int>(k);
        double f_mc, f_data;
        if (cat == static_cast<int>(nwp) - 1) {          // passes tightest
          f_mc = eff_mc[nwp - 1];
          f_data = eff_data[nwp - 1];
        } else if (cat < 0) {                             // fails loosest
          f_mc = 1 - eff_mc[0];
          f_data = 1 - eff_data[0];
        } else {                                          // passes wp[cat], fails wp[cat+1]
          f_mc = eff_mc[cat] - eff_mc[cat + 1];
          f_data = eff_data[cat] - eff_data[cat + 1];
        }
        p_mc *= f_mc;
        p_data *= f_data;
      }
      if (!(p_mc >= 0)) neg_mc = true;      // as ak.all(p_MC >= 0)
      if (!(p_data >= 0)) neg_data = true;
      out[e] = p_data / p_mc;
    }
    if (j != pt.size()) throw std::runtime_error("BtagWeight: sum(counts) != number of jets");
    if (!ignore_warnings && neg_mc)
      throw std::runtime_error("Negative MC efficiency product (heavy " + variation_heavy +
                               ", light " + variation_light + ")");
    if (!ignore_warnings && neg_data)
      throw std::runtime_error("Negative DATA efficiency product (heavy " + variation_heavy +
                               ", light " + variation_light + ")");
    for (auto& w : out) {
      if (!std::isfinite(w)) w = 1;   // replace inf or NaN with 1
      if (!ignore_warnings && w < 0)
        throw std::runtime_error("Negative btag SF (heavy " + variation_heavy + ", light " +
                                 variation_light + ")");
    }
    return out;
  }

 private:
  std::unique_ptr<const correction::CorrectionSet> sf_set_, eff_set_;
  correction::Correction::Ref sf_heavy_, sf_light_, eff_;
  std::vector<std::string> wps_;
  std::vector<double> thr_;
};
