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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TTree.h"

class TrigAnalyzer : public edm::one::EDAnalyzer<>
{
public:
	explicit TrigAnalyzer(const edm::ParameterSet &);
	~TrigAnalyzer() override {};

	static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
	void beginJob() override;
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	void endJob() override;

	Bool_t ThreeBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id, const Int_t daughter3id) const;
	Bool_t TwoBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const;

	edm::EDGetTokenT<reco::GenParticleCollection> GenPartToken_;
	edm::EDGetTokenT<edm::TriggerResults> triggerToken_;

	unsigned long nEventsProcessed_;
	unsigned long nEventsLepton_;
	unsigned long nEventsGen_;
	unsigned long nEventsTriggered_;

	std::map<TString, int> triggerCounts_;
	std::map<TString, int> triggerTotal_;
};

TrigAnalyzer::TrigAnalyzer(const edm::ParameterSet &iConfig)
{
	edm::InputTag GenPartTag_("prunedGenParticles");
	GenPartToken_ = consumes<reco::GenParticleCollection>(GenPartTag_);

	edm::InputTag TriggerBitsTag_("TriggerResults", "", "HLT");
	triggerToken_ = consumes<edm::TriggerResults>(TriggerBitsTag_);
}

void TrigAnalyzer::beginJob()
{
	nEventsProcessed_ = 0;
	nEventsLepton_ = 0;
	nEventsGen_ = 0;
}

void TrigAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	nEventsProcessed_++;

	Int_t Gen_H_idx = 0;
	Int_t Gen_Dsst_idx = 0;
	Int_t Gen_W_idx = 0;

	Int_t num_Gen_H = 0;
	Int_t num_Gen_Dsst = 0;
	Int_t num_Gen_W = 0;
	Int_t num_Gen_mu = 0;
	Int_t num_Gen_e = 0;
	Int_t num_Gen_Ds = 0;
	Int_t num_Gen_gamma = 0;
	Int_t num_Gen_Kp = 0;
	Int_t num_Gen_Km = 0;
	Int_t num_Gen_pi = 0;

	//==================================================//
	//--------------------------------------------------//
	// Gen
	//--------------------------------------------------//
	//==================================================//
	edm::Handle<reco::GenParticleCollection> GenPartHandle;
	iEvent.getByToken(GenPartToken_, GenPartHandle);

	for (Int_t i = 0; i < static_cast<Int_t>(GenPartHandle->size()); i++)
	{
		const reco::GenParticle &gp = (*GenPartHandle)[i];
		if (!gp.isHardProcess())
			continue;
		if (gp.pdgId() == 25)
		{
			Gen_H_idx = i;
			num_Gen_H++;
		}
		else if ((abs(gp.pdgId()) == 433))
		{
			Gen_Dsst_idx = i;
			num_Gen_Dsst++;
		}
		else if (abs(gp.pdgId()) == 24)
		{
			Gen_W_idx = i;
			num_Gen_W++;
		}
		else if (abs(gp.pdgId()) == 13)
		{
			num_Gen_mu++;
		}
		else if (abs(gp.pdgId()) == 11)
		{
			num_Gen_e++;
		}
	}
	if (!((num_Gen_H == 1) && (num_Gen_Dsst == 1) && (num_Gen_W == 1) && (num_Gen_mu == 1) && (num_Gen_e == 0)))
		return;
	nEventsLepton_++;

	const reco::GenParticle &H_gp = (*GenPartHandle)[Gen_H_idx];
	const reco::GenParticle &Dsst_gp = (*GenPartHandle)[Gen_Dsst_idx];
	const reco::GenParticle &W_gp = (*GenPartHandle)[Gen_W_idx];

	const Int_t Dsst_charge = Dsst_gp.charge();

	for (Int_t i = 0; i < static_cast<Int_t>(GenPartHandle->size()); i++)
	{
		const reco::GenParticle &gp = (*GenPartHandle)[i];
		if ((gp.status() == 2) && gp.pdgId() == (Dsst_charge * 431))
		{
			if (TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()) && (gp.numberOfDaughters() > 2))
			{
				Bool_t daughter1_flag = false;
				Bool_t daughter2_flag = false;
				Bool_t daughter3_flag = false;
				for (size_t j = 0; j < gp.numberOfDaughters(); j++)
				{
					if (gp.daughter(j)->pdgId() == 321)
						daughter1_flag = true;
					if (gp.daughter(j)->pdgId() == -321)
						daughter2_flag = true;
					if (gp.daughter(j)->pdgId() == (Dsst_charge * 211))
						daughter3_flag = true;
				}
				if (daughter1_flag && daughter2_flag && daughter3_flag)
					num_Gen_Ds++;
			}
		}
		else if (gp.status() == 1)
		{
			if (gp.pdgId() == (Dsst_charge * 321))
			{
				if (ThreeBodyDecay(gp, Dsst_charge * 431, 321, -321, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
					num_Gen_Kp++;
			}
			else if (gp.pdgId() == (-Dsst_charge * 321))
			{
				if (ThreeBodyDecay(gp, Dsst_charge * 431, 321, -321, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
					num_Gen_Km++;
			}
			else if (gp.pdgId() == (Dsst_charge * 211))
			{
				if (ThreeBodyDecay(gp, Dsst_charge * 431, 321, -321, Dsst_charge * 211) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
					num_Gen_pi++;
			}
			else if (gp.pdgId() == 22)
			{
				if ((gp.mother(0)->pdgId() == Dsst_gp.pdgId()) && TwoBodyDecay(gp, Dsst_gp.pdgId(), Dsst_charge * 431, 22) && TwoBodyDecay(gp, 25, Dsst_gp.pdgId(), W_gp.pdgId()))
					num_Gen_gamma++;
			}
		}
	}

	if (!((num_Gen_Ds == 1) && (num_Gen_Kp == 1) && (num_Gen_Km == 1) && (num_Gen_pi == 1) && (num_Gen_gamma == 1)))
		return;

	//==================================================//
	//--------------------------------------------------//
	// Triggers
	//--------------------------------------------------//
	//==================================================//
	edm::Handle<edm::TriggerResults> triggerHandle;
	iEvent.getByToken(triggerToken_, triggerHandle);

	const edm::TriggerNames &names = iEvent.triggerNames(*triggerHandle);

	for (Size_t i = 0; i < names.size(); i++)
	{
		TString trigName = names.triggerName(i);
		if (!trigName.Contains("HLT"))
			continue;
		// if (!trigName.Contains("Mu") && !trigName.Contains("mu") && !trigName.Contains("MET"))
		// 	continue;
		Bool_t pass = triggerHandle->accept(i);
		triggerTotal_[trigName]++;
		if (pass)
			triggerCounts_[trigName]++;
	}
}

void TrigAnalyzer::endJob()
{
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

Bool_t TrigAnalyzer::ThreeBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id, const Int_t daughter3id) const
{
	if (gp.numberOfMothers() == 0)
		return false;

	for (size_t i = 0; i < gp.numberOfMothers(); i++)
	{
		const reco::Candidate *mother = gp.mother(i);
		if (!mother)
			continue;

		if (mother->pdgId() == motherid && mother->numberOfDaughters() >= 3)
		{
			Bool_t daughter1_flag = false;
			Bool_t daughter2_flag = false;
			Bool_t daughter3_flag = false;
			for (size_t j = 0; j < mother->numberOfDaughters(); j++)
			{
				if (mother->daughter(j)->pdgId() == daughter1id)
					daughter1_flag = true;
				if (mother->daughter(j)->pdgId() == daughter2id)
					daughter2_flag = true;
				if (mother->daughter(j)->pdgId() == daughter3id)
					daughter3_flag = true;
			}

			if (daughter1_flag && daughter2_flag && daughter3_flag)
				return true;
		}

		if (ThreeBodyDecay(*(gp.motherRef(i)), motherid, daughter1id, daughter2id, daughter3id))
			return true;
	}

	return false;
}

Bool_t TrigAnalyzer::TwoBodyDecay(const reco::GenParticle &gp, const Int_t motherid, const Int_t daughter1id, const Int_t daughter2id) const
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

		if (TwoBodyDecay(*(gp.motherRef(i)), motherid, daughter1id, daughter2id))
			return true;
	}

	return false;
}

void TrigAnalyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(TrigAnalyzer);
