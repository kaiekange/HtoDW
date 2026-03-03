
Bool_t FullReco::JMECorrFilter(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	const Bool_t usePuppi = JecApplication::usesPuppiMet(year_);

	// Inputs
	edm::Handle<pat::JetCollection> srcJetsHandle;
	iEvent.getByToken(srcJetsToken_, srcJetsHandle);

	edm::Handle<Double_t> rhoHandle;
	iEvent.getByToken(rhoToken_, rhoHandle);
	const Double_t rho = rhoHandle.isValid() ? *rhoHandle : 0.0;

	edm::Handle<pat::METCollection> metHandle;
	if (usePuppi)
		iEvent.getByToken(metPuppiToken_, metHandle);
	else
		iEvent.getByToken(metPFToken_, metHandle);

	if (!metHandle.isValid() || metHandle->empty())
		throw edm::Exception(edm::errors::ProductNotFound) << "JecAppliedMetProducer: MET collection is invalid or empty (" << (usePuppi ? srcMetPuppiTag_.encode() : srcMetPFTag_.encode()) << ").";

	auto CorrectedJets = std::make_unique<pat::JetCollection>();
	CorrectedJets->reserve(srcJetsHandle->size());

	// Build systematics payload
	JecApplication::SystematicOptions systOpts{};

	// Per-event Applier (handles JEC/JER handles and run-based residual logic)
	JecApplication::Applier applier = [&]() -> JecApplication::Applier
	{
		if (isData_)
		{
			std::optional<Double_t> runNumber;
			if (JecApplication::requiresRunBasedResidual(year_))
				runNumber = static_cast<Double_t>(iEvent.id().run());
			return JecApplication::Applier::DataAK4(*jecCfg_, year_, era_, isDebug_, runNumber);
		}
		else
			return JecApplication::Applier::McAK4(*jecCfg_, year_, isDebug_);
	}();

	const pat::MET &metIn = metHandle->front();

	const Double_t metRawPt = getMetRawPt(metIn);
	const Double_t metRawPhi = getMetRawPhi(metIn);

	// Build jetsForMet (RAW AK4 BEFORE any modifications)
	std::vector<JecApplication::JetForMet> jetsForMet;
	std::vector<JecApplication::JerInputs> jersForMet;
	jetsForMet.reserve(srcJetsHandle->size());
	jersForMet.reserve(srcJetsHandle->size());

	// Loop jets
	for (const auto &jet : *srcJetsHandle)
	{
		pat::Jet outJet = jet;

		const Double_t rawFactor = computeRawFactorFromMiniAOD(jet);

		JecApplication::JesInputs jesIn{jet.pt(), jet.eta(), jet.phi(), jet.jetArea(), rho, rawFactor};

		// --- JES nominal factor (via Applier) ---
		const Double_t jesFactor = applier.jesFactorNominal(jesIn);

		// Apply JES nominal
		setP4Scaled(outJet, jesFactor);

		// Gen-matching if available
		JecApplication::JerInputs jerIn{};
		jerIn.event = iEvent.id().event();
		jerIn.rho = rho;

		if (const auto *gj = jet.genJet())
		{
			jerIn.hasGen = true;
			jerIn.genPt = gj->pt();
			jerIn.genEta = gj->eta();
			jerIn.genPhi = gj->phi();
			jerIn.maxDr = 0.2; // AK4
		}

		// --- JER (MC only) ---
		double jerFactor = 1.0;
		if (!isData_)
		{
			// State after JES (and optional JES syst)
			JecApplication::JesInputs jAfter{outJet.pt(), outJet.eta(), outJet.phi(), jet.jetArea(), rho, 0.0};

			jerFactor = applier.jerFactor(jAfter, jerIn, systOpts);
			setP4Scaled(outJet, jerFactor);
		}

		// Attach userFloats for downstream bookkeeping
		double jesSystFactor = 1.0;
		const Double_t total = jesFactor * jesSystFactor * jerFactor;
		addUserFloats(outJet, jesFactor, jesSystFactor, jerFactor, total);

		CorrectedJets->push_back(std::move(outJet));

		JecApplication::JetForMet v;
		v.phi = jet.phi();
		v.eta = jet.eta();
		v.area = jet.jetArea();
		v.rawPt = getRawPtFromMiniAOD(jet);
		v.muonSubtrFactor = getMuonSubtrFraction(jet, muonSubtrUserFloat_);
		v.chEmEf = jet.chargedEmEnergyFraction();
		v.neEmEf = jet.neutralEmEnergyFraction();
		jetsForMet.push_back(v);

		jersForMet.push_back(jerIn);
	}

	// Call the shared Type-1 MET routine through Applier
	JecApplication::MetInputs metInputs{metRawPt, metRawPhi};
	const TLorentzVector p4CorrectedMET = applier.correctedMet(metInputs, jetsForMet, jersForMet, rho, systOpts);

	auto CorrectedMet = std::make_unique<pat::METCollection>();
	pat::MET metOut = makeCorrectedMetFrom(metIn, p4CorrectedMET);
	metOut.addUserFloat("JERC:MET_raw_pt_used", metRawPt);
	metOut.addUserFloat("JERC:MET_raw_phi_used", metRawPhi);
	CorrectedMet->push_back(std::move(metOut));

	iEvent.put(std::move(CorrectedJets), "CorrectedAK4");
	iEvent.put(std::move(CorrectedMet), "CorrectedMet");

	Bool_t passEvent = true;
	try
	{
		if (srcJetsHandle.isValid())
		{
			for (const auto &jIt : *srcJetsHandle)
			{
				if (jvmChecker_.checkJetInVetoRegion(jIt.eta(), jIt.phi(), jIt.pt(), 6, jIt.chargedEmEnergyFraction(), jIt.neutralEmEnergyFraction()))
				{
					passEvent = false;
					break; // No need to loop over the remaining jets.
				}
			}
		}
	}
	catch (std::exception const &e)
	{
		std::cout << "[JvmAppliedEventFilter] exception in filter -> " << e.what() << std::endl;
		// Be conservative on exception: do not veto the event.
		return true;
	}
	return passEvent;
}