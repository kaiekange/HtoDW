// Combinatorial background yield: WJets run with doTruth=False, so every
// candidate is a fake by construction. The analysis keeps ONE candidate per
// event (highest-pT Ds, then highest-pT gamma within it), so the background
// yield is simply the number of events that yield a selected candidate.
//
//   PF only : events with >=1 candidate whose legs are all packedPFCandidates
//   PF + LT : events with >=1 candidate from the merged pool
//
// n_lost==0 candidates are identical to a PF-only run (verified on signal), so
// one run gives both numbers.
//
//   root -l -b -q 'bkg.C("/eos/home-k/kakang/H2DW/Analysis/MC/2022/WJets_LostTrackReco_phipimu")'

void bkg(const char* dir, const char* proc="LTReco"){
  TChain t(Form("%s/Events",proc));
  int nf = t.Add(Form("%s/output_*.root",dir));

  ULong64_t ev; UInt_t run,lumi; Int_t cand,nlost; Float_t dspt,gpt; Bool_t hasg;
  t.SetBranchAddress("event_id",&ev); t.SetBranchAddress("run_id",&run); t.SetBranchAddress("lumi_id",&lumi);
  t.SetBranchAddress("cand_id",&cand); t.SetBranchAddress("n_lost",&nlost);
  t.SetBranchAddress("Ds_pt",&dspt);
  t.SetBranchAddress("gamma_pt",&gpt);  t.SetBranchAddress("has_gamma",&hasg);

  ULong64_t nEvtPF=0, nEvtALL=0;              // events with a selected candidate
  ULong64_t nEvtPF_g=0, nEvtALL_g=0;          // ... and the selected one has a gamma
  double nCandPF=0, nCandALL=0;               // deduped candidate multiplicity

  // per event: best (Ds_pt, then gamma_pt)
  double bPF=-1, bALL=-1, bPFg=-1, bALLg=-1;
  bool gotPF=false, gotALL=false, bPFhasG=false, bALLhasG=false;
  std::set<int> cPF, cALL;
  ULong64_t pev=0; UInt_t prun=0, plumi=0; bool first=true;

  auto flush=[&](){
    if(!gotPF && !gotALL) return;
    nCandPF += cPF.size(); nCandALL += cALL.size();
    if(gotPF){  nEvtPF++;  if(bPFhasG)  nEvtPF_g++;  }
    if(gotALL){ nEvtALL++; if(bALLhasG) nEvtALL_g++; }
    bPF=bALL=bPFg=bALLg=-1; gotPF=gotALL=false; bPFhasG=bALLhasG=false;
    cPF.clear(); cALL.clear();
  };

  Long64_t N=t.GetEntries();
  for(Long64_t i=0;i<N;i++){
    t.GetEntry(i);
    if(first||ev!=pev||run!=prun||lumi!=plumi){ if(!first) flush(); pev=ev;prun=run;plumi=lumi;first=false; }
    double g = hasg ? gpt : -1;
    cALL.insert(cand);
    if(dspt>bALL || (dspt==bALL && g>bALLg)){ bALL=dspt; bALLg=g; bALLhasG=hasg; }
    gotALL=true;
    if(nlost==0){
      cPF.insert(cand);
      if(dspt>bPF || (dspt==bPF && g>bPFg)){ bPF=dspt; bPFg=g; bPFhasG=hasg; }
      gotPF=true;
    }
  }
  flush();

  printf("\n=== WJets combinatorial background (%d files, %lld rows) ===\n",nf,N);
  printf("\nselected events (one best candidate per event):\n");
  printf("  PF only : %8llu\n",nEvtPF);
  printf("  PF + LT : %8llu    %+lld   (%+.2f%%)\n",nEvtALL,(long long)nEvtALL-(long long)nEvtPF,
         nEvtPF? 100.0*((double)nEvtALL-nEvtPF)/nEvtPF : 0);
  printf("  stat err on the increase: +-%.2f%%\n",
         nEvtPF? 100.0*sqrt((double)(nEvtALL-nEvtPF))/nEvtPF : 0);

  printf("\nselected events whose best candidate has a gamma (Ds* leg):\n");
  printf("  PF only : %8llu\n",nEvtPF_g);
  printf("  PF + LT : %8llu    %+lld   (%+.2f%%)\n",nEvtALL_g,(long long)nEvtALL_g-(long long)nEvtPF_g,
         nEvtPF_g? 100.0*((double)nEvtALL_g-nEvtPF_g)/nEvtPF_g : 0);

  printf("\ndeduped candidate multiplicity (before best-candidate selection):\n");
  printf("  PF only : %.0f\n",nCandPF);
  printf("  PF + LT : %.0f    x%.4f\n",nCandALL, nCandPF? nCandALL/nCandPF : 0);
  printf("\n");
}
