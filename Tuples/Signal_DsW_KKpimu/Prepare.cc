Int_t Prepare()
{
    ROOT::EnableImplicitMT();

    ROOT::RDataFrame mydf("RecoAllDsWKKpimu/Events", "/eos/home-k/kakang/H2DW/Analysis/MC/2022/SIG_DsW_KKpimu/SIG_DsW_KKpimu.root");
    mydf.Filter("event_match").Snapshot("Events", "/eos/home-k/kakang/H2DW/Analysis/MVA/DsW_KKpimu/signal.root");
    
    return 0;
}