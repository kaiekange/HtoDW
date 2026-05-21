Int_t Prepare()
{
    ROOT::EnableImplicitMT();

    ROOT::RDataFrame mydf("DsWKKpimu/Events", "/eos/home-k/kakang/H2DW/Analysis/MC/2022/SIG_DsW_KKpimu/SIG_DsW_KKpimu.root");
    mydf.Filter("event_match").Snapshot("Events", "/eos/home-k/kakang/H2DW/Analysis/MVA/DsW_KKpimu/Signal.root");
    
    return 0;
}