Int_t Prepare()
{
    ROOT::EnableImplicitMT();

    ROOT::RDataFrame mydf("DsstWKKpimu/Events", "/eos/home-k/kakang/H2DW/Analysis/MC/2022/SIG_DsstW_KKpimu/SIG_DsstW_KKpimu.root");
    mydf.Filter("event_match").Snapshot("Events", "/eos/home-k/kakang/H2DW/Analysis/MVA/DsstW_KKpimu/Signal.root");
    
    return 0;
}