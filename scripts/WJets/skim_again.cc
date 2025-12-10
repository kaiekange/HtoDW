int skim_again(){

    /* ROOT::EnableImplicitMT(); */

    ROOT::RDataFrame myDF("RecoAnalyzer/Events", "/user/kakang/Analysis/CMSSW_10_6_30_patch1/src/tuples/WJets.root");

    myDF.Filter("num_reco_Ds > 0").Snapshot("Events", "/user/kakang/Analysis/CMSSW_10_6_30_patch1/src/tuples/WJets_double_skim.root");

    return 0;
}
