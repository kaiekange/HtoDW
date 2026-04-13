Int_t AddVars()
{
    const TString ChannelSample[14] = {
        "Signal",
        "WJets",
        "W0J",
        "W1J",
        "W2J",
        "TT",
        "TW",
        "WW",
        "WZ",
        "WZto2L2Q",
        "ZZto2L2Q",
        "WWto2L2Nu_OS",
        "WWto2L2Nu_SS_EW",
        "WWto2L2Nu_SS_QCD"};

    // for (Int_t i = 1; i < 14; i++)
    // {
    //     TChain *mychain = new TChain("Events");
    //     mychain->Add("/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples/tuples/" + ChannelSample[i] + "_ReDefine.root");

    //     std::cout << mychain->GetEntries() << ", " << mychain->GetEntries("num_ss_mu == 0 && num_os_mu == 1 && Ds_pt > 15 && Kp_pt > 2 && Km_pt > 2 && pi_pt > 1.5 && H_mass < 125") << std::endl;
    // }

    ROOT::RDataFrame SigDF("FullReco/Events", "/eos/home-k/kakang/H2DW/Analysis/MC/2022/Samples/tuples/Signal.root");
    SigDF.Define("W_masst", "sqrt( 2 * mu_pt * met_pt * ( 1 - cos( mu_phi - met_phi ) ) )")
        .Define("H_px", "Ds_px + mu_px")
        .Define("H_py", "Ds_py + mu_py")
        .Define("H_pz", "Ds_pz + mu_pz")
        .Define("H_pt", "sqrt( H_px * H_px + H_py * H_py )")
        .Define("H_p", "sqrt( H_px * H_px + H_py * H_py + H_pz * H_pz )")
        .Define("H_mass", "sqrt( ( sqrt( Ds_mass * Ds_mass + Ds_p * Ds_p ) + sqrt( 0.1056583755 * 0.1056583755 + mu_p * mu_p) ) * ( sqrt( Ds_mass * Ds_mass + Ds_p * Ds_p ) + sqrt( 0.1056583755 * 0.1056583755 + mu_p * mu_p) ) - H_p * H_p )")
        .Define("K1_pt", "Kp_pt > Km_pt ? Kp_pt : Km_pt")
        .Define("K2_pt", "Kp_pt > Km_pt ? Km_pt : Kp_pt")
        .Define("K1_ip", "Kp_pt > Km_pt ? Kp_ip : Km_ip")
        .Define("K2_ip", "Kp_pt > Km_pt ? Km_ip : Kp_ip")
        .Define("K1_iperr", "Kp_pt > Km_pt ? Kp_iperr : Km_iperr")
        .Define("K2_iperr", "Kp_pt > Km_pt ? Km_iperr : Kp_iperr")
        .Define("K1_ipsig", "Kp_pt > Km_pt ? Kp_ipsig : Km_ipsig")
        .Define("K2_ipsig", "Kp_pt > Km_pt ? Km_ipsig : Kp_ipsig")
        .Define("dR_K1_pi", "Kp_pt > Km_pt ? dR_Kp_pi : dR_Km_pi")
        .Define("dR_K2_pi", "Kp_pt > Km_pt ? dR_Km_pi : dR_Kp_pi")
        .Define("dR_K1_phi", "Kp_pt > Km_pt ? dR_Kp_phi : dR_Km_phi")
        .Define("dR_K2_phi", "Kp_pt > Km_pt ? dR_Km_phi : dR_Kp_phi")
        .Define("dR_K1_Ds", "Kp_pt > Km_pt ? dR_Kp_Ds : dR_Km_Ds")
        .Define("dR_K2_Ds", "Kp_pt > Km_pt ? dR_Km_Ds : dR_Kp_Ds")
        .Define("dR_K1_mu", "Kp_pt > Km_pt ? dR_Kp_mu : dR_Km_mu")
        .Define("dR_K2_mu", "Kp_pt > Km_pt ? dR_Km_mu : dR_Kp_mu")
        .Snapshot("Events", "/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/Signal_ReDefine.root");

    for (Int_t i = 1; i < 14; i++)
    {
        ROOT::RDataFrame BkgDF("FullReco/Events", "/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/" + ChannelSample[i] + ".root");

        BkgDF.Define("W_masst", "sqrt( 2 * mu_pt * met_pt * ( 1 - cos( mu_phi - met_phi ) ) )")
            .Define("H_px", "Ds_px + mu_px")
            .Define("H_py", "Ds_py + mu_py")
            .Define("H_pz", "Ds_pz + mu_pz")
            .Define("H_pt", "sqrt( H_px * H_px + H_py * H_py )")
            .Define("H_p", "sqrt( H_px * H_px + H_py * H_py + H_pz * H_pz )")
            .Define("H_mass", "sqrt( ( sqrt( Ds_mass * Ds_mass + Ds_p * Ds_p ) + sqrt( 0.1056583755 * 0.1056583755 + mu_p * mu_p) ) * ( sqrt( Ds_mass * Ds_mass + Ds_p * Ds_p ) + sqrt( 0.1056583755 * 0.1056583755 + mu_p * mu_p) ) - H_p * H_p )")
            .Define("K1_pt", "Kp_pt > Km_pt ? Kp_pt : Km_pt")
            .Define("K2_pt", "Kp_pt > Km_pt ? Km_pt : Kp_pt")
            .Define("K1_ip", "Kp_pt > Km_pt ? Kp_ip : Km_ip")
            .Define("K2_ip", "Kp_pt > Km_pt ? Km_ip : Kp_ip")
            .Define("K1_iperr", "Kp_pt > Km_pt ? Kp_iperr : Km_iperr")
            .Define("K2_iperr", "Kp_pt > Km_pt ? Km_iperr : Kp_iperr")
            .Define("K1_ipsig", "Kp_pt > Km_pt ? Kp_ipsig : Km_ipsig")
            .Define("K2_ipsig", "Kp_pt > Km_pt ? Km_ipsig : Kp_ipsig")
            .Define("dR_K1_pi", "Kp_pt > Km_pt ? dR_Kp_pi : dR_Km_pi")
            .Define("dR_K2_pi", "Kp_pt > Km_pt ? dR_Km_pi : dR_Kp_pi")
            .Define("dR_K1_phi", "Kp_pt > Km_pt ? dR_Kp_phi : dR_Km_phi")
            .Define("dR_K2_phi", "Kp_pt > Km_pt ? dR_Km_phi : dR_Kp_phi")
            .Define("dR_K1_Ds", "Kp_pt > Km_pt ? dR_Kp_Ds : dR_Km_Ds")
            .Define("dR_K2_Ds", "Kp_pt > Km_pt ? dR_Km_Ds : dR_Kp_Ds")
            .Define("dR_K1_mu", "Kp_pt > Km_pt ? dR_Kp_mu : dR_Km_mu")
            .Define("dR_K2_mu", "Kp_pt > Km_pt ? dR_Km_mu : dR_Kp_mu")
            .Snapshot("Events", "/eos/home-k/kakang/H2DW/Analysis/MC/2022/old_samples/tuples/" + ChannelSample[i] + "_ReDefine.root");
    }
    return 0;
}