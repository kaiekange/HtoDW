#!/usr/bin/env python3
"""Merged custom-NanoAOD samples: signal (DsW+DsstW), W+jets GenTag (0J/1J/2J; OS and same-sign-pi SS),
ttbar GenTag (2L2Nu + LNu2Q, nominal + ext1; OS and same-sign-pi SS), DY GenTag (0J/1J/2J; OS and same-sign-pi SS), and Muon data OS / SS (Run2022 C-G), each merged over 2022+2022EE into one file.
NanoAOD-pipeline analogue of the old tuplizer-based make_samples.py.

Two modes, driven by condor_submit_make_samples_nano.py:

  worker <list.txt> <out.root> [--correct]
      One condor job = one chunk (<=50 nano files of ONE piece, listed in list.txt).
      Per file: read nano + its apply_corrections.py friend (row i <-> row i), assert
      the join on run/luminosityBlock/event, select, weight, tag; the whole chunk is
      written ONCE (local temp, then copy), so the merged file ends up with one basket
      per branch per chunk instead of one per nano file. Then <out>.json is written
      last: its existence is the success marker the merge checks.
      --correct: first run the Corrections repo's apply_corrections.py on every file of
      the chunk (subprocess, unchanged script, writes the usual friend), then skim. One
      submission instead of two; the friends are still written because MC systematics
      need unselected rows and the skim gets re-run far more often than the corrections.

  merge <group>   (Signal_merged | WtoLNu2Jets_merged | WtoLNu2Jets_SS_merged | TT_merged | TT_SS_merged | DY_merged | DY_SS_merged | Data_OS_merged | Data_SS_merged)
      Reads the chunk lists the submit script wrote (never re-globs), requires every
      chunk's .root+.json, checks that the chunks cover exactly the files on disk and
      share one branch list, hadds them in list order, writes <group>_meta.json (per
      piece sums) and asserts merged entries == sum of selected rows.
      Run with plain python3 on lxplus (no cmsenv: system hadd + user-site uproot), e.g.
      python3 make_samples_nano.py merge DY_merged --region=zmumu
      --year=2023 (default 2022): the 2023 groups (Data_OS_merged, Data_SS_merged), written to
      Samples/2023/ (Samples/<year>/ since 2026-11-12; workers derive the year from the era).

SELECTION (friend columns), plus the phi window:
    pass_mu_nom & pass_muveto & ~jme_jvm_veto & (jme_n_bjet == 0) & (|phi_mass - M_PHI| < PHI_WIN)
2026-11-04: jet veto map and b-veto on the JME jets (jme.full_columns; before: nano jets,
jvm_veto_nom / n_bjet_nom). MC rows that fail the nominal b-veto but pass it in at least one
of the 24 JME jet variations (jme_n_bjet_<v> == 0) go to a second tree EventsSyst (same
branches; `weight` is NaN there), so a variation can move events in as well as out. The jet
veto map is nominal only (no per-variation column). Data has no EventsSyst.
The phi window (|M(KK) - m_phi| < 20 MeV on the post-fit phi_mass branch, on top of the
producer's M(KK) < m_phi + 50 MeV) was added 2026-10-16 (user: keeps 0.949 of signal and
~0.5 of every non-resonant background). The b-tag efficiency maps are measured inside the
same window. Its value is written into every chunk .json and the merged <group>_meta.json.

WEIGHT (background, i.e. W+jets, ttbar and DY):
    weight = genWeight * norm * w_pu_nom * w_mu_nom * jme_w_btag_nom
weight_nobtag (MC) = the same without the b-tag factor; a JME variation's weight is
weight_nobtag * jme_w_btag_<v> (Events and EventsSyst alike).
`norm` is the piece's normalization.json value (xsec*L / (sum genWeight * acceptance
correction), Corrections repo normalization.py).

WEIGHT (signal) -- changed 2026-09-17, genWeight dropped:
    weight = norm * w_pu_nom * w_mu_nom * jme_w_btag_nom
Found 7/869,385 post-selection signal rows (49/1.6M pre-selection) with genWeight =
-23.1645 instead of +23.1645 -- a technical two-valued weight with no cross-section
meaning for this LO-ish EvtGen+Pythia generation (unlike the real amcatnloFXFX NLO
subtraction term W+jets has). It crashed DNN.py's WeightedRandomSampler (torch.multinomial
rejects any negative probability) the first time stage-1 was retrained on this sample.
`norm` for signal now divides by genEventCount (plain N_gen), not sum(genWeight), to stay
consistent with genWeight being absent from the numerator too -- see normalization.py.

SIGNAL IS ON THE SAME ABSOLUTE SCALE AS BACKGROUND here -- unlike the old Samples/
(weight_mode="signal_era", relative, anchored at 1.0; see eras.py SIGNAL_ERA_WEIGHT)
summing `weight` gives an expected yield.

DATA: no genWeight/normalization/SFs, so `weight` = 1.0 (present only so MC and data read
the same way), no `norm`. OS data is BLINDED here: rows with |Ds_mass - M_DS| <
SR_HALF_WIDTH (eras.py, 1.96835 +- 0.04) are dropped after the selection, the same
definition as highest_pt_dedup.apply_row_cuts(blind=True). SS data is not blinded.

PROVENANCE BRANCHES: era (eras.py wire format, 1=2022, 2=2022EE, 3=2023, 4=2023BPix, data too); signal file: channel (0=DsW, 1=DsstW);
W+jets files (OS and SS): jet_bin (0/1/2);
ttbar files (OS and SS): tt_decay (0=TTto2L2Nu, 1=TTtoLNu2Q; nominal and ext1 carry the same tag); DY files (OS and SS): jet_bin (0/1/2). Data has no extra tag -- `run` identifies the run period.

CARRIED: every custom analyzer branch (everything not in a stock NanoAOD collection this
analysis never reads, see _STD_PREFIXES) + run/luminosityBlock/event/genWeight/nTrueInt;
from the friend: w_pu_nom/w_mu_nom/jme_w_btag_nom, `norm`, and the corrected nominal
kinematics renamed *_corr (FRIEND_CORR_KIN; jet/MET ones are JME). The raw nano branches are
kept too. MC: the JME variation columns jme_{num_jets_pt30,jet_leadpt,n_bjet,w_btag}_<v> and
jme_met_{pt,phi}_<v>_prerecoil (uncl: met only).
NOT CARRIED: the other up/down systematic columns (PU, muon, b-tag SF); they stay in the friends.

RECOIL (2026-11-12): met_{pt,phi}_corr and W_tmass_corr = the full MET chain (Type-1 JES+JER, xy,
hadronic-recoil calibration) for DY / W+jets / signal (OS and SS; reference boson Z / W / H,
make_samples_nano.gen_boson_cols), = the pre-recoil values for tt and data. The pre-recoil ones are
kept as met_{pt,phi}_corr_prerecoil / W_tmass_corr_prerecoil (input of recoil_zmumu.py), and the
corrected samples carry gen_V_{px,py} / gen_Vvis_{px,py}. Maps: CR/zmumu_recoil/recoil_<era>.npz
(sha1 in the chunk/merge meta). The JME variation MET columns are NOT recoil-corrected yet.
2026-11-12: the nano-jet *_nano columns are gone (user).

ttcr: refused since 2026-11-04 (not migrated to the JME columns).

Runs inside cmssw-el9 + CMSSW_15_0_5 (uproot, awkward, hadd).
"""
import glob
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eras import M_DS, SR_HALF_WIDTH, TAG_TO_ERA  # noqa: E402

MC_BASE = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/MC"
DATA_BASE = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/Data"
# 2026-11-12 (user): one directory per data-taking year, Samples/<year>/, group and file names
# unchanged inside (2022 = 2022 + 2022EE, 2023 = 2023 + 2023BPix). Default year 2022.
SAMPLES_BASE = "/eos/home-k/kakang/H2DW/Analysis/NanoMVA/Samples"
ERA_YEAR = {"2022": "2022", "2022EE": "2022", "2023": "2023", "2023BPix": "2023"}
YEARS = ("2022", "2023")
OUT_DIR = f"{SAMPLES_BASE}/2022"   # the default year's directory
LIST_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "condor_logs", "make_samples_nano", "lists")
FILES_PER_JOB = 50
M_PHI, PHI_WIN = 1.019460, 0.020  # GeV, phi window of the selection (see docstring)
APPLY_CORRECTIONS = "/afs/cern.ch/work/k/kakang/H2DW/Corrections/corrections/apply_corrections.py"

# Control regions (2026-10-19, plan cr-z-mumu-w-jets-ttbar). A CR writes under
# OUT_DIR/CR/<region>/ and reuses the
# SR chunk lists (same files, same chunking). CR rows are never blinded and carry extra-lepton
# columns (EXTRA_LEP_DOC). Selections, all on top of pass_mu_nom & ~jme_jvm_veto
# (2026-11-04: JME jets, b-veto union -> EventsSyst as in the SR):
#   zmumu: ~pass_muveto & jme_n_bjet == 0 & exactly one extra muon, opposite charge to the
#          best muon and passing tight PF iso like it (xmu_pfIsoId >= 4, 2026-11-04; the
#          count itself stays the loose-iso veto definition); NO phi window (applied, if wanted, at plot time). b-veto removes ttbar
#          and keeps it disjoint from the ttbar mumu category; ~pass_muveto keeps it disjoint
#          from the SR. MC weight as in the SR (w_btag valid at n_bjet == 0).
#   ttcr:  n_bjet_nom >= 1 & phi window; categories (1l / mumu / emu) are cut downstream from
#          the lepton columns. MC weight x w_btag_any_nom (friend column valid at n_bjet >= 1;
#          the first pass 2026-10-19 had no b-tag SF, meta btag_sf_applied was False).
#          REFUSED since 2026-11-04 (nano-jet columns; not migrated to JME).
# Extra muons: corrected pt (MuonScaRe, as prompt_muon_veto) + ID x loose-iso SF (w_xmu,
# stored for ttcr) and ID x tight-iso SF (w_xmu_tight, multiplied into the zmumu weight).
# Extra electrons: raw, no SF yet.
REGIONS = ("sr", "zmumu", "ttcr")
HF_MATCH_DR = 0.3   # ttcr charm-origin match, reco Ds direction vs GenPart c hadron (see hf_origin_cols)
EXTRA_LEP_DOC = """n_xmu, xmu_{pt,eta,phi,charge}: muons other than the best one with tightId,
|eta|<2.4, pfIsoId>=2 (loose PF iso), CORRECTED pt>=15 -- exactly the pass_muveto definition
(n_xmu>=1 == ~pass_muveto, asserted); leading one (by corrected pt) stored, xmu_pt is the
corrected pt, xmu_pt_raw the nano one, xmu_pfIsoId its Muon_pfIsoId; w_xmu = its ID x
loose-iso SF, w_xmu_tight = its ID x tight-iso SF (MC, else 1; meaningful only if
xmu_pfIsoId >= 4). n_xele, xele_*: Electron_cutBased>=4 (tight, iso included), pt>15,
|eta|<2.5. m_mumu / m_emu: best muon (corrected pt) + leading extra muon / electron.
-777 where absent."""


def region_dir(region, year="2022"):
    base = f"{SAMPLES_BASE}/{year}"
    return base if region == "sr" else f"{base}/CR/{region}"

SIGNAL_CHANNEL = {"Signal_DsW_H_phipimu": 0, "Signal_DsstW_H_phipimu": 1}
WJETS_BIN = {
    "WtoLNu2Jets_0J_H_phipimu_GenTag": 0,
    "WtoLNu2Jets_1J_H_phipimu_GenTag": 1,
    "WtoLNu2Jets_2J_H_phipimu_GenTag": 2,
}
# same-sign-pi W+jets (control region for the OS-minus-SS comparison), same jet_bin tag
WJETS_SS_BIN = {f"{s}_SS": b for s, b in WJETS_BIN.items()}
# ttbar decay-filtered datasets -> tt_decay tag (normalized separately, each with its own xsec)
TT_DECAY = {
    "TTto2L2Nu_H_phipimu_GenTag": 0,
    "TTtoLNu2Q_H_phipimu_GenTag": 1,
    # _ext1: same process, more events; normalization.py gives it the nominal's shared norm
    "TTto2L2Nu_ext1_H_phipimu_GenTag": 0,
    "TTtoLNu2Q_ext1_H_phipimu_GenTag": 1,
}
# same-sign-pi ttbar, same tt_decay tag
TT_SS_DECAY = {f"{s}_SS": d for s, d in TT_DECAY.items()}
DY_BIN = {
    "DYto2L2Jets_0J_H_phipimu_GenTag": 0,
    "DYto2L2Jets_1J_H_phipimu_GenTag": 1,
    "DYto2L2Jets_2J_H_phipimu_GenTag": 2,
}
# same-sign-pi DY, same jet_bin tag
DY_SS_BIN = {f"{s}_SS": b for s, b in DY_BIN.items()}
# era -> OS data CRAB tasks (the SS task appends _SS), in merge order. 2023: two primary datasets
# (Muon0/Muon1) x the 22Sep2023_v1..v4 (C) / v1..v2 (D) datasets, task names with _vN.
DATA_TASKS = {
    "2022": [f"Muon_Run2022{r}_NanoH_phipimu" for r in "CD"],
    "2022EE": [f"Muon_Run2022{r}_NanoH_phipimu" for r in "EFG"],
    "2023": [f"Muon{pd}_Run2023C_v{v}_NanoH_phipimu" for pd in (0, 1) for v in (1, 2, 3, 4)],
    "2023BPix": [f"Muon{pd}_Run2023D_v{v}_NanoH_phipimu" for pd in (0, 1) for v in (1, 2)],
}


def data_task(task, same_sign):
    return task + ("_SS" if same_sign else "")


def data_pd(task):
    """Primary-dataset directory of a data task (Data/<era>/<pd>/<task>/...): Muon, Muon0, Muon1."""
    return task.split("_Run")[0]


# data task -> blind? (OS only)
DATA_TASK_BLIND = {data_task(t, ss): not ss
                   for ts in DATA_TASKS.values() for t in ts for ss in (False, True)}

# group -> ordered pieces (era tag, sample); this order is the row order of the merged file
GROUPS = {   # 2022 (the default year)
    "Signal_merged": [(e, s) for e in ("2022", "2022EE") for s in SIGNAL_CHANNEL],
    "WtoLNu2Jets_merged": [(e, s) for e in ("2022", "2022EE") for s in WJETS_BIN],
    "WtoLNu2Jets_SS_merged": [(e, s) for e in ("2022", "2022EE") for s in WJETS_SS_BIN],
    "TT_merged": [(e, s) for e in ("2022", "2022EE") for s in TT_DECAY],
    "TT_SS_merged": [(e, s) for e in ("2022", "2022EE") for s in TT_SS_DECAY],
    "DY_merged": [(e, s) for e in ("2022", "2022EE") for s in DY_BIN],
    "DY_SS_merged": [(e, s) for e in ("2022", "2022EE") for s in DY_SS_BIN],
    "Data_OS_merged": [(e, data_task(t, False)) for e in ("2022", "2022EE") for t in DATA_TASKS[e]],
    "Data_SS_merged": [(e, data_task(t, True)) for e in ("2022", "2022EE") for t in DATA_TASKS[e]],
}
GROUPS_BY_YEAR = {
    "2022": GROUPS,
    "2023": {   # data only so far (no 2023 MC)
        "Data_OS_merged": [(e, data_task(t, False)) for e in ("2023", "2023BPix") for t in DATA_TASKS[e]],
        "Data_SS_merged": [(e, data_task(t, True)) for e in ("2023", "2023BPix") for t in DATA_TASKS[e]],
    },
}


_STD_PREFIXES = (
    "HLT_", "L1_", "L1Reco_step", "L1simulation_step", "Muon_", "Electron_", "Tau_",
    "FatJet_", "Photon_", "Jet_", "LowPtElectron_", "Flag_", "boostedTau_", "SV_",
    "GenJet", "GenPart", "GenVisTau", "GenDressedLepton", "GenIsolatedPhoton",
    "GenProton", "GenVtx", "SubGenJet", "SubJet_", "SoftActivityJet", "CorrT1METJet",
    "IsoTrack", "OtherPV", "PSWeight", "LHE", "PV_", "Pileup_", "RawPuppiMET", "RawMET",
    "PuppiMET", "MET_", "ChsMET", "DeepMETResolutionTune", "DeepMETResponseTune",
    "TrigObj", "btagWeight", "fixedGridRho", "Rho_", "TkMET", "CaloMET", "FsrPhoton_",
    "HTXS_", "Generator_", "HLTriggerFinalPath", "HLTriggerFirstPath", "BeamSpot_",
    "GenMET_", "bunchCrossing", "genTtbarId", "cand_id", "nCorrT1METJet", "nFatJet",
    "nFsrPhoton", "nGenDressedLepton", "nGenIsolatedPhoton", "nGenJet", "nGenJetAK8",
    "nGenPart", "nGenProton", "nGenVisTau", "nIsoTrack", "nLHEPart", "nLHEPdfWeight",
    "nLHEReweightingWeight", "nLHEScaleWeight", "nLowPtElectron", "nOtherPV", "nPSWeight",
    "nPUInt", "nPhoton", "nSV", "nSoftActivityJet", "nSubGenJetAK8", "nSubJet", "nTau",
    "nTrigObj", "nboostedTau",
    "PPSLocalTrack_", "Proton_", "nPPSLocalTrack", "nProton_",  # data only
)
_ALWAYS = ("run", "luminosityBlock", "event", "genWeight", "nTrueInt")

FRIEND_JOIN_COLS = ["run", "luminosityBlock", "event"]
FRIEND_WEIGHT_COLS = ["w_pu_nom", "w_mu_nom", "jme_w_btag_nom"]
FRIEND_CUT_COLS = ["pass_mu_nom", "pass_muveto", "jme_jvm_veto", "jme_n_bjet"]
# 2026-11-04: jet/MET *_corr are the JME columns (jme.full_columns). 2026-11-12: the nano-jet
# *_nano columns are gone; the friend's MET (jme_met_*, jme_W_tmass) is the input of the recoil
# calibration and is carried as met_{pt,phi}_corr_prerecoil / W_tmass_corr_prerecoil, while
# met_{pt,phi}_corr / W_tmass_corr are recoil-corrected for the RECOIL samples (see worker).
FRIEND_CORR_KIN = {
    "mu_pt_nom": "mu_pt_corr",
    "jme_jet_ht_pt20": "jet_ht_pt20_corr",
    "jme_jet_leadpt": "jet_leadpt_corr",
    "jme_met_pt": "met_pt_corr_prerecoil",
    "jme_met_phi": "met_phi_corr_prerecoil",
    "jme_num_jets_pt30": "num_jets_pt30_corr",
    "jme_W_tmass": "W_tmass_corr_prerecoil",
    "jme_n_bjet": "n_bjet_corr",
    "Dsmu_pt_nom": "Dsmu_pt_corr", "Dsmu_eta_nom": "Dsmu_eta_corr",
    "Dsmu_phi_nom": "Dsmu_phi_corr", "Dsmu_mass_nom": "Dsmu_mass_corr",
    "Dsmu_Ht_nom": "Dsmu_Ht_corr",
    "Dsstmu_pt_nom": "Dsstmu_pt_corr", "Dsstmu_eta_nom": "Dsstmu_eta_corr",
    "Dsstmu_phi_nom": "Dsstmu_phi_corr", "Dsstmu_mass_nom": "Dsstmu_mass_corr",
    "Dsstmu_Ht_nom": "Dsstmu_Ht_corr",
}
# JME systematic variation columns (MC), carried under their friend names jme_<x>_<v>.
# Jet variations have all of JME_VAR_JET + JME_VAR_MET, uncl{Up,Down} only JME_VAR_MET.
JME_VAR_JET = ("num_jets_pt30", "jet_leadpt", "n_bjet", "w_btag")
JME_VAR_MET = ("met_pt", "met_phi")
N_JME_VARIATIONS = 26   # 11 JES sources x 2 + JER x 2 + uncl x 2
# Hadronic-recoil calibration (2026-11-12, recoil_zmumu.py): applied to the nominal MET of the
# samples with a reference boson (gen_boson_kind != 0: DY, W+jets, signal; OS and SS); not to tt
# or data. The JME variation MET columns are written WITHOUT recoil, renamed
# jme_met_{pt,phi}_<v>_prerecoil: per-variation maps come later (JME POG: repeat the derivation
# under the JEC/unclustered variations), so they are NOT comparable with the nominal met_*_corr.
RECOIL_DIR = {y: f"{SAMPLES_BASE}/{y}/CR/zmumu_recoil" for y in YEARS}   # maps per year
GEN_V_COLS = ("gen_V_px", "gen_V_py", "gen_Vvis_px", "gen_Vvis_py")

_IDX_RE = re.compile(r"(?:output|nano)_(\d+)\.root$")


def nano_branch_list(keys):
    return sorted({k for k in keys
                   if k in _ALWAYS or not any(k.startswith(p) for p in _STD_PREFIXES)})


def piece_nano_files(era, sample):
    if sample in DATA_TASK_BLIND:
        # CRAB layout <task>/<timestamp>/<0000|0001|...>/nano_N.root; the friend dir is flat
        files = glob.glob(f"{DATA_BASE}/{era}/{data_pd(sample)}/{sample}/*/*/nano_*.root")
        names = [os.path.basename(f) for f in files]
        if len(set(names)) != len(names):
            sys.exit(f"{era}/{sample}: duplicate nano_N basenames across CRAB subdirectories")
    else:
        files = glob.glob(f"{MC_BASE}/{era}/{sample}/output_*.root")
    return sorted(files, key=lambda p: int(_IDX_RE.search(p).group(1)))


def chunk_name(era, sample, k):
    return f"{era}_{sample}_chunk{k}"


def piece_chunks(era, sample):
    """Chunk list files of one piece as written by the submit script, in chunk order."""
    lists = glob.glob(f"{LIST_DIR}/{era}_{sample}_chunk*.txt")
    return sorted(lists, key=lambda p: int(re.search(r"_chunk(\d+)\.txt$", p).group(1)))


def skim_path(era, sample, k, region="sr"):
    return f"{region_dir(region, ERA_YEAR[era])}/skim/{era}/{sample}/chunk_{k}.root"


def parse_nano_path(path):
    """-> (era, sample, friend path, is_data); same friend layout as apply_corrections.py."""
    parts = os.path.abspath(path).split(os.sep)
    if "Data" in parts:
        i = len(parts) - 1 - parts[::-1].index("Data")
        era, task = parts[i + 1], parts[i + 3]
        if (era not in TAG_TO_ERA or task not in DATA_TASK_BLIND or parts[i + 2] != data_pd(task)
                or len(parts) != i + 7 or not _IDX_RE.search(parts[-1])):
            sys.exit(f"not a nano data path: {path}")
        return era, task, os.sep.join(parts[:i + 2] + ["Corrections", task, parts[-1]]), True
    era, sample = parts[-3], parts[-2]
    if parts[-4] != "MC" or era not in TAG_TO_ERA or not _IDX_RE.search(parts[-1]):
        sys.exit(f"not a nano MC path: {path}")
    if not any(sample in t for t in (SIGNAL_CHANNEL, WJETS_BIN, WJETS_SS_BIN, TT_DECAY, TT_SS_DECAY, DY_BIN, DY_SS_BIN)):
        sys.exit(f"sample {sample!r} belongs to no group")
    friend = os.sep.join(parts[:-2] + ["Corrections", sample, parts[-1]])
    return era, sample, friend, False


# ---------------------------------------------------------------- worker
def extra_lepton_cols(nano_path, mu_pt_corr, era, is_data):
    """EXTRA_LEP_DOC columns for every row of one nano file, in nano row order. The best
    muon is re-identified with muo.locate_selected_muon (same code the friend
    used; it exits if mu_pt/eta/charge are not reproduced). Extra muons get the corrected
    pt of muo.prompt_muon_veto -- the same MuonScaRe calls in the same order
    (scale; + resolution smearing on MC) -- so n_xmu >= 1 must equal ~pass_muveto row by
    row; the caller asserts it."""
    import awkward as ak
    import numpy as np
    import uproot
    sys.path.insert(0, os.path.dirname(APPLY_CORRECTIONS))
    import apply_corrections as acorr
    import muo

    with uproot.open(nano_path) as fn:
        a = fn["Events"].arrays(
            ["event", "luminosityBlock", "mu_pt", "mu_eta", "mu_phi", "mu_charge", "Muon_pt",
             "Muon_eta", "Muon_phi", "Muon_charge", "Muon_tightId", "Muon_pfIsoId",
             "Muon_nTrackerLayers", "Electron_pt", "Electron_eta", "Electron_phi",
             "Electron_charge", "Electron_cutBased"])
    # Same correction files and MuonScaRe.cc as apply_corrections.main() (both via muo.py).
    muo.load_scalesmearing(acorr.CAMPAIGN[era])
    cset_z = muo.sf_set(acorr.CAMPAIGN[era])
    idx = muo.locate_selected_muon(a)
    null = acorr.NULL
    n_ev = len(a)

    # ---- extra muons: candidates as in prompt_muon_veto, corrected pt, >= 15 ----
    cand = ((a["Muon_tightId"] == 1) & (abs(a["Muon_eta"]) < muo.MU_ETA_MAX)
            & (a["Muon_pfIsoId"] >= 2) & (ak.local_index(a["Muon_pt"]) != idx))
    evt, lumi = ak.to_numpy(a["event"]), ak.to_numpy(a["luminosityBlock"])
    ev_of = np.repeat(np.arange(n_ev), ak.to_numpy(ak.num(a["Muon_pt"], axis=1)))
    flat = {f: ak.to_numpy(ak.flatten(a[f"Muon_{f}"])) for f in
            ("pt", "eta", "phi", "charge", "nTrackerLayers", "pfIsoId")}
    cflat = ak.to_numpy(ak.flatten(cand))
    pt_corr = np.full(len(cflat), -1.0)
    for j in np.where(cflat)[0]:
        i = ev_of[j]
        pt_corr[j] = muo.pt_corrected(is_data, flat["pt"][j], flat["eta"][j], flat["phi"][j],
                                      flat["charge"][j], flat["nTrackerLayers"][j], evt[i], lumi[i])
    xsel = cflat & (pt_corr >= muo.VETO_MU_PT_MIN)
    n_xmu = np.bincount(ev_of[xsel], minlength=n_ev).astype(np.int32)
    lead = np.full(n_ev, -1)
    order = np.argsort(np.where(xsel, -pt_corr, np.inf), kind="stable")   # highest corrected pt first
    for j in order[:xsel.sum()][::-1]:
        lead[ev_of[j]] = j
    has = lead >= 0
    out = {"n_xmu": n_xmu}
    out["xmu_pt"] = np.where(has, pt_corr[lead], null)
    out["xmu_pt_raw"] = np.where(has, flat["pt"][lead], null)
    out["xmu_eta"] = np.where(has, flat["eta"][lead], null)
    out["xmu_phi"] = np.where(has, flat["phi"][lead], null)
    out["xmu_charge"] = np.where(has, flat["charge"][lead], int(null)).astype(np.int32)
    out["xmu_pfIsoId"] = np.where(has, flat["pfIsoId"][lead], int(null)).astype(np.int32)
    # ID x loose-iso / tight-iso SF of the leading extra muon (MC; 1 elsewhere). Domain
    # |eta|<2.4, pt>=15 for both iso keys: exactly the extra-muon selection, no clamping needed.
    w, w_tight = np.ones(n_ev), np.ones(n_ev)
    if not is_data and has.any():
        e_, p_ = out["xmu_eta"][has], out["xmu_pt"][has]
        w_id = cset_z[muo.MUO_ID_KEY].evaluate(e_, p_, "nominal")
        w[has] = w_id * cset_z[muo.MUO_LOOSE_ISO_KEY].evaluate(e_, p_, "nominal")
        w_tight[has] = w_id * cset_z[muo.MUO_ISO_KEY].evaluate(e_, p_, "nominal")
    out["w_xmu"] = w
    out["w_xmu_tight"] = w_tight

    # ---- extra electrons: raw (no electron corrections/SF yet) ----
    xele = ((a["Electron_cutBased"] >= 4) & (a["Electron_pt"] > 15.0)
            & (abs(a["Electron_eta"]) < 2.5))
    n = ak.to_numpy(ak.sum(xele, axis=1)).astype(np.int32)
    im = ak.argmax(ak.where(xele, a["Electron_pt"], -1.0), axis=1, keepdims=True)
    out["n_xele"] = n
    for f in ("pt", "eta", "phi", "charge"):
        v = ak.to_numpy(ak.flatten(ak.fill_none(a[f"Electron_{f}"][im], null))).astype(
            np.int32 if f == "charge" else np.float64)
        out[f"xele_{f}"] = np.where(n > 0, v, int(null) if f == "charge" else null)

    mu_eta, mu_phi = ak.to_numpy(a["mu_eta"]), ak.to_numpy(a["mu_phi"])
    for tag, mname in (("xmu", "m_mumu"), ("xele", "m_emu")):
        h = out[f"n_{tag}"] > 0   # absent rows hold -777: keep them out of cosh (overflow)
        dphi = np.abs(np.where(h, out[f"{tag}_phi"], mu_phi) - mu_phi)
        dphi = np.where(dphi > np.pi, 2 * np.pi - dphi, dphi)
        deta = np.where(h, out[f"{tag}_eta"], mu_eta) - mu_eta
        m2 = 2 * mu_pt_corr * np.where(h, out[f"{tag}_pt"], 0.0) * (np.cosh(deta) - np.cos(dphi))
        out[mname] = np.where(h, np.sqrt(np.maximum(m2, 0)), null)
    return out


def hf_origin_cols(nano_path, is_data):
    """cand_hf_from_b per nano row: does the candidate's true charm hadron come from a b hadron?
    GenTag gives only the hadron's species (cand_gen_common_hf_pdgId, set when all three
    tracks share one c/b-hadron ancestor); the nano GenPart keeps every c/b hadron with its
    mother index (not their hadronic daughters). Match: the GenPart with the same |pdgId|
    nearest in dR to the reco Ds direction, dR < HF_MATCH_DR; then walk genPartIdxMother up
    for a b hadron (|pdgId| 500-599 or 5000-5999).
    -1 = no truth (data, signal without GenTag, no common CHARM ancestor, or no match); 0 = prompt
    charm; 1 = from b. cand_hf_dR = the matched dR (-777 if none)."""
    import awkward as ak
    import numpy as np
    import uproot

    with uproot.open(nano_path) as fn:
        t = fn["Events"]
        n = t.num_entries
        if is_data or "cand_gen_common_hf_pdgId" not in t.keys():
            return {"cand_hf_from_b": np.full(n, -1, np.int32), "cand_hf_dR": np.full(n, -777.0)}
        a = t.arrays(["cand_gen_common_hf_pdgId", "Ds_eta", "Ds_phi", "GenPart_pdgId",
                      "GenPart_genPartIdxMother", "GenPart_eta", "GenPart_phi"])
    hf = np.abs(ak.to_numpy(a["cand_gen_common_hf_pdgId"]))
    deta_, dphi_ = ak.to_numpy(a["Ds_eta"]), ak.to_numpy(a["Ds_phi"])
    origin = np.full(n, -1, np.int32)
    dr_out = np.full(n, -777.0)
    # charm ancestors only: a common b-hadron ancestor (no charm in between) has no
    # "prompt vs from-b" question, it stays -1
    for i in np.where(((hf >= 400) & (hf < 500)) | ((hf >= 4000) & (hf < 5000)))[0]:
        pid = np.abs(ak.to_numpy(a["GenPart_pdgId"][i]))
        cand = np.where(pid == hf[i])[0]
        if len(cand) == 0:
            continue
        eta = ak.to_numpy(a["GenPart_eta"][i])[cand]
        phi = ak.to_numpy(a["GenPart_phi"][i])[cand]
        dphi = np.abs((phi - dphi_[i] + np.pi) % (2 * np.pi) - np.pi)
        dr = np.hypot(eta - deta_[i], dphi)
        k = int(np.argmin(dr))
        if dr[k] >= HF_MATCH_DR:
            continue
        dr_out[i] = dr[k]
        mo = ak.to_numpy(a["GenPart_genPartIdxMother"][i])
        j, from_b = mo[cand[k]], 0
        while j >= 0:
            if 500 <= pid[j] < 600 or 5000 <= pid[j] < 6000:
                from_b = 1
                break
            j = mo[j]
        origin[i] = from_b
    return {"cand_hf_from_b": origin, "cand_hf_dR": dr_out}


# Gen boson for the hadronic-recoil calibration (2026-11-12, recoil_zmumu.py). Reference boson per
# sample (gen_boson_kind): DY = Z, W+jets = W, signal = H (user decision: the Ds is an H decay
# product, not recoil); other samples (tt) are not recoil-corrected.
#   V     DY/W+jets: sum of the hard-process leptons = status-1 e/mu/nu with fromHardProcess + the
#         last-copy tau with fromHardProcess (no Z/W entry needed: missing in ~1-4% of the events);
#         signal: the last-copy H.
#   V_vis V minus the hard-process neutrinos (status 1, |pdgId| 12/14/16, fromHardProcess or
#         isHardProcessTauDecayProduct).
# Recoil u = -(MET + V_vis), projected on the V direction.
GEN_BOSON_BRANCHES = ["GenPart_pdgId", "GenPart_status", "GenPart_statusFlags", "GenPart_pt", "GenPart_phi"]


def gen_boson_kind(sample):
    """23 / 24 / 25 for the samples the recoil calibration corrects (DY / W+jets / signal), else 0."""
    for prefix, kind in (("DYto2L", 23), ("WtoLNu", 24), ("Signal_", 25)):
        if sample.startswith(prefix):
            return kind
    return 0


def gen_boson_cols(a, kind):
    """{gen_V_px, gen_V_py, gen_Vvis_px, gen_Vvis_py} (float64) from the GenPart branches of `a`.
    Fails loud if an event does not have the expected hard-process content
    (Z: 2 charged leptons/taus, no direct nu; W: 1 + 1 direct nu; H: exactly one last-copy H)."""
    import awkward as ak
    import numpy as np

    pid, st, fl = abs(a["GenPart_pdgId"]), a["GenPart_status"], a["GenPart_statusFlags"]
    hp, hp_tau, last = ((fl >> 8) & 1) == 1, ((fl >> 9) & 1) == 1, ((fl >> 13) & 1) == 1
    is_nu = (st == 1) & ((pid == 12) | (pid == 14) | (pid == 16))
    px = a["GenPart_pt"] * np.cos(a["GenPart_phi"])
    py = a["GenPart_pt"] * np.sin(a["GenPart_phi"])

    def vsum(m):
        return (ak.to_numpy(ak.sum(ak.where(m, px, 0.0), axis=1)).astype(np.float64),
                ak.to_numpy(ak.sum(ak.where(m, py, 0.0), axis=1)).astype(np.float64))

    def count(m):
        return ak.to_numpy(ak.sum(m, axis=1))

    if kind == 25:
        v = (pid == 25) & last
        bad = count(v) != 1
    elif kind in (23, 24):
        charged = ((st == 1) & ((pid == 11) | (pid == 13)) & hp) | ((pid == 15) & last & hp)
        nu_direct = is_nu & hp
        v = charged | nu_direct
        n_ch, n_nu = count(charged), count(nu_direct)
        bad = ((n_ch != 2) | (n_nu != 0)) if kind == 23 else ((n_ch != 1) | (n_nu != 1))
    else:
        sys.exit(f"gen_boson_cols: no reference boson for kind {kind}")
    if bad.any():
        sys.exit(f"gen_boson_cols (kind {kind}): {int(bad.sum())}/{len(bad)} events without the expected "
                 f"hard-process content, first index {int(np.where(bad)[0][0])}")
    vx, vy = vsum(v)
    nx, ny = vsum(is_nu & (hp | hp_tau))
    return {"gen_V_px": vx, "gen_V_py": vy, "gen_Vvis_px": vx - nx, "gen_Vvis_py": vy - ny}


def recoil_project(qx, qy, vx, vy, met_px, met_py):
    """(q_T, x_par, u_perp): u = -(MET + V_vis) projected on q_hat; x_par = u_par + q_T."""
    import numpy as np
    qt = np.hypot(qx, qy)
    if not np.all(qt > 0):
        sys.exit("recoil_project: q_T = 0")
    hx, hy = qx / qt, qy / qt
    ux, uy = -(met_px + vx), -(met_py + vy)
    return qt, ux * hx + uy * hy + qt, -ux * hy + uy * hx


def recoil_load(path):
    import numpy as np
    with np.load(path) as f:
        return {k: f[k] for k in f.files}


def recoil_correct(P, njet, qx, qy, vx, vy, met_px, met_py):
    """MET (px, py) after the recoil calibration of recoil_zmumu.py (CDF quantile mapping, nominal only).
    P = recoil_load(recoil_<era>.npz): per cell (N_jet class, q_T bin) the CDF of the simulation (src) and of
    data - tt (tgt) on u_edges, for x_par = u_par + q_T and for u_perp. x' = CDF_tgt^-1(CDF_src(x));
    values outside u_edges are left unchanged. MET' = -(u' + V_vis)."""
    import numpy as np
    qt, xpar, uperp = recoil_project(qx, qy, vx, vy, met_px, met_py)
    hx, hy = qx / qt, qy / qt
    cls = np.full(len(qt), -1)
    for k, (lo, hi) in enumerate(zip(P["nj_lo"], P["nj_hi"])):
        cls[(njet >= lo) & (njet <= hi)] = k
    if (cls < 0).any():
        sys.exit(f"recoil_correct: {int((cls < 0).sum())} rows outside the N_jet classes")
    edges, ub = P["qt_edges"], P["u_edges"]
    iq = np.clip(np.digitize(qt, edges) - 1, 0, len(edges) - 2)
    out = {}
    for c, x in (("xpar", xpar), ("uperp", uperp)):
        y = x.copy()
        inside = (x > ub[0]) & (x < ub[-1])
        for k in range(len(P["nj_lo"])):
            for i in range(len(edges) - 1):
                s = (cls == k) & (iq == i) & inside
                cs, ct = P[f"src_{c}"][k, i], P[f"tgt_{c}"][k, i]
                keep = np.concatenate([[True], np.diff(ct) > 0])
                y[s] = np.interp(np.interp(x[s], ub, cs), ct[keep], ub[keep])
        out[c] = y
    upar = out["xpar"] - qt
    ux = upar * hx - out["uperp"] * hy
    uy = upar * hy + out["uperp"] * hx
    return -(ux + vx), -(uy + vy)


def worker(list_path, out_path, correct=False, region="sr"):
    # numpy/uproot only here and in merge(): the submit script imports this module from
    # lxplus's native python, which has neither.
    import numpy as np
    import uproot

    with open(list_path) as f:
        nano_files = [line.strip() for line in f if line.strip()]
    if not nano_files:
        sys.exit(f"empty list {list_path}")
    pieces = {parse_nano_path(p)[:2] for p in nano_files}
    if len(pieces) != 1:
        sys.exit(f"{list_path} mixes pieces: {sorted(pieces)}")
    era, sample, _, is_data = parse_nano_path(nano_files[0])

    if region not in REGIONS:
        sys.exit(f"unknown region {region!r}, choose from {REGIONS}")
    if region == "ttcr":
        sys.exit("ttcr is not migrated to the JME columns (2026-11-04): its n_bjet >= 1 "
                 "selection has no per-variation b-tag weight in the friend")
    tags = {"era": TAG_TO_ERA[era]}
    if is_data:
        norm = xsec_pb = lumi_pb = None
        blind = DATA_TASK_BLIND[sample] and region == "sr"
        friend_cols = FRIEND_JOIN_COLS + FRIEND_CUT_COLS + list(FRIEND_CORR_KIN)
    else:
        with open(f"{MC_BASE}/{era}/Corrections/{sample}/normalization.json") as f:
            nd = json.load(f)
        norm, xsec_pb, lumi_pb = nd["norm"], nd["xsec_pb"], nd["lumi_pb"]
        blind = False
        friend_cols = (FRIEND_JOIN_COLS + FRIEND_WEIGHT_COLS + FRIEND_CUT_COLS
                       + list(FRIEND_CORR_KIN))
        if sample in SIGNAL_CHANNEL:
            tags["channel"] = SIGNAL_CHANNEL[sample]
        elif sample in TT_DECAY:
            tags["tt_decay"] = TT_DECAY[sample]
        elif sample in TT_SS_DECAY:
            tags["tt_decay"] = TT_SS_DECAY[sample]
        elif sample in DY_BIN:
            tags["jet_bin"] = DY_BIN[sample]
        elif sample in DY_SS_BIN:
            tags["jet_bin"] = DY_SS_BIN[sample]
        elif sample in WJETS_SS_BIN:
            tags["jet_bin"] = WJETS_SS_BIN[sample]
        else:
            tags["jet_bin"] = WJETS_BIN[sample]
    kind = 0 if is_data else gen_boson_kind(sample)
    recoil_path = recoil_sha1 = recoil_P = None
    if kind:
        recoil_path = f"{RECOIL_DIR[ERA_YEAR[era]]}/recoil_{era}.npz"
        if not os.path.exists(recoil_path):
            sys.exit(f"missing recoil map {recoil_path}")
        with open(recoil_path, "rb") as f:
            recoil_sha1 = hashlib.sha1(f.read()).hexdigest()
        recoil_P = recoil_load(recoil_path)

    branches = None
    jme_vars = jet_vars = var_cols = None   # MC: fixed by the first friend, required identical after
    parts, nom_masks = [], []
    n_total = n_sel_total = n_syst_total = n_blinded = 0
    wsum = 0.0
    for nano_path in nano_files:
        friend_path = parse_nano_path(nano_path)[2]
        if correct:
            subprocess.run([sys.executable, APPLY_CORRECTIONS, nano_path], check=True)
        with uproot.open(nano_path) as fn:
            tn = fn["Events"]
            if branches is None:
                branches = nano_branch_list(tn.keys())
            nano = tn.arrays(branches, library="np")
            gen = tn.arrays(GEN_BOSON_BRANCHES) if kind else None
        with uproot.open(friend_path) as ff:
            tf = ff["Events"]
            if not is_data:
                keys = set(tf.keys())
                v_here = sorted(m.group(1) for k in keys
                                for m in [re.match(r"jme_met_pt_(.+)$", k)] if m)
                if jme_vars is None:
                    jme_vars = v_here
                    jet_vars = [v for v in jme_vars if f"jme_n_bjet_{v}" in keys]
                    var_cols = [f"jme_{x}_{v}" for v in jme_vars
                                for x in (JME_VAR_JET if v in jet_vars else ()) + JME_VAR_MET]
                    if len(jme_vars) != N_JME_VARIATIONS or len(jet_vars) != N_JME_VARIATIONS - 2:
                        sys.exit(f"{friend_path}: {len(jme_vars)} JME variations "
                                 f"({len(jet_vars)} with jets), expected {N_JME_VARIATIONS} "
                                 f"({N_JME_VARIATIONS - 2})")
                elif v_here != jme_vars:
                    sys.exit(f"{friend_path}: JME variation set differs from the chunk's first friend")
            fr = tf.arrays(friend_cols + (var_cols or []), library="np")
        n = len(nano["event"])
        n_total += n
        if len(fr["event"]) != n:
            sys.exit(f"row count mismatch {nano_path} ({n}) vs {friend_path} ({len(fr['event'])})")
        for c in FRIEND_JOIN_COLS:
            if not np.array_equal(nano[c], fr[c]):
                sys.exit(f"join mismatch on {c}: {nano_path} vs {friend_path}")

        # everything but the b-veto; the b-veto decides Events (nominal) vs EventsSyst
        if region == "sr":
            base = (fr["pass_mu_nom"] & fr["pass_muveto"] & ~fr["jme_jvm_veto"]
                    & (np.abs(nano["phi_mass"] - M_PHI) < PHI_WIN))
        else:  # zmumu
            lep = extra_lepton_cols(nano_path, fr["mu_pt_nom"], era, is_data)
            if not np.array_equal(lep["n_xmu"] >= 1, ~fr["pass_muveto"]):
                sys.exit(f"{nano_path}: extra-muon count disagrees with pass_muveto on "
                         f"{int(((lep['n_xmu'] >= 1) != ~fr['pass_muveto']).sum())} rows")
            base = (fr["pass_mu_nom"] & ~fr["jme_jvm_veto"] & ~fr["pass_muveto"]
                    & (lep["n_xmu"] == 1) & (lep["xmu_charge"] == -nano["mu_charge"])
                    & (lep["xmu_pfIsoId"] >= 4))
        bveto = fr["jme_n_bjet"] == 0
        any_bveto = bveto.copy()
        for v in jet_vars or []:   # MC: keep the union over the JME jet variations
            any_bveto |= fr[f"jme_n_bjet_{v}"] == 0
        sel = base & any_bveto
        if blind:   # data only, where sel == base & bveto
            in_sr = sel & (np.abs(nano["Ds_mass"] - M_DS) < SR_HALF_WIDTH)
            n_blinded += int(in_sr.sum())
            sel &= ~in_sr
        n_keep = int(sel.sum())
        if n_keep == 0:
            continue
        is_nom = bveto[sel]
        if is_data:
            w = np.ones(n_keep, dtype=np.float64)
        else:
            # Signal drops genWeight (see module docstring, WEIGHT (signal)): the
            # ~23.1645-magnitude technical weight isn't a real cross-section factor for
            # this generation, and a rare sign flip (7/869,385 post-selection) crashes
            # torch's WeightedRandomSampler downstream. normalization.py's `norm` for
            # signal is computed to match (genEventCount denominator, not sum(genWeight)).
            gw = (np.ones(n_keep, dtype=np.float64) if sample in SIGNAL_CHANNEL
                  else nano["genWeight"][sel].astype(np.float64))
            w_nob = gw * norm * fr["w_pu_nom"][sel] * fr["w_mu_nom"][sel]
            if region == "zmumu":  # every zmumu row has exactly one extra muon, tight iso
                w_nob = w_nob * lep["w_xmu_tight"][sel]
            if not np.all(np.isfinite(w_nob)):
                sys.exit(f"non-finite weight_nobtag in {nano_path} after selection")
            # jme_w_btag_nom is NaN by design at jme_n_bjet > 0, i.e. on every EventsSyst row
            w = w_nob * fr["jme_w_btag_nom"][sel]
            if not np.all(np.isfinite(w[is_nom])):
                sys.exit(f"non-finite weight in {nano_path} after selection")
        out = {b: nano[b][sel] for b in branches}
        out["weight"] = w
        if not is_data:
            out["weight_nobtag"] = w_nob
            out["norm"] = np.full(n_keep, norm, dtype=np.float64)
            for c in var_cols:   # JME variation MET: no recoil yet (see RECOIL_DIR comment)
                out[f"{c}_prerecoil" if c.startswith("jme_met_") else c] = fr[c][sel]
        for fcol, ocol in FRIEND_CORR_KIN.items():
            out[ocol] = fr[fcol][sel]
        # recoil calibration of the nominal MET; W_tmass rebuilt with the friend's formula (closure)
        mu_pt_c, mu_phi = out["mu_pt_corr"].astype(np.float64), nano["mu_phi"][sel].astype(np.float64)

        def w_tmass(pt, phi):
            return np.sqrt(2.0 * mu_pt_c * pt * (1.0 - np.cos(mu_phi - phi)))
        mpt0 = out["met_pt_corr_prerecoil"].astype(np.float64)
        mphi0 = out["met_phi_corr_prerecoil"].astype(np.float64)
        ref = out["W_tmass_corr_prerecoil"]
        if (np.abs(w_tmass(mpt0, mphi0) - ref) > 1e-6 * np.maximum(ref, 1.0)).any():
            sys.exit(f"{nano_path}: W_tmass formula does not reproduce the friend's jme_W_tmass")
        if kind:
            g = gen_boson_cols(gen[sel], kind)
            px, py = recoil_correct(recoil_P, out["num_jets_pt30_corr"], g["gen_V_px"], g["gen_V_py"],
                                    g["gen_Vvis_px"], g["gen_Vvis_py"], mpt0 * np.cos(mphi0), mpt0 * np.sin(mphi0))
            mpt, mphi = np.hypot(px, py), np.arctan2(py, px)
            for c in GEN_V_COLS:
                out[c] = g[c].astype(np.float32)
        else:
            mpt, mphi = mpt0, mphi0
        out["met_pt_corr"], out["met_phi_corr"], out["W_tmass_corr"] = mpt, mphi, w_tmass(mpt, mphi)
        for c in ("met_pt_corr", "met_phi_corr", "W_tmass_corr"):
            if not np.all(np.isfinite(out[c])):
                sys.exit(f"{nano_path}: non-finite {c} after the recoil step")
        for k, v in tags.items():
            out[k] = np.full(n_keep, v, dtype=np.int32)
        if region != "sr":
            for k, v in lep.items():
                out[k] = v[sel]
        parts.append(out)
        nom_masks.append(is_nom)
        n_sel_total += int(is_nom.sum())
        n_syst_total += n_keep - int(is_nom.sum())
        wsum += float(w[is_nom].sum())

    meta = dict(era=era, sample=sample, list=os.path.basename(list_path),
                files=[os.path.basename(p) for p in nano_files], n_files=len(nano_files),
                n_rows_total=n_total, n_rows_selected=n_sel_total,
                n_rows_selected_syst=n_syst_total, n_rows_blinded=n_blinded,
                weight_sum=wsum,
                norm=norm, xsec_pb=xsec_pb, lumi_pb=lumi_pb, phi_win=PHI_WIN,
                branches_sha1=None, root_written=False,
                # 2026-11-04: jet/MET *_corr, the jet veto map, n_bjet and w_btag are the JME
                # columns of the friend (jme.full_columns); models trained before must not
                # score these samples.
                jme_applied=True, jme_variations=jme_vars,
                # 2026-11-12: recoil calibration of the nominal MET (None = sample not corrected)
                recoil_kind=kind, recoil_map=recoil_path, recoil_map_sha1=recoil_sha1)
    if region != "sr":
        meta.update(region=region, phi_win=None,
                    btag_sf_applied=True, extra_muon_corrected=True,
                    extra_muon_sf_in_weight=True, extra_electron_sf_applied=False,
                    zmumu_xmu_iso="tight")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    if parts:
        cols = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
        nom = np.concatenate(nom_masks)
        trees = {"Events": {k: v[nom] for k, v in cols.items()}}
        if not is_data:   # MC: rows passing the b-veto only in some JME variation
            trees["EventsSyst"] = {k: v[~nom] for k, v in cols.items()}
        meta["branches_sha1"] = hashlib.sha1("\n".join(sorted(cols)).encode()).hexdigest()
        with tempfile.TemporaryDirectory() as tmp:
            local = os.path.join(tmp, "chunk.root")
            with uproot.recreate(local) as f:
                for name, tc in trees.items():
                    f[name] = tc
            shutil.copyfile(local, out_path)
        with uproot.open(out_path) as f:
            for name, n_exp in (("Events", n_sel_total), ("EventsSyst", n_syst_total)):
                if name not in trees:
                    continue
                t = f[name]
                back = t.arrays(["event", "weight"], library="np")
                if (t.num_entries != n_exp
                        or not np.array_equal(back["event"], trees[name]["event"])
                        or not np.array_equal(back["weight"], trees[name]["weight"],
                                              equal_nan=True)):
                    sys.exit(f"read-back of {out_path}:{name} does not match")
        meta["root_written"] = True
    elif os.path.exists(out_path):
        os.remove(out_path)
    with open(out_path[:-len(".root")] + ".json", "w") as f:
        json.dump(meta, f, indent=1)
    print(f"{era}/{sample} {os.path.basename(list_path)}: files={len(nano_files)} "
          f"rows={n_total} sel={n_sel_total} sel_syst={n_syst_total} blinded={n_blinded} "
          f"weight_sum={wsum:.6g} -> {out_path}")


# ---------------------------------------------------------------- merge
def merge(group, region="sr", year="2022"):
    import uproot

    if region not in REGIONS:
        sys.exit(f"unknown region {region!r}, choose from {REGIONS}")

    if year not in GROUPS_BY_YEAR or group not in GROUPS_BY_YEAR[year]:
        sys.exit(f"unknown year/group {year}/{group!r}, choose from "
                 f"{ {y: sorted(g) for y, g in GROUPS_BY_YEAR.items()} }")
    roots, piece_meta, shas, problems = [], [], set(), []
    recoil_shas = {}   # era -> recoil map sha1s used by the chunks (must be one)
    for era, sample in GROUPS_BY_YEAR[year][group]:
        lists = piece_chunks(era, sample)
        if not lists:
            problems.append(f"{era}/{sample}: no chunk lists in {LIST_DIR}")
            continue
        listed = []
        agg = dict(era=era, sample=sample, n_chunks=len(lists), n_files=0, n_rows_total=0,
                   n_rows_selected=0, n_rows_selected_syst=0, n_rows_blinded=0, weight_sum=0.0,
                   norm=None, xsec_pb=None, lumi_pb=None)
        ks = [int(re.search(r"_chunk(\d+)\.txt$", lp).group(1)) for lp in lists]
        if ks != list(range(len(lists))):
            problems.append(f"{era}/{sample}: chunk numbers {ks} are not 0..{len(lists) - 1}")
            continue
        for k, lp in zip(ks, lists):
            with open(lp) as f:
                names = [os.path.basename(x.strip()) for x in f if x.strip()]
            listed += names
            rp = skim_path(era, sample, k, region)
            jp = rp[:-len(".root")] + ".json"
            if not os.path.exists(jp):
                problems.append(f"missing {jp}")
                continue
            with open(jp) as f:
                m = json.load(f)
            if m["files"] != names:
                problems.append(f"{jp}: file list differs from {lp}")
                continue
            if not m.get("jme_applied") or "n_rows_selected_syst" not in m:
                problems.append(f"{jp}: written before the 2026-11-04 JME skim, re-run it")
                continue
            if "recoil_map_sha1" not in m:
                problems.append(f"{jp}: written before the 2026-11-12 recoil skim, re-run it")
                continue
            recoil_shas.setdefault(era, set()).add(m["recoil_map_sha1"])
            if m["root_written"]:
                if not os.path.exists(rp):
                    problems.append(f"missing {rp}")
                    continue
                roots.append(rp)
                shas.add(m["branches_sha1"])
            for key in ("n_files", "n_rows_total", "n_rows_selected", "n_rows_selected_syst",
                        "n_rows_blinded", "weight_sum"):
                agg[key] += m[key]
            for key in ("norm", "xsec_pb", "lumi_pb", "phi_win", "jme_applied", "recoil_map_sha1"):
                agg[key] = m.get(key)
        on_disk = [os.path.basename(p) for p in piece_nano_files(era, sample)]
        if sorted(listed) != sorted(on_disk) or len(set(listed)) != len(listed):
            problems.append(f"{era}/{sample}: chunk lists ({len(listed)}) do not cover the "
                            f"{len(on_disk)} nano files on disk exactly")
        piece_meta.append(agg)
    if len(shas) > 1:
        problems.append(f"chunks disagree on the branch list ({len(shas)} variants)")
    for era, s in recoil_shas.items():
        if len(s) > 1:
            problems.append(f"{era}: chunks used {len(s)} different recoil maps, re-run the stale ones")
    if problems:
        sys.exit("merge refused:\n  " + "\n  ".join(problems))

    target = f"{region_dir(region, year)}/{group}.root"
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
        local = os.path.join(tmp, f"{group}.root")
        subprocess.run(["hadd", "-f", local] + roots, check=True)
        shutil.copyfile(local, target)
    n_expect = sum(p["n_rows_selected"] for p in piece_meta)
    n_expect_syst = sum(p["n_rows_selected_syst"] for p in piece_meta)
    is_data_group = group.startswith("Data_")
    with uproot.open(target) as f:
        n_got = f["Events"].num_entries
        has_syst = "EventsSyst" in f
        n_got_syst = f["EventsSyst"].num_entries if has_syst else 0
    if n_got != n_expect:
        sys.exit(f"{target}: {n_got} entries, chunk JSONs sum to {n_expect}")
    if has_syst == is_data_group or n_got_syst != n_expect_syst:
        sys.exit(f"{target}: EventsSyst present={has_syst} with {n_got_syst} entries, "
                 f"expected {'none' if is_data_group else n_expect_syst}")
    with open(f"{region_dir(region, year)}/{group}_meta.json", "w") as f:
        json.dump(piece_meta, f, indent=1)
    for p in piece_meta:
        xl = "" if p["xsec_pb"] is None else f" xsec*L={p['xsec_pb'] * p['lumi_pb']:.6g}"
        print(f"  {p['era']}/{p['sample']}: chunks={p['n_chunks']} files={p['n_files']} "
              f"rows={p['n_rows_total']} sel={p['n_rows_selected']} "
              f"sel_syst={p['n_rows_selected_syst']} "
              f"blinded={p['n_rows_blinded']} weight_sum={p['weight_sum']:.6g}{xl}")
    print(f"wrote {target} ({n_got} entries, EventsSyst {n_got_syst}, from {len(roots)} chunks)")


def main():
    argv = [a for a in sys.argv[1:] if not a.startswith(("--region=", "--year="))]
    regions = [a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--region=")]
    region = regions[-1] if regions else "sr"
    years = [a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--year=")]
    year = years[-1] if years else "2022"
    if len(argv) in (3, 4) and argv[0] == "worker" and argv[3:] in ([], ["--correct"]):
        worker(argv[1], argv[2], correct=len(argv) == 4, region=region)
    elif len(argv) == 2 and argv[0] == "merge":
        merge(argv[1], region, year)
    else:
        sys.exit("usage: make_samples_nano.py worker <list.txt> <out.root> [--correct] [--region=R]\n"
                 "       make_samples_nano.py merge <" + "|".join(GROUPS) + "> [--region=R] [--year=Y]\n"
                 f"       R in {REGIONS}, default sr; Y in {YEARS}, default 2022 (2023: Data_* only)")


if __name__ == "__main__":
    main()
