"""Gen-level event counting on MiniAOD, before the JME jet-veto-map filter.

The tuplizer counters (nEventsProcessed_/nEventsLepton_/nEventsGen_ in the Info
tree) are all filled *after* jvmAppliedEventFilter has already dropped events, so
they cannot be used as a gen-level normalisation denominator. This reads the PAT
files directly and replicates the gen filter of
EDAnalyzers/Reconstruction/plugins/H_phipimu.cc, with no veto and no reco.

The sample contains both W->munu and W->enu (roughly half each). The tuplizer only
ever keeps the munu half (its hard-process cut is 1mu + 0e), so both flavours are
counted here. The Ds -> phi pi -> KK pi truth-match is independent of the lepton
flavour -- the decay chain only involves the Ds/W pdgIds -- so it is evaluated once
per event and booked against whichever flavour the event has.

Counts five numbers per shard:
  total   - MiniAOD events read
  mu      - hard process = 1H + 1Ds(*) + 1W + 1mu + 0e   (W->munu)
  mu_gen  - mu,  and additionally truth-matched Ds -> phi pi -> KK pi
  e       - hard process = 1H + 1Ds(*) + 1W + 1e  + 0mu  (W->enu)
  e_gen   - e,   and additionally truth-matched Ds -> phi pi -> KK pi
            (EvtGen does not force every Ds down this chain)

Usage: python3 gencount.py <DsW|DsstW> <JOBID> <GROUPSIZE> <OUTFILE>
File indices follow the tuplizer convention: JOBID*GROUPSIZE+1 ... +GROUPSIZE.
"""
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.FWLiteEnabler.enable()
from DataFormats.FWLite import Events, Handle

CHANNEL = sys.argv[1]
JOBID = int(sys.argv[2])
GROUPSIZE = int(sys.argv[3])
OUTFILE = sys.argv[4]

SAMPLES = {
    "DsW":   ("ggHtoDsW_phipi",   431, 1000, set()),
    "DsstW": ("ggHtoDsstW_phipi", 433, 1000, {188}),  # output_188 does not exist
}
DATASET, RESID, NFILES, SKIP = SAMPLES[CHANNEL]
BASE = ("root://maite.iihe.ac.be//store/user/kakang/Analysis/"
        "{}/20260624/Simulation/PAT/2022".format(DATASET))


def TwoBodyDecay(gp, motherid, daughter1id, daughter2id):
    """Recursive mother-walk; mirrors H_phipimu::TwoBodyDecay."""
    for i in range(gp.numberOfMothers()):
        mother = gp.mother(i)
        if not mother:
            continue
        if mother.pdgId() == motherid and mother.numberOfDaughters() >= 2:
            f1 = f2 = False
            for j in range(mother.numberOfDaughters()):
                pdgid = mother.daughter(j).pdgId()
                if pdgid == daughter1id:
                    f1 = True
                if pdgid == daughter2id:
                    f2 = True
            if f1 and f2:
                return True
        if TwoBodyDecay(mother, motherid, daughter1id, daughter2id):
            return True
    return False


start = JOBID * GROUPSIZE + 1
end = start + GROUPSIZE
infiles = ["{}/output_{}.root".format(BASE, i)
           for i in range(start, end) if i <= NFILES and i not in SKIP]

nEventsTotal = 0
nEventsMu = nEventsMuGen = 0
nEventsE = nEventsEGen = 0

if infiles:
    GenPartHandle = Handle("std::vector<reco::GenParticle>")

    for event in Events(infiles):
        nEventsTotal += 1
        event.getByLabel("prunedGenParticles", GenPartHandle)
        GenPart = GenPartHandle.product()

        Gen_H_idx = Gen_R_idx = Gen_W_idx = 0
        num_Gen_H = num_Gen_R = num_Gen_W = num_Gen_mu = num_Gen_e = 0
        for i, gp in enumerate(GenPart):
            if not gp.isHardProcess():
                continue
            pdgid = gp.pdgId()
            if pdgid == 25:
                Gen_H_idx = i
                num_Gen_H += 1
            elif abs(pdgid) == RESID:
                Gen_R_idx = i
                num_Gen_R += 1
            elif abs(pdgid) == 24:
                Gen_W_idx = i
                num_Gen_W += 1
            elif abs(pdgid) == 13:
                num_Gen_mu += 1
            elif abs(pdgid) == 11:
                num_Gen_e += 1
        if not (num_Gen_H == 1 and num_Gen_R == 1 and num_Gen_W == 1):
            continue
        if num_Gen_mu == 1 and num_Gen_e == 0:
            flavour = "mu"
            nEventsMu += 1
        elif num_Gen_e == 1 and num_Gen_mu == 0:
            flavour = "e"
            nEventsE += 1
        else:
            continue

        R_gp = GenPart[Gen_R_idx]          # Ds (DsW) or Ds* (DsstW)
        W_gp = GenPart[Gen_W_idx]
        R_pdgid = R_gp.pdgId()
        W_pdgid = W_gp.pdgId()
        R_charge = R_gp.charge()

        num_Gen_Ds = num_Gen_gamma = 0
        num_Gen_Kp = num_Gen_Km = num_Gen_pi = 0

        if CHANNEL == "DsW":
            for gp in GenPart:
                if gp.status() != 1:
                    continue
                pdgid = gp.pdgId()
                if pdgid == R_charge * 321:
                    if (TwoBodyDecay(gp, 333, 321, -321)
                            and TwoBodyDecay(gp, R_pdgid, 333, R_charge * 211)
                            and TwoBodyDecay(gp, 25, R_pdgid, W_pdgid)):
                        num_Gen_Kp += 1
                elif pdgid == -R_charge * 321:
                    if (TwoBodyDecay(gp, 333, 321, -321)
                            and TwoBodyDecay(gp, R_pdgid, 333, R_charge * 211)
                            and TwoBodyDecay(gp, 25, R_pdgid, W_pdgid)):
                        num_Gen_Km += 1
                elif pdgid == R_charge * 211:
                    if (TwoBodyDecay(gp, R_pdgid, 333, R_charge * 211)
                            and TwoBodyDecay(gp, 25, R_pdgid, W_pdgid)):
                        num_Gen_pi += 1
            if num_Gen_Kp == 1 and num_Gen_Km == 1 and num_Gen_pi == 1:
                if flavour == "mu":
                    nEventsMuGen += 1
                else:
                    nEventsEGen += 1

        else:  # DsstW
            for gp in GenPart:
                status = gp.status()
                pdgid = gp.pdgId()
                if status == 2 and pdgid == R_charge * 431:
                    if (TwoBodyDecay(gp, R_pdgid, R_charge * 431, 22)
                            and TwoBodyDecay(gp, 25, R_pdgid, W_pdgid)
                            and gp.numberOfDaughters() > 1):
                        d1 = d2 = False
                        for j in range(gp.numberOfDaughters()):
                            dpdgid = gp.daughter(j).pdgId()
                            if dpdgid == 333:
                                d1 = True
                            if dpdgid == R_charge * 211:
                                d2 = True
                        if d1 and d2:
                            num_Gen_Ds += 1
                elif status == 1:
                    if pdgid == R_charge * 321:
                        if (TwoBodyDecay(gp, 333, 321, -321)
                                and TwoBodyDecay(gp, R_charge * 431, 333, R_charge * 211)
                                and TwoBodyDecay(gp, R_pdgid, R_charge * 431, 22)
                                and TwoBodyDecay(gp, 25, R_pdgid, W_pdgid)):
                            num_Gen_Kp += 1
                    elif pdgid == -R_charge * 321:
                        if (TwoBodyDecay(gp, 333, 321, -321)
                                and TwoBodyDecay(gp, R_charge * 431, 333, R_charge * 211)
                                and TwoBodyDecay(gp, R_pdgid, R_charge * 431, 22)
                                and TwoBodyDecay(gp, 25, R_pdgid, W_pdgid)):
                            num_Gen_Km += 1
                    elif pdgid == R_charge * 211:
                        if (TwoBodyDecay(gp, R_charge * 431, 333, R_charge * 211)
                                and TwoBodyDecay(gp, R_pdgid, R_charge * 431, 22)
                                and TwoBodyDecay(gp, 25, R_pdgid, W_pdgid)):
                            num_Gen_pi += 1
                    elif pdgid == 22:
                        if (gp.mother(0).pdgId() == R_pdgid
                                and TwoBodyDecay(gp, R_pdgid, R_charge * 431, 22)
                                and TwoBodyDecay(gp, 25, R_pdgid, W_pdgid)):
                            num_Gen_gamma += 1
            if (num_Gen_Ds == 1 and num_Gen_Kp == 1 and num_Gen_Km == 1
                    and num_Gen_pi == 1 and num_Gen_gamma == 1):
                if flavour == "mu":
                    nEventsMuGen += 1
                else:
                    nEventsEGen += 1

with open(OUTFILE, "w") as f:
    f.write("{} {} {} {} {}\n".format(nEventsTotal, nEventsMu, nEventsMuGen,
                                      nEventsE, nEventsEGen))

print("{} JOBID={} files={}..{}  total={}  mu={} mu_gen={}  e={} e_gen={}".format(
    CHANNEL, JOBID, start, min(end - 1, NFILES),
    nEventsTotal, nEventsMu, nEventsMuGen, nEventsE, nEventsEGen))
