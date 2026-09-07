import sys
import subprocess
import os

BKG = "root://maite.iihe.ac.be//store/user/kakang/Analysis/BKG/phipimu/20260624"
SIG = "/eos/home-k/kakang/H2DW/Analysis/MC/2022"

# Backgrounds: (sample, Lumi_reweight). Input file and in-file directory follow
# from the sample name. Lumi_reweight = equivalent-luminosity weight.
BKG_SAMPLES = [
    ("WtoLNu_2Jets",                  "2.8833E+02"),
    ("TTtoLNu2Q",                     "1.7325E+00"),
    ("WGtoLNuG_1Jets",                "4.0053E+01"),
    ("TQbartoLNu",                    "4.6699E+00"),
    ("TbarQtoLNu",                    "4.5078E+00"),
    ("WWtoLNu2Q",                     "2.9648E-01"),
    ("VBFtoLNu",                      "1.3431E+00"),
    ("TbarWplustoLNu2Q",              "1.2800E+00"),
    ("TWminustoLNu2Q",                "1.2167E+00"),
    ("WZtoLNu2Q",                     "2.9351E-01"),
    ("WWGtoLNu2QG",                   "1.1708E-01"),
    ("TTZ_ZtoQQ_TTtoLNu2Q",           "7.3358E-02"),

    ("WZto3LNu",                      "1.3753E-01"),
    ("WZtoL3Nu",                      "1.0083E+00"),

    ("TTto2L2Nu",                     "5.1889E+00"),
    ("TbarWplusto2L2Nu",              "2.4676E+00"),
    ("TWminusto2L2Nu",                "2.4011E+00"),
    ("WWto2L2Nu",                     "3.0185E-01"),
    ("ggWWto2L2Nu_LL",                "2.6335E-01"),
    ("ggWWto2L2Nu_LT",                "8.9708E-02"),
    ("ggWWto2L2Nu_TL",                "9.4232E-02"),
    ("ggWWto2L2Nu_TT",                "3.3844E+00"),
    ("ZZto2L2Nu",                     "1.0748E-02"),

    ("DYto2L_MLL4to10",               "4.5483E+02"),
    ("DYto2L_MLL10to50",              "4.6451E+01"),
    ("DYto2L_MLL50",                  "7.7605E+00"),
    ("DYGto2LG_MLL4to50_PTG10to100",  "3.9076E-01"),
    ("DYGto2LG_MLL4to50_PTG100to200", "1.2644E-02"),
    ("DYGto2LG_MLL4to50_PTG200",      "3.3039E-03"),
    ("DYGto2LG_MLL50",                "1.6773E+00"),
    ("VBFto2L_MLL50",                 "7.1390E-01"),
    ("WZto2L2Q",                      "2.9231E-01"),
    ("ZZto2L2Q",                      "7.4613E-02"),
    ("ZGto2LG",                       "3.7396E-01"),
]

# (sample, Lumi_reweight, input file, in-file directory)
SAMPLES = [(name, weight, f"{BKG}/{name}.root", "Hphipimu") for name, weight in BKG_SAMPLES] + [
    # _2022 = Run3Summer22 (preEE); both signal productions are preEE only.
    ("HtoDsW_2022",   "6.3308E-06", f"{SIG}/Signal_H_phipimu/Signal_H_phipimu.root",             "Hphipimu"),
    ("HtoDsstW_2022", "1.0781E-05", f"{SIG}/Signal_DsstW_H_phipimu/Signal_DsstW_H_phipimu.root", "Hphipimu"),
]

def submit_jobs(only):

    outdir="/eos/home-k/kakang/H2DW/Analysis/MVA/Samples"

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(outdir+"/log", exist_ok=True)

    # Scheduler logs go in condor_logs/ rather than next to the script: one per
    # (era, bin, sign) piles up fast in the working directory.
    os.makedirs("condor_logs", exist_ok=True)
    scheduler_log = f"condor_logs/HTCondor.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    for sample, weight, infile, subdir in SAMPLES:

        if only and sample not in only:
            continue

        submit_description = f"""
executable = run.sh
arguments = {sample} {weight} {infile} {subdir}
log = {scheduler_log}
JobBatchName = Samples_phipimu
request_cpus = 1
request_memory = 8G
request_disk = 100M
+JobFlavour = "workday"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue 1
"""
        subprocess.run(["condor_submit"], input=submit_description.encode())

if __name__ == "__main__":
    sys.exit(submit_jobs(sys.argv[1:]))
