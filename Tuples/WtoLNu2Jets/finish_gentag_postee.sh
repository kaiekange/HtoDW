#!/bin/bash
# postEE GenTag 收尾流水线：等 265/661 补交 -> 取舍 -> 重新 hadd 0J -> 验证
#                          -> make_samples (cut+dedup+两个era合并) -> DNN predict
#
# 每步验证通过才走下一步，任何一步失败立即 exit 并说明原因 —— 宁可停在中途让人来看，
# 也不要在一个坏的中间产物上继续跑完整条链。
#
# 停止方式: kill $(cat condor_logs/finish_gentag.pid) 然后 pgrep -af finish_gentag 确认。
set -u
SRC=/afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src
D=/eos/home-k/kakang/H2DW/Analysis/MC/2022EE/WtoLNu2Jets_0J_H_phipimu_GenTag
SAMPLES=/eos/home-k/kakang/H2DW/Analysis/MVA/Samples
CLUSTER=16097844
cd $SRC/Tuples/WtoLNu2Jets || exit 1
echo $$ > condor_logs/finish_gentag.pid

say() { echo ""; echo "=========== $* @ $(date '+%m-%d %H:%M:%S') ==========="; }
die() { echo ""; echo "!!!!!!!!!! 停止: $* @ $(date '+%m-%d %H:%M:%S') !!!!!!!!!!"; rm -f condor_logs/finish_gentag.pid; exit 1; }
PY() { source /afs/cern.ch/work/k/kakang/miniconda3/etc/profile.d/conda.sh; conda activate ds_dnn; python3 "$@"; }

say "1/6 等待 265/661 补交 (cluster $CLUSTER)"
for i in $(seq 1 240); do
  n=$(condor_q $CLUSTER -af ClusterId 2>/dev/null | wc -l)
  [ "$n" = "0" ] && break
  sleep 60
done
[ "$n" = "0" ] || die "等了 4 小时补交仍未结束 (剩 $n 个)"
echo "补交结束 (等了 $i 分钟)"

say "2/6 265/661 取舍: 新旧比事例数，新的不如旧的就还原备份"
# 补交是有下行风险的: run.sh 开头 rm -f 输出，失败的重跑会留下比原来更少的事例。
# 所以先前做了备份，这里按事例数择优 —— 不看文件大小，看 nEventsProcessed_。
for j in 265 661; do
  RES=$(PY - <<EOF
import uproot
def n(p):
    try: return int(uproot.open(p)["Hphipimu/Info"].arrays(["nEventsProcessed_"],library="np")["nEventsProcessed_"].sum())
    except Exception: return -1
new=n("$D/output_$j.root"); old=n("$D/backup_before_refill/output_$j.root")
print(f"{new} {old} {'KEEP' if new>old else 'RESTORE'}")
EOF
)
  set -- $RES
  echo "  JOBID $j: 新=$1 旧=$2 -> $3"
  [ "$3" = "RESTORE" ] && { cp -p $D/backup_before_refill/output_$j.root $D/output_$j.root || die "还原 $j 失败"; echo "    已还原备份"; }
done

say "3/6 重新 hadd 0J postEE"
MERGED=$D/WtoLNu2Jets_0J_GenTag_2022EE.root
( source /cvmfs/cms.cern.ch/cmsset_default.sh && cd $SRC && eval $(scramv1 runtime -sh) && \
  hadd -f -k $MERGED $D/output_*.root ) > condor_logs/hadd_0J_postee.log 2>&1 || die "hadd 失败, 见 condor_logs/hadd_0J_postee.log"
echo "hadd 完成: $(stat -c %s $MERGED) 字节"

say "4/6 验证 hadd: job 数与事例数必须与盘上单文件一致"
PY - <<EOF || die "hadd 验证失败"
import uproot, glob, sys
m=uproot.open("$MERGED")["Hphipimu/Info"]
ms=int(m.arrays(["nEventsProcessed_"],library="np")["nEventsProcessed_"].sum())
tot=0; nf=0
for f in glob.glob("$D/output_*.root"):
    try:
        tot+=int(uproot.open(f)["Hphipimu/Info"].arrays(["nEventsProcessed_"],library="np")["nEventsProcessed_"].sum()); nf+=1
    except Exception: pass
print(f"  合并文件: {m.num_entries} 个 job, sum={ms:,}")
print(f"  盘上单文件: {nf} 个, sum={tot:,}")
if m.num_entries!=663 or nf!=663: sys.exit(f"job 数不是 663")
if ms!=tot: sys.exit(f"事例数不一致: 差 {ms-tot}")
print("  一致 OK")
EOF

say "5/6 make_samples: cut + dedup + 两个 era 合并"
cd $SRC/analysis/macros || die "cd 失败"
PY make_samples.py --force --only WtoLNu2Jets_GenTag || die "make_samples 失败"
PY - <<EOF || die "Samples 产物验证失败"
import uproot, collections, sys
f=uproot.open("$SAMPLES/WtoLNu2Jets_GenTag.root"); t=f["Events"]
a=t.arrays(["era","jet_bin","Lumi_reweight","cand_gen_mode"],library="np")
eras=collections.Counter(a["era"].tolist()); bins=collections.Counter(a["jet_bin"].tolist())
print(f"  rows={t.num_entries:,} branches={len(t.keys())}")
print(f"  era 分布: {dict(sorted(eras.items()))}")
print(f"  jet_bin 分布: {dict(sorted(bins.items()))}")
print(f"  加权产额 sum(Lumi_reweight)={a['Lumi_reweight'].sum():,.0f}")
if set(eras)!={1,2}: sys.exit(f"era 不是 {{1,2}} 而是 {set(eras)}")
if "cand_gen_mode" not in t.keys(): sys.exit("缺 cand_gen_mode")
print("  OK")
EOF

say "6/6 DNN predict"
cd $SRC/analysis/DNN || die "cd 失败"
PY predict.py ssdata --samples WtoLNu2Jets_GenTag || die "predict 失败"

say "全部完成"
rm -f $SRC/Tuples/WtoLNu2Jets/condor_logs/finish_gentag.pid
