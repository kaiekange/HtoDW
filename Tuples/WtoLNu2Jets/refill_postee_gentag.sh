#!/bin/bash
# postEE GenTag 补交：0J 缺的 500-662，加两个墙钟空壳(0J 181, 1J 81)。
#
# 截断文件故意不补(0J 16 个 / 1J 3 个 / 2J 2 个)：普通版生产两周前撞在几乎同一批
# JOBID 上、同样报 FileReadError，说明是输入 MiniAOD 源头的坏 basket，重交必然再失败。
# 截断文件本身归一化不偏(sumGenWeights_ 逐文件自洽)，只损失统计量。
#
# 分批 + 批间等待，因为 967 个一次性交会让 AAA redirector 超时：09-03 实测 47.5% 的
# job 报 FileOpenError，而 100 个一批时是 0%。
#
# 进度写在 condor_logs/(AFS)，不写 /tmp —— 上一版脚本连同日志被 /tmp 清理杀掉，
# 三天没被发现，最后两批 163 个 job 从来没交成。
# 要中途停止：kill $(cat condor_logs/refill_postee.pid) 之后必须 pgrep 确认，
# TaskStop 只杀 harness 包装层，杀不掉这一层(09-04 因此把 JOBID 300-399 交了两次)。
cd /afs/cern.ch/work/k/kakang/H2DW/CMSSW_15_0_5/src/Tuples/WtoLNu2Jets || exit 1
mkdir -p condor_logs
PIDFILE=condor_logs/refill_postee.pid
REF=condor_logs/refill_postee.ref
D0=/eos/home-k/kakang/H2DW/Analysis/MC/2022EE/WtoLNu2Jets_0J_H_phipimu_GenTag/log
PAUSE=720
echo $$ > $PIDFILE

# 判据是 log 里的 Begin processing / FileOpenError，且只认比本批提交时刻更新的 log。
# 不看输出文件大小：cmsRun 结束时才写输出，运行中的 job 天然是空壳(09-03 差点因此
# 把 88% 正常的一批当成全失败)；也不看旧 log(09-03 读到前一天的残留报出假的 10/10)。
wait_batch() {  # $1=log目录  $2..=JOBID列表
  local d=$1; shift; local ids="$*"; local n=$(echo $ids | wc -w)
  for i in $(seq 1 90); do
    sleep 60
    ok=0; bad=0; pend=0
    for j in $ids; do
      f=$d/run_$j.log
      if [ -f "$f" ] && [ "$f" -nt "$REF" ]; then
        if   grep -q "Begin processing" "$f" 2>/dev/null; then ok=$((ok+1))
        elif grep -q "FileOpenError"   "$f" 2>/dev/null; then bad=$((bad+1))
        else pend=$((pend+1)); fi
      else pend=$((pend+1)); fi
    done
    [ $((ok+bad)) -ge $n ] && break
  done
  echo "  -> 正常 $ok / FileOpenError $bad / 未定 $pend (等 ${i} 分钟)"
}

echo "=========== 补交开始 @ $(date '+%m-%d %H:%M:%S')  PID=$$ ==========="
touch $REF; sleep 1

echo "--- 批次 A: 0J JOBID 181 + 500-599 (101 个) @ $(date '+%H:%M:%S') ---"
python3 condor_submit.py 2022EE 0J --gentag --jobs 181,$(seq -s, 500 599) 2>&1 | grep -E "submitted to cluster"
echo "--- 同批提交: 1J JOBID 81 (1 个) ---"
python3 condor_submit.py 2022EE 1J --gentag --jobs 81 2>&1 | grep -E "submitted to cluster"
wait_batch $D0 181 $(seq 500 599)
if [ "$bad" -gt 10 ]; then echo "!! 批次 A 失败 $bad 个 (>10%)，停止"; rm -f $PIDFILE; exit 1; fi

echo "--- 间隔等待 ${PAUSE}s @ $(date '+%H:%M:%S') ---"
sleep $PAUSE
touch $REF; sleep 1

echo "--- 批次 B: 0J JOBID 600-662 (63 个) @ $(date '+%H:%M:%S') ---"
python3 condor_submit.py 2022EE 0J --gentag --jobs $(seq -s, 600 662) 2>&1 | grep -E "submitted to cluster"
wait_batch $D0 $(seq 600 662)

echo "=========== 提交完毕 @ $(date '+%m-%d %H:%M:%S') ==========="
condor_q -af JobBatchName JobStatus 2>&1 | grep GenTag | sort | uniq -c
rm -f $PIDFILE
