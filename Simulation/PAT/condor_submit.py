import sys
import re
import subprocess
import os

def list_jobids(server, path):
    # 用 -l 拿文件大小，0 字节的文件（之前抓到过：job 空跑失败但仍以 exit 0 收场，
    # 甚至连文件都没有）不算「已完成」，否则会被当成有效输入喂给 PAT。
    result = subprocess.run(["xrdfs", server, "ls", "-l", path], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"警告: 无法列出 {path} ({result.stderr.strip()})，当作空目录处理")
    ids = set()
    for line in result.stdout.splitlines():
        m = re.search(r'output_(\d+)\.root$', line)
        if not m:
            continue
        parts = line.split()
        try:
            size = int(parts[3])
        except (IndexError, ValueError):
            size = 0
        if size > 0:
            ids.add(int(m.group(1)))
    return ids

def running_jobids(outdir):
    # 排除 condor 队列里已经在跑/排队、正在写这个 outdir 的 JOBID，避免同一个 JOBID
    # 被重复提交、两个 job 同时写同一个输出文件。按 outdir 而不是 mode 匹配——
    # AODSIM 和 PAT 共用同一个 mode 字符串，按 mode 匹配会把正在跑的 AODSIM job
    # 也误判成「PAT 正在跑」，outdir 才是每个步骤唯一的。
    result = subprocess.run(["condor_q", "-af", "Args"], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"警告: 无法查询 condor_q ({result.stderr.strip()})，当作没有正在跑的 job 处理")
    ids = set()
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 4 and parts[2] == outdir:
            try:
                ids.add(int(parts[-1]))
            except ValueError:
                pass
    return ids

def submit_for_mode(mode, ver, year, NJOBS, server):
    indir = f"/store/user/kakang/Analysis/{mode}/{ver}/Simulation/AODSIM/{year}"
    outdir = f"/store/user/kakang/Analysis/{mode}/{ver}/Simulation/PAT/{year}"

    os.makedirs("./tmp", exist_ok=True)

    subprocess.run(["xrdfs", server, "mkdir", "-p", f"{outdir}/log"])
    subprocess.run(["./mkconfig.sh", mode])

    # 只对 AODSIM 已经真正落地、且 PAT 还没跑过/没在跑的 JOBID 提交，
    # 避免 AODSIM 还没跑完时对着不存在的输入文件空跑，也避免重复处理已完成/正在跑的 PAT。
    expected = set(range(1, NJOBS + 1))
    aod_done = list_jobids(server, indir) & expected
    pat_done = list_jobids(server, outdir) & expected
    in_flight = running_jobids(outdir) & expected
    todo = sorted(aod_done - pat_done - in_flight)

    print(f"[{mode}] AODSIM 已完成(非空): {len(aod_done)}/{NJOBS}  缺失: {sorted(expected - aod_done)}")
    print(f"[{mode}] PAT 已完成: {len(pat_done)}  正在跑: {len(in_flight)}  本次待提交: {len(todo)}")

    if not todo:
        print(f"[{mode}] 没有需要提交的 job。")
        return

    joblist = f"./tmp/joblist_{mode}.txt"
    with open(joblist, "w") as f:
        f.write("\n".join(str(j) for j in todo) + "\n")

    scheduler_log = f"condor_logs/HTCondor_{mode}.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    submit_description = f"""
executable = run.sh
arguments = {mode} {indir} {outdir} $(JOBID)
log = {scheduler_log}
output = condor_logs/{mode}_$(JOBID).out
error = condor_logs/{mode}_$(JOBID).err
JobBatchName = MC_{mode}_PAT_{year}_{ver}
request_cpus = 1
request_memory = 8G
request_disk = 10M
+JobFlavour = "workday"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue JOBID from {joblist}
"""
    subprocess.run(["condor_submit"], input=submit_description.encode())

def submit_jobs():

    ver = 20260717
    year = 2022

    # 每个 channel 的总 job 数，跟 AODSIM(20260717 产线)保持一致。
    # 每次运行会自动处理下面列出的所有 channel，不需要手动切换。
    # 现在只跑 DsW，DsstW 先不管——以后要加回来，取消下面这行注释即可。
    NJOBS_MAP = {
        "ggHtoDsW_phipi":   500,
        # "ggHtoDsstW_phipi": 1000,
    }

    server = "maite.iihe.ac.be"

    for mode, NJOBS in NJOBS_MAP.items():
        submit_for_mode(mode, ver, year, NJOBS, server)

if __name__ == "__main__":
    sys.exit(submit_jobs())
