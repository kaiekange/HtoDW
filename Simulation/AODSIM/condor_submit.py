import sys
import re
import subprocess
import os

def list_jobids(server, path):
    # 用 -l 拿文件大小，0 字节的文件（之前抓到过：job 空跑失败但仍以 exit 0 收场）
    # 不算「已完成」，否则会被永久跳过，或者被当成有效输入喂给下一步。
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
    # AODSIM 和 PAT 共用同一个 mode 字符串（比如 ggHtoDsW_phipi），按 mode 匹配会把
    # 正在跑的 AODSIM job 也误判成「PAT 正在跑」，outdir 才是每个步骤唯一的。
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

def submit_for_mode(MODE, ver, YEAR, CONDITIONS, NJOBS, timestamp, server):
    indir = f"/store/user/kakang/Analysis/{MODE}/{ver}/Simulation/LHEGS/{YEAR}/{timestamp}/0000"
    outdir = f"/store/user/kakang/Analysis/{MODE}/{ver}/Simulation/AODSIM/{YEAR}"

    os.makedirs("./tmp", exist_ok=True)

    subprocess.run(["xrdfs", server, "mkdir", "-p", f"{outdir}/log"])
    subprocess.run(["./mkconfig.sh", MODE, CONDITIONS])

    # 只对 LHEGS 已经真正落地、且 AODSIM 还没跑过/没在跑的 JOBID 提交，
    # 避免 LHEGS 还没跑完时对着不存在的输入文件空跑，也避免重复处理已完成/正在跑的 AODSIM。
    expected = set(range(1, NJOBS + 1))
    lhegs_done = list_jobids(server, indir) & expected
    aod_done = list_jobids(server, outdir) & expected
    in_flight = running_jobids(outdir) & expected
    todo = sorted(lhegs_done - aod_done - in_flight)

    print(f"[{MODE}] LHEGS 已完成: {len(lhegs_done)}/{NJOBS}  缺失: {sorted(expected - lhegs_done)}")
    print(f"[{MODE}] AODSIM 已完成(非空): {len(aod_done)}  正在跑: {len(in_flight)}  本次待提交: {len(todo)}")

    if not todo:
        print(f"[{MODE}] 没有需要提交的 job。")
        return

    joblist = f"./tmp/joblist_{MODE}.txt"
    with open(joblist, "w") as f:
        f.write("\n".join(str(j) for j in todo) + "\n")

    scheduler_log = f"condor_logs/HTCondor_{MODE}.log"
    if os.path.exists(scheduler_log):
        os.remove(scheduler_log)

    submit_description = f"""
executable = run.sh
arguments = {MODE} {indir} {outdir} $(JOBID)
log = {scheduler_log}
JobBatchName = MC_{MODE}_AODSIM_{YEAR}_{ver}
request_cpus = 1
request_memory = 8G
request_disk = 10M
+JobFlavour = "tomorrow"
notify_user = kai.kang@cern.ch
notification = error
max_retries = 1
should_transfer_files = NO
queue JOBID from {joblist}
"""
    subprocess.run(["condor_submit"], input=submit_description.encode())

def submit_jobs():

    ver = 20260717
    YEAR = "2022"
    CONDITIONS = "124X_mcRun3_2022_realistic_v12"
    # YEAR="2022EE"
    # CONDITIONS="124X_mcRun3_2022_realistic_postEE_v1"

    # LHEGS crab 输出的时间戳子目录和实际 job 数，两者都跟 MODE 绑定。
    # 每次运行会自动处理下面列出的所有 channel，不需要手动切换。
    # 现在只跑 DsW，DsstW 先不管——以后要加回来，取消下面这行注释即可。
    LHEGS_INFO = {
        "ggHtoDsW_phipi":   dict(timestamp="260720_124856", njobs=500),
        # "ggHtoDsstW_phipi": dict(timestamp="260720_130703", njobs=1000),
    }

    server = "maite.iihe.ac.be"

    for MODE, info in LHEGS_INFO.items():
        submit_for_mode(MODE, ver, YEAR, CONDITIONS, info["njobs"], info["timestamp"], server)

if __name__ == "__main__":
    sys.exit(submit_jobs())
