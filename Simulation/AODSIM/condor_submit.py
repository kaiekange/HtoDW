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

def submit_for_mode(MODE, ver, YEAR, NJOBS, timestamp, server, jobid_range=None):
    indir = f"/store/user/kakang/Analysis/{MODE}/{ver}/Simulation/LHEGS/{YEAR}/{timestamp}/0000"
    outdir = f"/store/user/kakang/Analysis/{MODE}/{ver}/Simulation/AODSIM/{YEAR}"

    os.makedirs("./tmp", exist_ok=True)

    subprocess.run(["xrdfs", server, "mkdir", "-p", f"{outdir}/log"])
    # mkconfig.sh owns the GlobalTag/HLT-menu table; it only needs the era. Passing
    # the GlobalTag from here as well would give the two files a way to disagree.
    subprocess.run(["./mkconfig.sh", MODE, YEAR], check=True)

    # 只对 LHEGS 已经真正落地、且 AODSIM 还没跑过/没在跑的 JOBID 提交，
    # 避免 LHEGS 还没跑完时对着不存在的输入文件空跑，也避免重复处理已完成/正在跑的 AODSIM。
    expected = set(range(1, NJOBS + 1))
    lhegs_done = list_jobids(server, indir) & expected
    aod_done = list_jobids(server, outdir) & expected
    in_flight = running_jobids(outdir) & expected
    todo = sorted(lhegs_done - aod_done - in_flight)

    print(f"[{MODE}] LHEGS 已完成: {len(lhegs_done)}/{NJOBS}  缺失: {sorted(expected - lhegs_done)}")
    print(f"[{MODE}] AODSIM 已完成(非空): {len(aod_done)}  正在跑: {len(in_flight)}  可提交: {len(todo)}")

    # 分批提交：只提交落在 jobid_range 内的。范围外的留到后续批次，
    # 发现机制本身不受影响（下次跑会自动跳过已完成/正在跑的）。
    if jobid_range is not None:
        lo, hi = jobid_range
        todo = [j for j in todo if lo <= j <= hi]
        print(f"[{MODE}] 本批范围 {lo}-{hi}，本次待提交: {len(todo)}")

    if not todo:
        print(f"[{MODE}] 没有需要提交的 job。")
        return

    joblist = f"./tmp/joblist_{MODE}_{YEAR}.txt"
    with open(joblist, "w") as f:
        f.write("\n".join(str(j) for j in todo) + "\n")

    scheduler_log = f"condor_logs/HTCondor_{MODE}_{YEAR}.log"
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

# ver 跟 era 绑定（它进输出路径），所以放在按 YEAR 查的表里，不能当常数写在外面：
# 写在外面就会出现"preEE 的 GlobalTag 挂在 postEE 的路径下"这种不报错的错配。
# GlobalTag / HLT menu 不在这里，由 mkconfig.sh 按同一个 YEAR 查表，避免两处各存一份。
ERA_INFO = {
    "2022":   dict(ver=20260717),
    "2022EE": dict(ver=20260825),
}

# LHEGS crab 输出的时间戳子目录和实际 job 数，跟 (MODE, YEAR) 绑定。
# 每次运行会处理当前 era 下列出的所有 channel。
# jobid_range: 分批提交用，只提交这个 JOBID 区间内的；设成 None 表示不限制（全提交）。
# DsstW preEE 从后往前分批跑（LHEGS 唯一没落地的 JOBID 108 在前面，最后一批再处理）：
#   第一批 (701, 1000) —— 已提交 cluster 5863923
#   第二批 (401, 700)  —— 已提交 cluster 5878697
#   第三批 (1, 400)    —— 已提交 cluster 5898993
#   清理批 None        —— 已提交 cluster 6013266，补交损坏(8)+缺log(92)+从未成功(17)共117个
LHEGS_INFO = {
    ("ggHtoDsW_phipi",   "2022"):   dict(timestamp="260720_124856", njobs=500,  jobid_range=None),
    ("ggHtoDsstW_phipi", "2022"):   dict(timestamp="260720_130703", njobs=1000, jobid_range=None),
    # DsW postEE 分批：
    #   第一批 (301, 500) —— cluster 16058790，189 个（当时 11 个 LHEGS 还没落地）
    #   第二批 (1, 300)   —— cluster 16060285，296 个（当时 4 个 LHEGS 还没落地）
    #   补交   None       —— LHEGS 已 500/500，补上上面两批漏掉的 15 个
    #   DsW postEE 收尾：AODSIM 497 好 + 3 个反复崩溃的（19/153/395）已放弃。
    #   那 3 个截断文件**故意留在原处不删**——它们非零，所以 list_jobids 认为已完成，
    #   todo 算出来是空的。删掉反而会让发现机制把它们当"缺失"重新提交。
    ("ggHtoDsW_phipi",   "2022EE"): dict(timestamp="260825_094550", njobs=500,  jobid_range=None),
    # DsstW postEE：
    #   第一批 (701, 1000) 300 个 —— 已提交 cluster 16074741
    #   其余   None       700 个 —— 用户决定不再分批，一次交完 1-700
    #   （None 时 todo = LHEGS 1000 − 已完成 − 正在跑，所以不会碰到第一批还在跑的那些）
    ("ggHtoDsstW_phipi", "2022EE"): dict(timestamp="260827_125408", njobs=1000, jobid_range=None),
}

# (MODE, YEAR) 在这里 = 那批 AODSIM 已经产完、验证过、并被主动删除（只留 PAT，腾空间）。
# 发现机制只看"文件在不在"，会把"故意删掉"读成"缺失"并重新生成整批——真的发生过一次，
# 一次性提交了 500 个本不该跑的多小时 job。
# 写成 (channel, era) 而不是把整个 channel 注释掉：postEE 的 DsW 必须能正常跑，
# 而 preEE 的 DsW 依然撞墙，且不依赖任何人记得"别取消那行注释"。
PERMANENTLY_DONE = {
    ("ggHtoDsW_phipi", "2022"),
    # postEE 的 DsW 同理：PAT 已产完（497/500），AODSIM 被主动删除腾空间，LHEGS 也只留了一个。
    # 不挡住的话，发现机制会拿"LHEGS 还剩几个、AODSIM 一个都没有"算出 todo 并重新提交——
    # 2026-09-01 真的因此误交了 1 个 job（JOBID 1），当场 condor_rm 撤回。
    ("ggHtoDsW_phipi", "2022EE"),
}

def submit_jobs():

    if len(sys.argv) != 2 or sys.argv[1] not in ERA_INFO:
        sys.exit(f"usage: {sys.argv[0]} <{'|'.join(ERA_INFO)}>")
    YEAR = sys.argv[1]
    ver = ERA_INFO[YEAR]["ver"]

    server = "maite.iihe.ac.be"

    keys = [k for k in LHEGS_INFO if k[1] == YEAR]
    if not keys:
        sys.exit(f"LHEGS_INFO 里没有 {YEAR} 的条目")

    for key in keys:
        MODE = key[0]
        if key in PERMANENTLY_DONE:
            print(f"[{MODE}] {YEAR}: AODSIM 已产完并主动删除，跳过（PERMANENTLY_DONE）")
            continue
        info = LHEGS_INFO[key]
        submit_for_mode(MODE, ver, YEAR, info["njobs"], info["timestamp"], server,
                        jobid_range=info["jobid_range"])

if __name__ == "__main__":
    sys.exit(submit_jobs())
