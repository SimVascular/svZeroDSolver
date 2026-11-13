"""Run generated 3D cases on a Slurm cluster (e.g. Yale's Bouchet) over SSH.

The cases are generated locally (run_study.py / run_active_study.py without
--run, or cylinder_3d.py generate), uploaded, run with svMultiPhysics in an
Apptainer image of the Docker image used locally, watched, and fetched back
for the usual local post-processing.

Settings (environment):
  SVZEROD_CLUSTER          SSH host (default "bouchet"); log in once with SSH
                           multiplexing (ControlMaster) so that commands need
                           no interactive (MFA) login
  SVZEROD_CLUSTER_WORK     remote work directory (scratch), required
  SVZEROD_CLUSTER_ACCOUNT  Slurm account (optional)

Commands:
  image                    build WORK/svmp-solver.sif (cluster/build_image.sbatch)
  submit CASE_DIR...       upload the cases and submit one job each
  watch CASE_DIR...        stream the remote solver logs through the checks of
                           cylinder_3d.run; cancel a job on failure
  fetch CASE_DIR...        download results/, svmultiphysics.log, wall_time.txt
  status                   the user's jobs

Remote layout: WORK/<case path relative to this folder>, the image
WORK/svmp-solver.sif, and the svZeroD interface library (coupled cases)
WORK/svzerod/build-docker/src/interface/, bound to /svzerod.
"""
import argparse
import json
import os
import shlex
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cylinder_3d as C3  # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
LIB = os.path.join("build-docker", "src", "interface", "libsvzero_interface.so")


def settings():
    work = os.environ.get("SVZEROD_CLUSTER_WORK")
    if not work:
        raise SystemExit("set SVZEROD_CLUSTER_WORK to the remote work directory")
    return dict(host=os.environ.get("SVZEROD_CLUSTER", "bouchet"), work=work.rstrip("/"),
                account=os.environ.get("SVZEROD_CLUSTER_ACCOUNT"))


def ssh(cmd, check=True):
    """Run a shell command on the cluster (non-interactive)."""
    s = settings()
    r = subprocess.run(["ssh", "-o", "BatchMode=yes", s["host"], cmd], capture_output=True, text=True)
    if check and r.returncode != 0:
        raise RuntimeError(f"ssh {s['host']} {cmd!r}: {r.stderr.strip()}")
    return r.stdout


def rsync(src, dst, *extra):
    subprocess.run(["rsync", "-a", "-e", "ssh -o BatchMode=yes", *extra, src, dst], check=True)


def remote_dir(case_dir):
    rel = os.path.relpath(os.path.abspath(case_dir), HERE)
    if rel.startswith(".."):
        raise ValueError(f"{case_dir} is not inside {HERE}")
    return f"{settings()['work']}/{rel}"


def job_script(case, rdir, np_, time_limit, partition, mem_per_cpu):
    s = settings()
    coupled = "coupling" in case
    bind = f"--bind {s['work']}/svzerod:/svzerod " if coupled else ""
    account = f"#SBATCH --account={s['account']}\n" if s["account"] else ""
    # --cleanenv: the image's MPI must not see the host's Slurm variables (it
    # would try to launch through srun); one node, mpirun inside the image
    return f"""#!/bin/bash
#SBATCH --job-name=svmp-{os.path.basename(rdir)}
#SBATCH --partition={partition}
{account}#SBATCH --nodes=1
#SBATCH --ntasks={np_}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu={mem_per_cpu}
#SBATCH --time={time_limit}
#SBATCH --output=slurm_%j.log
cd {shlex.quote(rdir)}
rm -rf results
start=$(date +%s)
apptainer exec --cleanenv {bind}{s['work']}/svmp-solver.sif \\
    mpirun -n {np_} {C3.EXE_IN_IMAGE} solver.xml > svmultiphysics.log 2>&1
status=$?
echo $(( $(date +%s) - start )) > wall_time.txt
exit $status
"""


def submit(case_dir, np_=16, time_limit="12:00:00", partition="day", mem_per_cpu="4G"):
    """Upload a generated case and submit it; returns the job id."""
    with open(os.path.join(case_dir, "case.json")) as f:
        case = json.load(f)
    rdir = remote_dir(case_dir)
    s = settings()
    ssh(f"mkdir -p {shlex.quote(rdir)}")
    rsync(os.path.abspath(case_dir) + "/", f"{s['host']}:{rdir}/",
          "--exclude", "results/", "--exclude", "post/", "--exclude", "*.log",
          "--exclude", "kinematics.json", "--exclude", "coupled.json", "--exclude", "wall_time.txt")
    if "coupling" in case:
        lib_dir = f"{s['work']}/svzerod/{os.path.dirname(LIB)}"
        ssh(f"mkdir -p {lib_dir}")
        rsync(os.path.join(REPO, LIB), f"{s['host']}:{lib_dir}/")
    script = job_script(case, rdir, np_, time_limit, partition, mem_per_cpu)
    # remove the output of earlier runs: `watch` must not read a stale log
    out = ssh(f"cd {shlex.quote(rdir)} && rm -rf results svmultiphysics.log wall_time.txt && "
              f"cat > run.sbatch << 'EOF'\n{script}EOF\nsbatch --parsable run.sbatch")
    job = out.strip().split(";")[0]
    with open(os.path.join(case_dir, "cluster_job.json"), "w") as f:
        json.dump(dict(job=job, host=s["host"], remote=rdir, np=np_, partition=partition,
                       submitted=time.strftime("%Y-%m-%d %H:%M:%S")), f, indent=1)
    return job


def job_state(job):
    out = ssh(f"sacct -j {job} -X -n -P --format=State", check=False).strip().split("\n")[0]
    return out.split()[0] if out else "UNKNOWN"


ACTIVE = ("PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "UNKNOWN")


def watch(case_dirs, echo_every=25, poll=30.0):
    """Follow submitted jobs: stream their logs through cylinder_3d.LogChecker
    (progress lines prefixed with the case), cancel a job on NaN or a diverging
    residual, report state changes and time steps that end above the Newton
    tolerance (svMultiPhysics goes on; a long job is not cancelled for them,
    also not when a restarted watch reads an old log from its start). Returns
    {case_dir: final Slurm state} when all have ended."""
    jobs = {}
    for d in case_dirs:
        with open(os.path.join(d, "case.json")) as f:
            case = json.load(f)
        with open(os.path.join(d, "cluster_job.json")) as f:
            job = json.load(f)
        name = os.path.relpath(os.path.abspath(d), HERE)
        check = C3.LogChecker(case, echo_every=echo_every, strict=False)
        check.echo = False  # progress printed here, with the case name
        jobs[d] = dict(job=job, case=case, check=check, seen=0, state=None, name=name, last=0)
    final = {}
    while len(final) < len(jobs):
        for d, j in jobs.items():
            if d in final:
                continue
            state = job_state(j["job"]["job"])
            if state != j["state"]:
                print(f"{j['name']}: job {j['job']['job']} {state}", flush=True)
                j["state"] = state
            if state == "PENDING":  # its log does not exist yet
                continue
            log = f"{j['job']['remote']}/svmultiphysics.log"
            text = ssh(f"tail -c +{j['seen'] + 1} {shlex.quote(log)} 2>/dev/null", check=False)
            complete = text[:text.rfind("\n") + 1]  # whole lines only
            j["seen"] += len(complete.encode())
            reason = None
            n_unconverged = len(j["check"].unconverged)
            for line in complete.splitlines():
                reason = j["check"](line)
                if reason:
                    break
            for step, res in j["check"].unconverged[n_unconverged:]:
                print(f"{j['name']}: WARNING time step {step} ended at residual {res:.1e} "
                      f"(Newton tolerance {j['case']['newton_tol']:.0e})", flush=True)
            if reason:
                ssh(f"scancel {j['job']['job']}", check=False)
                print(f"{j['name']}: FAILED (job cancelled): {reason}\n  "
                      + "\n  ".join(j["check"].tail), flush=True)
                final[d] = "CANCELLED"
                continue
            step, n = j["check"].step, j["case"]["n_steps"]
            if step and (step // echo_every > j["last"] // echo_every or step == n) and step != j["last"]:
                print(f"{j['name']}: step {step}/{n}", flush=True)
                j["last"] = step
            if state not in ACTIVE:
                final[d] = state
        if len(final) < len(jobs):
            time.sleep(poll)
    return final


def fetch(case_dir):
    """Download the results of a finished job."""
    with open(os.path.join(case_dir, "cluster_job.json")) as f:
        job = json.load(f)
    src = f"{job['host']}:{job['remote']}/"
    rsync(src, os.path.abspath(case_dir) + "/", "--include", "results/***",
          "--include", "svmultiphysics.log", "--include", "wall_time.txt",
          "--include", "slurm_*.log", "--exclude", "*")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("image")
    p = sub.add_parser("submit")
    p.add_argument("cases", nargs="+")
    p.add_argument("--np", type=int, default=16)
    p.add_argument("--time", default="12:00:00")
    p.add_argument("--partition", default="day")
    p.add_argument("--mem-per-cpu", default="4G")
    p = sub.add_parser("watch")
    p.add_argument("cases", nargs="+")
    p.add_argument("--echo-every", type=int, default=25)
    p = sub.add_parser("fetch")
    p.add_argument("cases", nargs="+")
    sub.add_parser("status")
    a = ap.parse_args()
    sys.stdout.reconfigure(line_buffering=True)
    if a.cmd == "image":
        s = settings()
        ssh(f"mkdir -p {s['work']}")
        rsync(os.path.join(HERE, "cluster", "build_image.sbatch"), f"{s['host']}:{s['work']}/")
        acct = f"--account={s['account']} " if s["account"] else ""
        print(ssh(f"cd {s['work']} && sbatch {acct}build_image.sbatch").strip())
    elif a.cmd == "submit":
        for c in a.cases:
            print(f"{c}: job {submit(c, a.np, a.time, a.partition, a.mem_per_cpu)}")
    elif a.cmd == "watch":
        final = watch(a.cases, a.echo_every)
        print("all jobs ended: " + ", ".join(f"{os.path.relpath(os.path.abspath(c), HERE)} {s}"
                                             for c, s in final.items()), flush=True)
        sys.exit(0 if all(s == "COMPLETED" for s in final.values()) else 1)
    elif a.cmd == "fetch":
        for c in a.cases:
            fetch(c)
            print(f"{c}: fetched")
    elif a.cmd == "status":
        print(ssh("squeue --me -o '%.10i %.12P %.30j %.8T %.10M %.6D %R'"))


if __name__ == "__main__":
    main()
