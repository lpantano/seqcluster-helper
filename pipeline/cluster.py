import os
import logger
from cluster_helper import cluster as ipc


resources = {"adapter": [4, 2],
             "align": [45, 8],
             "qc": [8, 1],
             "annotate": [16, 1],
             "report": [8, 1]}


def get_cluster_view(args):
    if not os.path.exists("ipython"):
        os.mkdir("ipython")
        os.mkdir("checkpoint")
    return ipc.cluster_view(args.scheduler, args.queue,
                          args.num_jobs, args.cores_per_job,
                          start_wait=args.timeout,
                          profile="ipython",
                          extra_params={"resources": args.resources,
                                        "mem": args.memory_per_job,
                                        "tag": "ts",
                                        "run_local": args.local})


def wait_until_complete(jobs):
    return [j.get() for j in jobs]


def is_done(step):
    if os.path.exists(os.path.join("checkpoint", step)):
        return True
    return False


def flag_done(step):
    with open(os.path.join("checkpoint", step), "w") as handle:
        handle.write("done")


def send_job(fn, data, args, step):
    """decide if send jobs with ipython or run locally"""
    res = []
    logger.my_logger.debug("doing %s" % step)
    if step not in resources:
        raise ValueError("step not in resources %s" % step)
    else:
        args.memory_per_job = resources[step][0]
        args.cores_per_job = resources[step][1]
    if args.parallel == "ipython":
        if not is_done(step):
            with get_cluster_view(args) as view:
                for sample in data:
                    res.append(view.apply_async(fn, sample, args))
                res = wait_until_complete(res)
            flag_done(step)
            return res
    for sample in data:
        res.append(fn(sample, args))
    return res
