from argparse import ArgumentParser
# from bcbio.install import get_cloudbiolinux, _default_deploy_args
from bcbio.pipeline import config_utils


# checl installed tools
def check_requirements():
    for tool in ['STAR', 'seqcluster', 'samtools',
                 'bedtools']:
        try:
            config_utils.get_program(tool, {})
            print "Tool installed: %s" % tool
        except Exception, e:
            print "Error: Tool not found: %s" % tool
            raise(e)
    print "Everything is fine."


if __name__ == "__main__":
    parser = ArgumentParser(description="Check installation")
    parser.add_argument("--check", action="store_true",
                        default=False, help="check")

    args = parser.parse_args()
    if args.check:
        check_requirements()


# install bcbio


# install tools
def upgrade_thirdparty_tools(args, remotes):
    """Install and update third party tools used in the pipeline.

    Creates a manifest directory with installed programs on the system.
    """
    s = {"fabricrc_overrides": {"system_install": args.tooldir,
                                "local_install": os.path.join(args.tooldir, "local_install"),
                                "distribution": args.distribution,
                                "use_sudo": args.sudo,
                                "edition": "minimal"}}
    s = _default_deploy_args(args)
    s["flavor"] = "seqcluster_flavor",
    s["actions"] = ["install_biolinux"]
    s["fabricrc_overrides"]["system_install"] = args.tooldir
    s["fabricrc_overrides"]["local_install"] = os.path.join(args.tooldir, "local_install")
    cbl = get_cloudbiolinux(remotes)
    sys.path.insert(0, cbl["dir"])
    cbl_deploy = __import__("cloudbio.deploy", fromlist=["deploy"])
    cbl_deploy.deploy(s)
    manifest_dir = os.path.join(_get_data_dir(), "manifest")
    print("Creating manifest of installed packages in %s" % manifest_dir)
    cbl_manifest = __import__("cloudbio.manifest", fromlist=["manifest"])
    if os.path.exists(manifest_dir):
        for fname in os.listdir(manifest_dir):
            if not fname.startswith("toolplus"):
                os.remove(os.path.join(manifest_dir, fname))
    cbl_manifest.create(manifest_dir, args.tooldir)

