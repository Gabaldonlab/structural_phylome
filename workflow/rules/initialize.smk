# copied from https://github.com/moiexpositoalonsolab/grenepipe

# import pandas as pd
import os, sys, pwd, re
import socket, platform
# import subprocess
from datetime import datetime
import logging

from snakemake_interface_executor_plugins.settings import ExecMode

# Ensure min Snakemake version
snakemake.utils.min_version("8.11")
basedir = workflow.basedir

# We are currently setting up our own extra log file, so that the below banner is shown.
# Snakemake currently only activates logging to the `.snakemake/log` files _after_ having
# processed all snakefiles, which is not really how logging should work...
# See https://github.com/snakemake/snakemake/issues/2974 for the issue.
# We need to distinguish between the main instance, and the instances of each rule job.
# if logger.mode == ExecMode.DEFAULT:
#     extra_logdir = "snakemake"
# else:
#     extra_logdir = "snakemake-jobs"
# os.makedirs(os.path.join("logs", extra_logdir), exist_ok=True)
# extra_logfile = os.path.abspath(
#     os.path.join(
#         "logs",
#         extra_logdir,
#         datetime.now().isoformat().replace(":", "") + ".log",
#     )
# )
# logger.logger.addHandler(logging.FileHandler(extra_logfile))


# After changing to our new scheme, we can verify the scheme to fit our expextation.
# snakemake.utils.validate(config, schema="../schemas/config.schema.yaml")

# Get a nicely formatted username and hostname
username = pwd.getpwuid(os.getuid())[0]
hostname = socket.gethostname()
hostname = hostname + ("; " + platform.node() if platform.node() != socket.gethostname() else "")

indent = 24

# Get some info on the platform and OS
pltfrm = platform.platform() + "\n" + (" " * indent) + platform.version()
try:
    # Not available in all versions, so we need to catch this
    ld = platform.linux_distribution()
    if len(ld):
        pltfrm += "\n" + (" " * indent) + ld
    del ld
except:
    pass
try:
    # Mac OS version comes back as a nested tuple?!
    # Need to merge the tuples...
    def merge_tuple(x, bases=(tuple, list)):
        for e in x:
            if type(e) in bases:
                for e in merge_tuple(e, bases):
                    yield e
            else:
                yield e

    mv = " ".join(merge_tuple(platform.mac_ver()))
    if not mv.isspace():
        pltfrm += "\n" + (" " * indent) + mv
    del mv, merge_tuple
except:
    pass


# Get the conda env name, if available.
# See https://stackoverflow.com/a/42660674/4184258
conda_env = os.environ["CONDA_DEFAULT_ENV"] + " (" + os.environ["CONDA_PREFIX"] + ")"
if conda_env == " ()":
    conda_env = "n/a"

# Get nicely wrapped command line
cmdline = sys.argv[0]
for i in range(1, len(sys.argv)):
    if sys.argv[i].startswith("--"):
        cmdline += "\n" + (" " * indent) + sys.argv[i]
    else:
        cmdline += " " + sys.argv[i]

# Get abs paths of all config files
cfgfiles = []
for cfg in workflow.configfiles:
    cfgfiles.append(os.path.abspath(cfg))
cfgfiles = "\n                        ".join(cfgfiles)

# Main grenepipe header, helping with debugging etc for user issues
logger.info("=====================================================================================")
logger.info("                               Structural phylome analysis                           ")
logger.info("")
logger.info("    Date:               " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
logger.info("    Platform:           " + pltfrm)
logger.info("    Host:               " + hostname)
logger.info("    User:               " + username)
logger.info("    Python:             " + str(sys.version.split(" ")[0]))
logger.info("    Snakemake:          " + str(snakemake.__version__))
logger.info("    Command:            " + cmdline)
logger.info("")
logger.info("    Base directory:     " + workflow.basedir)
logger.info("    Working directory:  " + os.getcwd())
logger.info("    Config file(s):     " + cfgfiles)
logger.info("")
logger.info("=====================================================================================")
logger.info("")

# No need to have these output vars available in the rest of the snakefiles
del indent
del pltfrm, hostname, username
del cmdline, cfgfiles
