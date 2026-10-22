# log_config.py

import logging
import os
import sys

_REDIRECTED = False


def redirect_to_snakemake_log():
    """Send this process's stdout and stderr to the rule's ``log:`` file.

    Only acts inside a ``script:`` rule that declares a log (``snakemake.log``). The file is
    opened with mode "w", so every (re)run of the job overwrites the previous log. The file
    descriptors themselves are redirected (``dup2``), so logging, ``print`` and the output
    of subprocesses all end up in the file. Safe to call more than once.
    """
    global _REDIRECTED
    if _REDIRECTED:
        return None
    smk = getattr(sys.modules.get("__main__"), "snakemake", None)
    log = getattr(smk, "log", None)
    if not log:
        return None

    path = os.fspath(log[0])
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    handle = open(path, "w", buffering=1)
    for stream in (sys.stdout, sys.stderr):
        stream.flush()
    os.dup2(handle.fileno(), 1)
    os.dup2(handle.fileno(), 2)
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(line_buffering=True)
    _REDIRECTED = True
    return path


def configure_logger(file_path=None):
    redirect_to_snakemake_log()

    if file_path is not None:
        logger_name = os.path.splitext(os.path.basename(file_path))[0]
    else:
        logger_name = __name__
    logger = logging.getLogger(logger_name)

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s [%(name)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    return logger
