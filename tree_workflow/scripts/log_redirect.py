"""Send a script: rule's stdout/stderr to its ``log:`` file (see metapackage_workflow/scripts/log_config.py)."""
import os
import sys

_REDIRECTED = False


def redirect_to_snakemake_log():
    """Redirect fds 1 and 2 to ``snakemake.log[0]`` (mode "w": overwritten on every run)."""
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
