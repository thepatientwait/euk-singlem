# Project notes for agents

Read these first, in this order, and do not re-explore the repo unless they are insufficient:

1. `MEMORY.md`: current state, decisions in force, gotchas.
2. `TODO.md`: open work.
3. `ARCHITECTURE.md`: which file does what.
4. `HISTORY.md`: what was done and why (read a section only when needed; it is long).

At the end of a working session: append a section to `HISTORY.md`, update the status table and decisions in
`MEMORY.md`, tick off or add items in `TODO.md`, and update `ARCHITECTURE.md` if files moved. Commit the docs.

Conventions: tests are `unittest` (see `ARCHITECTURE.md`); heavy jobs go through the batch queue; never commit
`config/secrets.yaml`, outputs, snapshots or rendered notebooks (all gitignored).

**Running workflows: the user always runs the Snakemake workflows themselves** (from their login shell, with
their profile and envs). Do not launch real Snakemake runs (download, generate, tree) from the sandbox, even
locally or in the background. Dry runs (`-n`) to check the DAG, unit tests, and small scratch trials of single
scripts are fine. When a workflow needs rerunning, say which command to run and what to expect.
