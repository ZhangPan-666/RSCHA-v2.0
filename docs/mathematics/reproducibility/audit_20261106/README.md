# Independent checks: 2026-11-06

Install `numpy`, `scipy>=1.15`, and `mpmath`, then run each Python script in
this directory. Each writes the JSON result next to itself. The three scripts
contain independent finite-parameter checks and explicit counterexamples;
they do not import the old handbook verification script and do not prove
every function-space statement. Exact parameter samples and tolerances are
in the code. `structure.json` and `visual_qa.json` record build-review results.
The analytical review and limits are in `notes/handbook_audit_20261106.md`.
