"""Output-only audit status for gate records (lane R, review response). It adds fields; it never changes a decision.

gate/gate.py is hash-frozen (scripts/check_frozen.py), so the status is computed here, from fields the frozen gate
already writes, and attached next to them. `pass_` and `reasons` are never read for anything but the mapping below and
are never modified.

status (one of):
  pass        the frozen gate passed the claim (pass_ is True).
  untestable  the claim failed, and the state offers no downstream transcriptional readout for it: none of the pair's
              listed target genes is measured, so the gate fell back to the receptor's own expression
              (record["response"] == {"receptor": R}), or the ligand/receptor is not measured at all
              (reasons == ["gene_not_in_state"]). A failure here is "no evidence either way", not evidence against.
  fail        the claim failed although a target readout was measured: the measured readout did not show the
              exposure-response relation the claim implies (or a claimed cell type is absent from the section).
              This is evidence against the claim AS READ THROUGH ITS LISTED TARGETS; it is not proof that no
              signalling occurs (a mechanism that bypasses the listed targets, as in the H4 variant A4, also fails).
readout (one of): targets | receptor_only | none   -- what the gate measured as the response.
"""

STATUSES = ("pass", "fail", "untestable")


def readout_of(rec):
    resp = rec.get("response") or {}
    if "targets" in resp:
        return "targets"
    if "receptor" in resp:
        return "receptor_only"
    return "none"


def status_of(rec):
    if rec.get("pass_") is True:
        return "pass"
    ro = readout_of(rec)
    if ro == "receptor_only":
        return "untestable"
    if ro == "none" and "gene_not_in_state" in (rec.get("reasons") or []):
        return "untestable"
    return "fail"


NOTES = {"pass": "state supports the claim through the measured readout",
         "fail": "measured readout does not support the claim",
         "untestable": "no downstream readout measured: no evidence either way"}


def annotate(rec):
    """Return a NEW dict: the record unchanged plus status, readout and status_note (existing keys untouched)."""
    out = dict(rec)
    for k in ("status", "readout", "status_note"):
        assert k not in rec, f"record already has {k!r}"
    s = status_of(rec)
    out["status"], out["readout"], out["status_note"] = s, readout_of(rec), NOTES[s]
    return out


def write_audit(recs, path):
    """Drop-in for gate.gate.write_audit that writes status-annotated records."""
    import json
    with open(path, "w") as f:
        for r in recs:
            f.write(json.dumps(annotate(r), default=float) + "\n")
