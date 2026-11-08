"""Lane R: the audit status is output-only. It maps existing record fields to pass / fail / untestable and must never
change any field the frozen gate wrote."""
import copy
from gate.status import annotate, status_of, readout_of, STATUSES

T = {"targets": ["TGT_1_a", "TGT_1_b"]}


def rec(pass_, response=None, reasons=None):
    r = dict(claim=dict(sender="A", receiver="B", ligand="L", receptor="R"), pass_=pass_,
             reasons=reasons if reasons is not None else ([] if pass_ else ["split_inconsistent"]))
    if response is not None:
        r["response"] = response
    return r


def test_mapping():
    assert status_of(rec(True, T)) == "pass"
    assert status_of(rec(True, {"receptor": "R"})) == "pass"            # a pass is a pass whatever the readout
    assert status_of(rec(False, T)) == "fail"                           # measured target readout, no support
    assert status_of(rec(False, {"receptor": "R"})) == "untestable"     # receptor-only fallback
    assert status_of(rec(False, None, ["gene_not_in_state"])) == "untestable"
    assert status_of(rec(False, None, ["cell_type_absent"])) == "fail"
    assert readout_of(rec(False, T)) == "targets" and readout_of(rec(False, None)) == "none"


def test_annotate_is_output_only():
    for r in (rec(True, T), rec(False, T), rec(False, {"receptor": "R"}), rec(False, None, ["gene_not_in_state"])):
        before = copy.deepcopy(r)
        a = annotate(r)
        assert r == before                                              # input not mutated
        assert {k: a[k] for k in before} == before                      # every existing field unchanged
        assert set(a) - set(before) == {"status", "readout", "status_note"}
        assert a["status"] in STATUSES
        assert (a["status"] == "pass") == (r["pass_"] is True)
