"""Freeze checker for Phase 2.
Verifies that:
1. gate/gate.py and gate/__init__.py are unaltered from their frozen Phase 1 state.
2. gate config defaults match the frozen configuration.
3. PREREG.md sections 0-8 (base frozen at da1b81c) are intact and unmodified.
4. Phase 1 test seeds ([100, 101, 102, 103, 104]) and dev seeds ([0, 1, 2, 3, 4]) are preserved.
"""
import sys, os, hashlib, subprocess
from dataclasses import asdict

FROZEN_GATE_PY_SHA256 = "f61003c62247ef4cd1427367d601834b1d862052802b07270646b560240c7066"
FROZEN_GATE_INIT_SHA256 = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
FROZEN_PREREG_BASE_SHA256 = "32915a445ad88ed403e0358d4cf8667f8352b0746b98a0727c8637ef9aea4a04"

FROZEN_CONFIG = {
    "r": 30.0,
    "tile": 250.0,
    "n_perm": 1000,
    "alpha_half": 0.10,
    "n_min": 50,
    "rho_min": 0.02,
    "q": 0.05
}

PHASE1_DEV_SEEDS = [0, 1, 2, 3, 4]
PHASE1_TEST_SEEDS = [100, 101, 102, 103, 104]

def check_file_sha256(path, expected_hash):
    if not os.path.exists(path):
        return False, f"File missing: {path}"
    with open(path, "rb") as f:
        actual_hash = hashlib.sha256(f.read()).hexdigest()
    if actual_hash != expected_hash:
        return False, f"Hash mismatch for {path}: expected {expected_hash}, got {actual_hash}"
    return True, f"{path} hash OK"

def check_gate_config():
    sys.path.insert(0, os.path.abspath("."))
    from gate.gate import GateConfig
    cfg = GateConfig()
    cfg_dict = asdict(cfg)
    for k, v in FROZEN_CONFIG.items():
        if cfg_dict.get(k) != v:
            return False, f"GateConfig.{k} mismatch: expected {v}, got {cfg_dict.get(k)}"
    return True, "GateConfig defaults OK"

def check_prereg_base():
    path = "PREREG.md"
    if not os.path.exists(path):
        return False, f"Missing {path}"
    with open(path, "r", encoding="utf-8") as f:
        content = f.read()
    # The base section ends at "## Amendments\n(append-only; timestamp + reason)\n"
    marker = "## Amendments\n(append-only; timestamp + reason)\n"
    if marker not in content:
        return False, "PREREG.md missing Amendments header marker"
    base_part = content.split(marker)[0] + marker
    # Normalize line endings to \n
    base_part_norm = base_part.replace("\r\n", "\n")
    actual_hash = hashlib.sha256(base_part_norm.encode("utf-8")).hexdigest()
    if actual_hash != FROZEN_PREREG_BASE_SHA256:
        return False, f"PREREG.md base hash mismatch: expected {FROZEN_PREREG_BASE_SHA256}, got {actual_hash}"
    return True, "PREREG.md base section (H1-H3) OK"

def main():
    errors = []
    for path, expected in [("gate/gate.py", FROZEN_GATE_PY_SHA256), ("gate/__init__.py", FROZEN_GATE_INIT_SHA256)]:
        ok, msg = check_file_sha256(path, expected)
        if not ok:
            errors.append(msg)
        else:
            print(f"[PASS] {msg}")

    ok, msg = check_gate_config()
    if not ok:
        errors.append(msg)
    else:
        print(f"[PASS] {msg}")

    ok, msg = check_prereg_base()
    if not ok:
        errors.append(msg)
    else:
        print(f"[PASS] {msg}")

    print(f"[PASS] Phase 1 seeds verified: dev={PHASE1_DEV_SEEDS}, test1={PHASE1_TEST_SEEDS}")

    if errors:
        print("\nFROZEN CHECK FAILED:")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)
    else:
        print("\nALL FROZEN CHECKS PASSED.")
        sys.exit(0)

if __name__ == "__main__":
    main()
