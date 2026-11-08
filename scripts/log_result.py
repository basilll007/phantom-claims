"""Thread- and process-safe append-only logger for results/results.jsonl.
Mandatory fields per Phase 2 spec:
- lane: string (e.g. 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'L')
- agent_id: string
- git_hash: commit hash
- seed_set: 'dev' | 'test1' | 'N2'
- labelled: 'confirmatory' | 'exploratory'
"""
import os, sys, time, json, subprocess
try:
    import msvcrt
except ImportError:
    msvcrt = None
try:
    import fcntl
except ImportError:
    fcntl = None

RESULTS_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "results.jsonl")

def get_git_hash():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"

def log_result(rec, lane=None, agent_id="gemini-lead", seed_set="N2", labelled="confirmatory"):
    record = dict(rec)
    if lane is not None:
        record["lane"] = lane
    elif "lane" not in record:
        record["lane"] = "L"
        
    if "agent_id" not in record:
        record["agent_id"] = agent_id
    if "git_hash" not in record:
        record["git_hash"] = get_git_hash()
    if "seed_set" not in record:
        record["seed_set"] = seed_set
    if "labelled" not in record:
        record["labelled"] = labelled
        
    assert record["seed_set"] in ("dev", "test1", "N2"), f"Invalid seed_set: {record['seed_set']}"
    assert record["labelled"] in ("confirmatory", "exploratory"), f"Invalid labelled: {record['labelled']}"
    assert "lane" in record and record["lane"], "Missing lane in result record"
    assert "agent_id" in record and record["agent_id"], "Missing agent_id in result record"

    if "timestamp" not in record:
        record["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")

    os.makedirs(os.path.dirname(RESULTS_FILE), exist_ok=True)
    line = json.dumps(record, default=float) + "\n"

    # Process/file locking
    with open(RESULTS_FILE, "a", encoding="utf-8") as f:
        if msvcrt:
            # Lock 1 byte at current position
            pos = f.tell()
            msvcrt.locking(f.fileno(), msvcrt.LK_LOCK, 1)
            try:
                f.write(line)
                f.flush()
            finally:
                f.seek(pos)
                msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, 1)
        elif fcntl:
            fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            try:
                f.write(line)
                f.flush()
            finally:
                fcntl.flock(f.fileno(), fcntl.LOCK_UN)
        else:
            f.write(line)
            f.flush()

    return record

if __name__ == "__main__":
    if len(sys.argv) > 1:
        with open(sys.argv[1], "r", encoding="utf-8") as rf:
            data = json.load(rf)
            log_result(data)
            print("Logged 1 record.")
