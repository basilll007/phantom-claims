"""Lane S Orchestrator:
Runs the remaining pipeline components in sequential order to guarantee memory safety:
1. breast_rep2 (if not already completed)
2. lung (all 10 seeds N3)
3. laneS_s3_eval.py (5K richer panels evaluation)
4. laneS_aggregate.py --log (hypothesis evaluation H4, H8, H9, figure generation, logging)
"""
import glob, os, subprocess, sys, time

PYTHON = os.path.abspath(".venv/Scripts/python.exe")
ENV = os.environ.copy()
ENV["PYTHONPATH"] = "."

def run_step(desc, cmd):
    print(f"\n{'='*60}", flush=True)
    print(f"STEP: {desc}", flush=True)
    print(f"CMD:  {' '.join(cmd)}", flush=True)
    print(f"{'='*60}\n", flush=True)
    t0 = time.time()
    res = subprocess.run(cmd, env=ENV)
    if res.returncode != 0:
        print(f"FAILED with return code {res.returncode}", flush=True)
        sys.exit(res.returncode)
    print(f"\n[DONE] {desc} in {time.time() - t0:.1f}s\n", flush=True)

if __name__ == "__main__":
    # Check if breast_rep2 is still running or needs completion
    run_step("Benchmark on breast_rep2 N3", [PYTHON, "-m", "scripts.laneS_run", "run", "breast_rep2", "N3"])
    
    # Run lung N3
    run_step("Benchmark on lung N3", [PYTHON, "-m", "scripts.laneS_run", "run", "lung", "N3"])
    
    # Run S3 Richer Panels (5K) evaluation
    run_step("S3 Richer Panels (5K) evaluation", [PYTHON, "scripts/laneS_s3_eval.py"])
    
    # Run Lane S Aggregation & Hypothesis Evaluation
    run_step("Lane S Aggregation and Hypotheses (H4, H8, H9)", [PYTHON, "-m", "scripts.laneS_aggregate", "--log"])
    
    print("\nALL LANE S EXPERIMENTS AND EVALUATIONS COMPLETED SUCCESSFULLY!", flush=True)

