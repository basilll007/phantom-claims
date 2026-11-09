"""Base-Rate Law Calculator for Cell-Cell Communication Candidate Lists (Lane S, S2).

Closed-form law:
    False-Pass Share = (sum_{b in Neg} pi_b * f_b) / (pi_T * r_T + sum_{b in Neg} pi_b * f_b)

where:
    pi_b : composition fraction of negative class b in candidate list
    f_b  : false-pass rate of negative class b under the verification filter
    pi_T : composition fraction of true claims in candidate list
    r_T  : retention / pass rate of true claims under the verification filter

Computes expected false-pass share and 95% CI via Monte Carlo simulation / delta method.
"""
import argparse, json, sys
import numpy as np


def compute_false_pass_share(pi_dict, rate_dict, rate_ci_dict=None, n_mc=10000, seed=42):
    """
    pi_dict: dict mapping class -> candidate fraction (must sum to 1.0)
             classes: 'true', 'niche', 'decoy_sender', 'phantom', 'random'
    rate_dict: dict mapping class -> pass rate
    rate_ci_dict: optional dict mapping class -> [lo_95, hi_95]
    """
    # Normalize pi
    total_pi = sum(pi_dict.values())
    assert abs(total_pi - 1.0) < 1e-4, f"Fractions must sum to 1.0, got {total_pi}"
    
    r_T = rate_dict.get("true", 0.0)
    pi_T = pi_dict.get("true", 0.0)
    
    neg_classes = [c for c in pi_dict if c != "true"]
    
    numerator = sum(pi_dict[c] * rate_dict.get(c, 0.0) for c in neg_classes)
    denominator = pi_T * r_T + numerator
    
    expected_fps = float(numerator / denominator) if denominator > 0 else 0.0
    
    # If CI bounds provided, compute Monte Carlo CI
    ci = None
    if rate_ci_dict is not None:
        rng = np.random.default_rng(seed)
        sim_rates = {}
        for c in pi_dict:
            mean = rate_dict.get(c, 0.0)
            bounds = rate_ci_dict.get(c, [mean, mean])
            # Estimate standard error from 95% CI (1.96 * SE ~ width / 2)
            se = max((bounds[1] - bounds[0]) / (2 * 1.96), 1e-6)
            sim_rates[c] = np.clip(rng.normal(mean, se, n_mc), 0.0, 1.0)
            
        sim_num = sum(pi_dict[c] * sim_rates[c] for c in neg_classes)
        sim_denom = pi_T * sim_rates["true"] + sim_num
        sim_fps = np.where(sim_denom > 0, sim_num / sim_denom, 0.0)
        ci = [float(np.percentile(sim_fps, 2.5)), float(np.percentile(sim_fps, 97.5))]
        
    return {
        "expected_false_pass_share": expected_fps,
        "ci_95": ci,
        "numerator_false_passes": float(numerator),
        "denominator_total_passes": float(denominator),
        "composition": pi_dict,
        "pass_rates": rate_dict,
    }


def main():
    parser = argparse.ArgumentParser(description="Base-Rate Law Calculator for CCC Candidate Lists")
    parser.add_argument("--json", type=str, help="Path to JSON file with composition and rates")
    parser.add_argument("--pi_true", "--pi-true", dest="pi_true", type=float, default=0.5, help="Fraction of TRUE claims in candidate list")
    parser.add_argument("--pi_niche", "--pi-niche", dest="pi_niche", type=float, default=0.25, help="Fraction of NICHE claims in candidate list")
    parser.add_argument("--pi_decoy", "--pi-decoy", dest="pi_decoy", type=float, default=0.15, help="Fraction of DECOY-SENDER claims")
    parser.add_argument("--pi_phantom", "--pi-phantom", dest="pi_phantom", type=float, default=0.05, help="Fraction of PHANTOM claims")
    parser.add_argument("--pi_random", "--pi-random", dest="pi_random", type=float, default=0.05, help="Fraction of RANDOM claims")
    
    # Rates
    parser.add_argument("--r_true", "--r-true", dest="r_true", type=float, default=0.85, help="Pass rate of TRUE claims (retention)")
    parser.add_argument("--f_niche", "--f-niche", dest="f_niche", type=float, default=0.10, help="False-pass rate of NICHE claims")
    parser.add_argument("--f_decoy", "--f-decoy", dest="f_decoy", type=float, default=0.14, help="False-pass rate of DECOY-SENDER claims")
    parser.add_argument("--f_phantom", "--f-phantom", dest="f_phantom", type=float, default=0.002, help="False-pass rate of PHANTOM claims")
    parser.add_argument("--f_random", "--f-random", dest="f_random", type=float, default=0.001, help="False-pass rate of RANDOM claims")
    
    args = parser.parse_args()
    
    if args.json:
        data = json.load(open(args.json))
        pi = data["composition"]
        rates = data["pass_rates"]
        cis = data.get("pass_rate_cis")
    else:
        # Normalize pi
        raw_pi = {
            "true": args.pi_true,
            "niche": args.pi_niche,
            "decoy_sender": args.pi_decoy,
            "phantom": args.pi_phantom,
            "random": args.pi_random,
        }
        tot = sum(raw_pi.values())
        pi = {k: v / tot for k, v in raw_pi.items()}
        rates = {
            "true": args.r_true,
            "niche": args.f_niche,
            "decoy_sender": args.f_decoy,
            "phantom": args.f_phantom,
            "random": args.f_random,
        }
        cis = None
        
    res = compute_false_pass_share(pi, rates, cis)
    print("\n=== Base-Rate Law Prediction ===")
    print(f"List Composition: {json.dumps(res['composition'], indent=2)}")
    print(f"Per-Class Pass Rates: {json.dumps(res['pass_rates'], indent=2)}")
    print(f"\nExpected False-Pass Share: {res['expected_false_pass_share']:.4f} ({res['expected_false_pass_share']*100:.1f}%)")
    if res["ci_95"]:
        print(f"95% Confidence Interval:  [{res['ci_95'][0]:.4f}, {res['ci_95'][1]:.4f}]")


if __name__ == "__main__":
    main()
