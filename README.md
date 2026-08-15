# leafblower

High-performance bounded survey-weight calibration: multiple solvers
(`oris`, `oris_soft`, raking, Sinkhorn, Greenkhorn, GREG, logit, Newton-KL)
sharing one C++17 core between R and Python, with `min_weight`/`max_weight`
bounds enforced on every returned weight.

On 100,000 observations across 5 margins, `oris_soft` calibrates in
**0.0427s** wall time at **max margin error 3.35e-05**, honouring the `[0,3]`
weight bound and retaining **67,489.4** effective observations — see
[docs/performance.md](docs/performance.md) for the full methodology, machine
provenance, and one-command reproduction.
