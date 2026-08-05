# samart-mahler

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21711884.svg)](https://doi.org/10.5281/zenodo.21711884)

Companion code and certificates for three papers by **Huimin Zheng**
(College of Information and Network Engineering, Anhui Science and
Technology University; `zhhm@ahstu.edu.cn`):

1. **Mahler measures at interior CM points: proofs of two conjectures of
   Samart** ([arXiv:2608.02255](https://arxiv.org/abs/2608.02255);
   `paper/formal.tex`, `paper/formal.pdf`) — proves, as
   identities of genuine Mahler measures,
   - Theorem A: `m((x+1/x)(y+1/y)(z+1/z)+1) = 4 L'(g7,0)` (Samart's
     conjecture for `k=1`), and
   - Theorem B: `n2((47±45√−7)/2) = (4/7)(54 M7 + d7)`
     (the conjugate pair of Samart's 2015 table).
2. **Samart's conjecture n4(81)=40M7: the exact CM evaluation and the
   two obstructions. A status report**
   ([arXiv:2608.02265](https://arxiv.org/abs/2608.02265);
   `paper/formal2.tex`, `paper/formal2.pdf`) — proves the exact
   series-side evaluation
   `EK4(τ2) = 40M7`, analyses the two obstructions (astroid critical
   image; wrong sheet of the U-series below `Im τ = 1/√2`), and records
   that the literal identity is numerically refuted:
   `n4(81) = 4.1655349907533676508(5)` versus `40M7 = 4.1068643111560…`.
3. **A certified continuation machine for Mahler measure identities at
   CM points, with twelve new proofs of conjectures of Samart**
   (in preparation; `methods/methods.tex`, `methods/methods.pdf`) —
   axiomatizes the differential-comparison continuation method and
   proves twelve entries of Samart's 2015 Table 6 (`n4` family):
   `s = 8292456±3132675√7`, `s = 3656±2600√2`, `s = −144`,
   `s = 143208±101574√2`,
   `s = 1207368+853632√2+697680√3+493272√6` and
   `s = 1207368+853632√2−697680√3−493272√6` (class number two),
   `s = (−192303−85995√5)/2` (Q(√−15), class number two), and
   `s = −893952±516096√3` (Q(√−21), class number four).
   The paper also states six new Heegner-point conjectures
   (discriminants `D = 11, 19, 27, 43, 67, 163`) with exact
   lattice-derived coefficients and a new ray-class constant type,
   and an umbrella conjecture (Section 9; scripts in
   `methods/gen_conj/`).

## Layout

```
paper/     the two papers (LaTeX sources and PDFs), figure data
cert/      certification scripts (the trusted base of the proofs),
           the frozen machine-readable certificate
           certificate_k064.json, and their complete run logs in
           cert/output/
numerics/  numerical-evidence scripts (cross-checks, diagnostics, and
           the direct torus integrations) and logs in numerics/output/
methods/   the methods paper (LaTeX source and PDF), its nineteen
           certification/verification scripts, and the conjecture
           scripts in methods/gen_conj/
```

## Requirements

- Python 3.12+
- [mpmath](https://mpmath.org) (developed and audited against 1.3.0;
  each script prints the version it runs under)
- numpy (only for `paper/gen_fig_n4.py` and `numerics/n4_81_final.py`)

No other dependencies.

## Quick start

The single most important script — the interval-arithmetic path
certificate for Theorem A:

```
cd cert
python cert2_path_k064.py
```

Expected: `ALL CERT-2 (cut [0,64]) CHECKS PASSED (fully rigorous
interval certificate)`, 99 blocks (16/18/37/15/13 across the five path
pieces), global minimum margin `6.009249e-5`, cap coefficient
`≤ -13.54928`. Compare with the archived log
`cert/output/cert2_path_k064.log`. The run also emits the frozen
machine-readable certificate `certificate_k064.json`; the independent
verifier (no adaptive search — it re-checks the frozen block list only)

```
python cert2_verify.py
```

should print `CERTIFICATE VERIFIED` (log: `cert/output/cert2_verify.log`).

## Script inventory

### cert/ — certification scripts (print PASS/FAIL per check)

| script | certifies | typical runtime |
|---|---|---|
| `cert2_path_k064.py` | path for Theorem A avoids the cut `[0,64]` (99 blocks + Taylor cap); emits `certificate_k064.json` | ~1 min |
| `cert2_verify.py` | independent re-verification of the frozen certificate | seconds |
| `s2_iv.py` | shared interval-evaluation core for `s2` (imported by the two scripts above) | — |
| `cert2_path.py` | companion shorter path (avoided set `[-8,8]`, 33 blocks) | ~1 min |
| `cert0_s2_eq_1.py` | exactness `s2(τ0)=1`: formal q-series identity, sextic factorization, interval locks `\|j(2τ0)+3375\| ≤ 1.8e-91`, `\|λ(2τ0)−λ0\| ≤ 6.3e-95` | ~1 min |
| `cert0_s2_n2pair.py` | exactness `s2(τ')=(47+45√−7)/2` and the conjugate (100 dps locks) | ~1 min |
| `verify_P1.py` | CM evaluation at `τ0`: theta identity to `q^60` (exact integers), `T1, T4` Poisson rows, FE assembly | ~1 min |
| `verify_P1_n2pair.py` | CM evaluation at `τw`: 43 checks (60 dps, worst diff `6.9e-52`) | ~1 min |

### numerics/ — numerical evidence (not part of the trusted base)

| script | purpose | typical runtime |
|---|---|---|
| `mahler_m.py` | true `m(f+1)` to 40 digits by direct torus integration | minutes |
| `lvalue_g7.py` | `L'(g7,0)` to 50 digits; functional-equation self-check | ~1 min |
| `find_tau0.py` | `τ0`, `s2(τ0)`, `EK(τ0)` to 60 digits | ~1 min |
| `samart_ek.py` | EK formula; three-way validation at `k=64` | ~1 min |
| `verify_n2_pair.py` | three-way numerical confirmation of Theorem B (50 dps) | ~1 min |
| `diag_ek_branch.py` | wrong-sheet diagnostic at `τ'` | ~1 min |
| `verify_P1_n4_81.py` | (P1) for the n4 note: `EK4(τ2)=40M7`, 45 checks | ~17 min |
| `diag_n4_astroid.py` | wrong-sheet / astroid diagnostics for the n4 note | minutes |
| `kink_resultant_W.py` + `n4_81_final.py` | exact kink detection and direct integration at `c=3`: `n4(81)=4.1655349907533676508(5)` | ~7.5 h |

### methods/ — the continuation-machine paper

`methods.tex` / `methods.pdf` are the paper; the scripts below are its
certifier (all print per-check PASS/FAIL and a final
`ALL CHECKS PASSED`):

| script | certifies | typical runtime |
|---|---|---|
| `cert0_s4_s7pair.py` | exactness `s4(√−7) = 8292456+3132675√7`, `s4(√−7/2) = 8292456−3132675√7` (Q(√7) pair lock) | seconds |
| `verify_P1_n4_s7pair.py` | Theorem A (`s = 8292456±3132675√7`), three tracks | <1 min |
| `cert0_s4_m144.py` | exactness `s4((1+√−3)/2) = −144`: formal q-series + interval lock | seconds |
| `cert2_path_n4_m144.py` | certified path from the anchor box to `τ1 = (1+√−3)/2` inside `V4` | <1 min |
| `verify_P1_n4_m144.py` | Theorem B (`s = −144`), three tracks | <1 min |
| `n4_m144_true_mahler.py` | true-Mahler side of `n4(−144)` by direct torus integration at `c = 2√3 e^{−iπ/4}` (outside the astroid) | seconds |
| `cert0_n4_p3_t1t2.py` | exactness `s4(√−2) = 3656+2600√2`, `s4((1+√−2)/2) = 3656−2600√2` | seconds |
| `verify_P1_n4_p3_t1t2.py` | Theorems C and D (`s = 3656+2600√2`, `s = 143208+101574√2`), three tracks | <1 min |
| `cert0_n4_p4_t3.py` | exactness of the class-number-two `s4` quartic in Q(√2,√3) | seconds |
| `verify_P1_n4_p4_t3.py` | Theorem E (`s = 1207368+853632√2+697680√3+493272√6`, Q(√−6), h=2), three tracks | <1 min |
| `n4_p4_t4_cert.py` | certified path from the anchor box to `τ4 = (1+√−2)/2` inside `V4` | <1 min |
| `n4_p4_t4_verify.py` | Theorem F (`s = 3656−2600√2`), three tracks | <1 min |
| `n5_line_cert.py` | vertical line `Re τ = 1/2`, `0.702 ≤ Im τ ≤ √21/2` inside the good component `W` (station S3 of Theorems G–K; also used for the status of `s = −3969`) | seconds |
| `cert0_n5_e56.py` | exactness of the `s4` pair in Q(√5): `X² + 192303X + 1185921` (Theorem I and blocked entry #5) | seconds |
| `cert0_n5_e78.py` | exactness of the `s4` pair in Q(√3): `X² + 1787904X + 84934656` (Theorems J, K) | seconds |
| `verify_P1_n5_e1.py` | Theorem G (`s = 143208−101574√2`, Q(i), h=1), three tracks, 49 checks | ~1 min |
| `verify_P1_n5_e2.py` | Theorem H (`s = 1207368+853632√2−697680√3−493272√6`, Q(√−6), h=2), three tracks, 55 checks | ~1 min |
| `verify_P1_n5_e6.py` | Theorem I (`s = (−192303−85995√5)/2`, Q(√−15), h=2), three tracks, 44 checks | ~1 min |
| `verify_P1_n5_e78.py` | Theorems J and K (`s = −893952±516096√3`, Q(√−21), h=4), three tracks, 73 checks | <1 min |

### methods/gen_conj/ — the Heegner-point conjecture scripts

The scripts behind Section 9 of the paper (Conjecture 9.2, the six new
Heegner-point identities; numerical evidence only, not iv-certified):

| script | purpose |
|---|---|
| `gen_conj_fit.py D` | per-discriminant fit machine; control run `D=7` reproduces Samart's `(10/7)(40M7+d7)` with residual 0 |
| `gen_conj_fit4.py D` | shifted B/G-row machine with exact tail correction; derives the coefficients of Conjecture 9.2 (`D = 11, 19, 43, 67, 163`) |
| `gen_conj_fit5.py` | the `D = 27` branch (orders of conductor 3 and 6) |
| `gen_conj_s4.py` | `s4` values and the 40-digit `\|s4\| > 256` ray scan |
| `gen_conj_minpoly.py` | 110-digit PSLQ search for minimal polynomials of the `s4` parameters (none of degree ≤ 4) |

Representative run logs (`gen_conj_fit4_*.log`, `gen_conj_fit5_27.log`,
`gen_conj_minpoly.log`) are included.

## Reproducibility

All bounds entering the proofs are produced by outward-rounded interval
operations (`mpmath.iv`); no machine-float arithmetic occurs in any
certificate path. Re-running the scripts in `cert/` regenerates every
certificate from scratch; the archived logs in `cert/output/` are the
exact outputs on the development machine.

## License

- Code: MIT License (see `LICENSE`).
- Papers (`paper/formal.*`, `paper/formal2.*`): © Huimin Zheng, 2026.
  All rights reserved.

## AI-assistance disclosure

The papers were written with AI assistance (Kimi, Moonshot AI), as
declared on their first pages. All mathematical content, including
every proof and every certified computation, has been checked by the
author, who takes full responsibility for correctness.
