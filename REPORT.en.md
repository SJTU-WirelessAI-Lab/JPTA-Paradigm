# JPTA 2022 Reproduction Report

[中文](REPORT.md) | English

Reproduced paper: V. V. Ratnam, J. Mo, A. Alammouri, B. L. Ng, J. Zhang, A. F. Molisch,
"Joint Phase-Time Arrays: A Paradigm for Frequency-Dependent Analog Beamforming in 6G",
*IEEE Access*, vol. 10, pp. 73364–73377, 2022.

Code: `matlab/` (MATLAB R2026b, no parallel toolbox, fully vectorized).
One-shot run: `cd matlab; run_all` (self-test + six figures, about 20 min); a single
figure with `run_all('fig4')` and so on.

---

## 1. Summary

| Figure | Content | Outcome |
|---|---|---|
| Fig. 4 | JPTA array-gain heat maps (behaviors 1/2, `N=M`, `κ=M`) | Same structure; behavior 1 gives `Fobj=0.9307` (paper ≈0.93) |
| Fig. 5 | `Fobj`–`N`, 4 curves × 2 behaviors | Matches point by point, max deviation 0.02 (0.03 for the heuristic curve) |
| Fig. 6 | `Fobj`–`κ`, 4 curves × 2 behaviors | Matches point by point; both the `κ=32` and the `κ≈1.6` saturation reproduced |
| Fig. 7 | Convergence percentile band (100 random realizations) | Shape and conclusion agree; the `i=1` percentile is about 0.07 high (sampling distribution assumed) |
| Fig. 9 | HBF-FC/PC array-gain heat maps | All four panels agree, including the two-segment PC beam and its weak mirror lobe |
| Fig. 10 | `Fobj`–`N_RF`, FC/PC with JPTA references | FC curve matches; both JPTA reference lines reproduced exactly |

The two heat-map figures (4 and 9) are qualitative reproductions; the three line-plot
figures (5, 6, 10) align point by point.
Every core claim of the paper is supported: JPTA approaches the ideal beam with a single
RF chain, whereas HBF-FC needs 22 and HBF-PC needs 32 RF chains to reach the same level.

---

## 2. Code structure

```
matlab/
  cfg_default.m          parameter table (Sec. VII), supports cfg_default('N',16,'kappa',32)
  selftest.m             self-test: norms, upper bound, eq. (8a) peak, convergence
  run_all.m              single entry point, dispatches per figure
  core/
    build_map.m          mapping matrix P of eq. (1)
    array_resp.m         array response a_k(θ) of eq. (2)
    gen_beam.m           desired beams {b_k} (behaviors 1/2)
    jpta_forward.m       forward model T*P*d_k (the only "physical" implementation)
    fobj.m               Fobj of eq. (12)
    array_gain.m         array gain |a_k(θ)†x_k|²
    hbf_fc.m             PE-AltMin of eqs. (13)(14) (fully connected)
    hbf_pc.m             alternating minimization under the block-diagonal constraint
  alg/
    solve_tau_ls.m       line search of eq. (8a) (FFT grid + parabolic interpolation)
    solve_tau_wls.m      phase unwrapping + weighted least squares of Lemma 1
    solve_phi.m          eq. (8b)
    update_alpha.m       eq. (11)
    alg1_jpta.m          Algorithm 1 (alternating optimization)
    alg2_heur1.m         Algorithm 2 (behavior-1 heuristic)
    alg3_heur2.m         Algorithm 3 (behavior-2 heuristic)
  figs/
    fig4_arraygain.m     Fig. 4
    fig5_fobj_vs_N.m     Fig. 5
    fig6_fobj_vs_kappa.m Fig. 6
    fig7_convergence.m   Fig. 7
    fig9_hbf_arraygain.m Fig. 9
    fig10_fobj_vs_nrf.m  Fig. 10
  out/                   PNGs and .mat caches (.mat is not tracked)
```

The forward model exists in exactly one place, `jpta_forward.m`, and every algorithm and
every figure calls it. This prevents subtle inconsistencies such as "Fig. 4 and Fig. 5 use
different array-response conventions".

---

## 3. Per-figure results

### 3.1 Fig. 4: array-gain heat maps

Parameters: `N=M=64`, `κ=M=64`, `maxiter=10`.

| Item | Reproduced | Paper |
|---|---|---|
| behavior 1 `Fobj` | **0.9307** | ≈0.93 |
| behavior 2 `Fobj` | **0.8047** | ≈0.80 |

All four panels agree in structure: (a) and (b) sweep the main lobe from 95 GHz/7.5° to
105 GHz/52.5°, (c) shows the step, and (d) can only switch gradually, with sidelobes.
The color scale is **linear 0–64** (=M), not dB.

### 3.2 Fig. 5: `Fobj` vs `N` (`κ=M`)

| N | eq. (8) | Paper | Lemma 1 | Paper | Heuristic | Paper |
|---|---|---|---|---|---|---|
| 1 | 0.1995 | 0.20 | 0.1989 | 0.20 | 0.0936 | 0.08 |
| 2 | 0.2736 | 0.27 | 0.2769 | 0.27 | 0.1571 | 0.15 |
| 4 | 0.3800 | 0.39 | 0.3852 | 0.39 | 0.2648 | 0.27 |
| 8 | 0.5086 | 0.50 | 0.5031 | 0.50 | 0.4263 | 0.43 |
| 16 | 0.7281 | 0.72 | 0.7276 | 0.72 | 0.6866 | 0.69 |
| 32 | 0.8844 | 0.88 | 0.8842 | 0.88 | 0.8150 | 0.82 |
| 64 | 0.9307 | 0.93 | 0.9306 | 0.93 | 0.8514 | 0.855 |

Behavior 2 (Fig. 5b): at `N=1` we get 0.5734/0.5733/0.3969 (paper 0.565/0.555/0.40), and
at `N=64` we get 0.8047/0.8049/0.5319 (paper 0.81/0.81/0.565).

Lemma 1 and eq. (8) coincide almost exactly (difference below 0.001), reproducing the
paper's claim that Lemma 1 cuts the complexity substantially with almost no loss in `Fobj`.

### 3.3 Fig. 6: `Fobj` vs `κ` (`N=M`)

| κ | eq. (8) | Paper | Lemma 1 | Paper | Heuristic | Paper |
|---|---|---|---|---|---|---|
| 1 | 0.2725 | 0.265 | 0.2062 | 0.20 | 0.1064 | 0.10 |
| 2 | 0.3392 | 0.325 | 0.2172 | 0.21 | 0.1392 | 0.135 |
| 4 | 0.4314 | 0.435 | 0.2405 | 0.25 | 0.2150 | 0.22 |
| 8 | 0.5120 | 0.485 | 0.4285 | 0.425 | 0.3818 | 0.385 |
| 16 | 0.7482 | 0.745 | 0.7429 | 0.735 | 0.6815 | 0.685 |
| 32 | 0.9307 | 0.93 | 0.9306 | 0.93 | 0.8514 | 0.855 |
| 64 | 0.9307 | 0.93 | 0.9306 | 0.93 | 0.8514 | 0.855 |

**The curves for `κ=32` and `κ=64` coincide exactly**, so the saturation point lies in
`(16,32]`, consistent with `κ_req≈22.3` derived in README §6.

Behavior 2 (Fig. 6b): at `κ=0.1` we get 0.6232/0.6122/0.4054 (paper 0.615/0.615/0.415);
for `κ≥1.6` the curves flatten at 0.8047/0.8049/0.5393, consistent with `3/W` from Remark 2.

### 3.4 Fig. 7: convergence speed

Mean and 10/90 percentiles of `Fobj(i)/Fobj(30)` over 100 random realizations:

| | behavior 1, i=1 | behavior 1, i=10 | behavior 2, i=1 | behavior 2, i=10 |
|---|---|---|---|---|
| 10th percentile | 0.576 | 0.985 | 0.920 | 0.9999 |
| mean | 0.802 | 0.994 | 0.962 | 0.9994 |
| 90th percentile | 1.000 | 1.000 | 1.000 | 1.000 |
| paper, i=1 (mean / 10th pct) | 0.73 / 0.47 | — | 0.93 / 0.85 | — |

**The conclusion agrees**: `max-iter=10` is sufficient, and after 10 iterations every
percentile is at least 0.985. The `i=1` percentiles are about 0.07 higher than in the
paper, because the paper does not state the distribution of the random parameters; this
implementation draws `N∈{1,2,4,8,16,32,64}`, `κ∈[1,M]`, `θ0∈[-π/3,π/3]`,
`Δθ∈[0,π/2]`, which is a milder distribution, so the first iteration is already good.
This figure therefore matches in shape and conclusion, with absolute values that depend on
the sampling assumption.

### 3.5 Fig. 9: HBF array-gain heat maps

| Panel | Configuration | Reproduced | Paper |
|---|---|---|---|
| (a) | FC, `N_RF=22`, behavior 1 | diagonal beam 7.5°→52.5° | same |
| (b) | FC, `N_RF=2`, behavior 2 | step from −45° to +30° | same |
| (c) | PC, `N_RF=32`, behavior 1 | main diagonal plus a weak mirror lobe below | same |
| (d) | PC, `N_RF=32`, behavior 2 | +15°/−45° at low frequencies, +30°/−30° at high frequencies | same |

Compared with Fig. 4, the beam-pointing transitions of HBF are sharper, which is the
"sharper transition" the paper describes and is visible by eye in the heat maps.

### 3.6 Fig. 10: `Fobj` vs `N_RF`

| `N_RF` | 1 | 2 | 4 | 8 | 16 | 32 | 64 |
|---|---|---|---|---|---|---|---|
| FC, behavior 1 (`θ0=π/6`) | 0.191 | 0.280 | 0.404 | 0.582 | 0.827 | 0.993 | 1.000 |
| paper | 0.18 | 0.30 | 0.42 | 0.60 | 0.80 | 0.99 | 1.0 |
| PC, behavior 1 | 0.187 | 0.286 | 0.394 | 0.554 | 0.767 | 0.948 | 1.000 |
| paper | 0.19 | 0.22 | 0.30 | 0.43 | 0.56 | 0.95 | 1.0 |
| FC, behavior 1 (`θ0=0`) | 0.062 | 0.096 | 0.174 | 0.307 | 0.481 | 0.748 | 1.000 |
| PC, behavior 1 (`θ0=0`) | 0.062 | 0.105 | 0.174 | 0.277 | 0.458 | 0.651 | 1.000 |
| FC, behavior 2 | 0.466 | 0.878 | 0.994 | 0.984 | 0.999 | 1.000 | 1.000 |
| paper | 0.67 | 0.885 | 0.99 | — | — | — | 1.0 |
| PC, behavior 2 | 0.466 | 0.685 | 0.704 | 0.709 | 0.769 | 0.810 | 1.000 |
| paper | 0.62 | 0.692 | 0.70 | 0.71 | 0.77 | 0.805 | 1.0 |

JPTA reference lines (`N=64`, `κ=64`):

| Scenario | Reproduced | Paper (read from figure) |
|---|---|---|
| behavior 1, `θ0=π/6`, `Δθ=π/4` | **0.9307** | 0.93 |
| behavior 1, `θ0=0`, `Δθ=2π/3` | **0.6661** | 0.666 |
| behavior 2 | **0.8047** | 0.805 |

All three reference lines match exactly. The "22 RF chains" of Proposition 1 corresponds to
`N_RF≈24` on the FC curve (0.827 at 16, 0.993 at 32), which agrees with the paper.

**Two deviations need comment**:

1. At `N_RF=1`, both FC and PC degenerate to a purely analog beam. A single beam can serve
   only one of `θ1`, `θ2`, so `Fobj≈0.5` (this implementation gives 0.466). The paper reads
   0.62–0.67, which is above that bound, and the paper itself has FC(1)≠PC(1) — yet the two
   feasible sets are identical at `N_RF=1`, so the two values must be equal. This
   implementation therefore takes the theoretical bound 0.466 and lowers the y-limit of
   Fig. 10(b) from 0.6 to 0.45 so the point is visible.
2. PC behavior 1 is about 0.2 above the paper for `N_RF=2…16`. The paper does not disclose
   its initialization; selecting solutions by the Frobenius residual of eq. (14) alone lands
   on a suboptimal point. This implementation selects by `Fobj` instead (see §5, item 4),
   which moves PC behavior 2 from 0.51 to 0.685 and matches the paper's 0.692, at the cost of
   being better than the paper in the mid range of behavior 1.

---

## 4. Typos and ambiguities in the paper

| Location | As printed | Should be | Why |
|---|---|---|---|
| Algorithm 3, `e_b` | both terms use `e^{jπm sinθ2}` | the second term should use `sinθ1` | otherwise it degenerates to a single beam and the step structure of behavior 2 is lost |
| Sec. VII | `Fobj` saturates after `κ=M sin(Δθ)` | `M sin(Δθ/2)` | Remark 2 gives M\|sin(Δθ/2)\|/W; M sin(Δθ)≈2×24.5, i.e. `Δθ/2` was replaced by `Δθ` |
| Sec. VII | behavior 2 saturates after `κ=2` | `κ=3` | Remark 2 gives `3/W`, the slope of Algorithm 3 corresponds to a span of 3, and the saturation in Fig. 6(b) lies between 2 and 3 |
| Proposition 1 | "22 RF chains" in the text | eq. (15) gives 23 (22.69 rounded) | the narrowband approximation of eq. (19) happens to give 22 |

Two further items are **not typos, but implementing them literally fails**:

- `0 ≤ τ_n ≤ κ/W` in Algorithm 1 should be read as the final hardware range, obtained by a
  centering shift. The line search must use the symmetric interval `[-κ/2W, κ/2W]` of
  eq. (4); otherwise the behavior-1 peak (negative τ) lies outside the interval and `Fobj`
  falls from 0.93 to 0.67 with oscillating iterations.
- The hook symbol in the Algorithm 3 pseudocode is the **phase operator `∠`**, not a
  multiplication or division sign. The PDF text layer extracts it as `̸`.

---

## 5. Implementation pitfalls (in the order they were hit)

1. **The search interval of eq. (8a)**: see above. The delays required by behavior 1 are
   negative (`τ_m≈-(m-1)sinθ0/(2f0)`, about −2.1 ns for m=64), so searching over `[0,κ/W]`
   stops on a sidelobe. Counter-evidence: Lemma 1 uses the symmetric interval of eq. (10a)
   and gives 0.9306 directly for the same parameters.
2. **Eq. (8a) can be evaluated over the whole candidate grid with one FFT**: the carrier
   part of `e^{-j2πf_kτ}` is a common phase independent of m and k and vanishes in the
   magnitude, leaving exactly the FFT bins. Nfft=32768 gives a 6.25 ps grid plus parabolic
   interpolation. The prerequisite is that `f0·K/W=20480` is an integer (the objective is
   periodic in τ with period K/W).
3. **The normal equations of Lemma 1 require frequency centering and normalization first**:
   `x=2π(f−mean(f))/max|·|`. Otherwise `2πf` reaches 6.6e11 while the in-band relative
   variation is only 10%, so slope and intercept are almost collinear, `RCOND≈1e-21`, and
   MATLAB reports a nearly singular matrix.
4. **HBF solutions must not be selected by the Frobenius residual alone**: under the PC
   structure, "each subarray points the same way" and "subarrays share the directions" give
   nearly identical residuals but differ by 0.2 in `Fobj` (0.51 versus 0.69 for behavior 2
   at `N_RF=2`). Selecting by residual alone lands consistently on the poor solution, which
   is exactly why Fig. 10(b) did not match at first. Selecting by `Fobj` fixes it.
5. **The legend order in Figs. 5 and 6 differs from the data-column order**: the data
   columns are `[eqn(8), Lemma1, heuristic]` while the style order is
   `[heuristic, eqn(8), Lemma1]`, so an index map is needed before assigning handles;
   otherwise the legend labels are mismatched with the curves.
6. **The Fig. 6 x-axis needs a decade of headroom**: `10^0–10^2` in (a) and `10^-1–10^1` in
   (b). If `xlim` is taken directly from the data range, the rightmost sample touches the
   frame and the `10^2` tick is not drawn.
7. **MATLAB environment details**: a `-batch` console is GBK, so Chinese text printed with
   `fprintf` is garbled — all console output is ASCII; `TiledChartLayout.Padding` only
   accepts `loose|compact|tight` in R2026b and rejects numeric vectors (use `subplot` when
   you need padding for a colorbar title); `clim` replaces the deprecated `caxis`.

---

## 6. Self-test and invariants

`selftest.m` checks the following on every run:

| Check | Result |
|---|---|
| `‖a_k(θ)‖²=M` | error 0 (exact) |
| `‖TPd_k‖=1` | error 0 (exact) |
| `Fobj` of the ideal beam | 1.0000 (exact) |
| `Fobj≤1` over 50 random parameter sets | all satisfied (Cauchy–Schwarz bound) |
| correctness of the eq. (8a) peak | no TTD improves after a ±80 ps perturbation |
| convergence | plateau after 3 iterations, 10 is enough |

Because `Fobj∈[0,1]`, it doubles as a unit test: any value above 1 is immediately a bug.

---

## 7. Open items and suggestions

**Not done**:
- Fig. 8 is an HBF architecture sketch, Table 1 is analytical and Table 2 is quoted data;
  none of them needs simulation.
- The multi-RF-chain JPTA discussed at the end of Sec. V (listed as future work) is out of
  scope.

**Two observations worth writing up** (obtained during the reproduction, not stated in the
paper):

1. Behavior 1 needs only about **2.2 ns** of TTD (`κ≈22`), whereas the paper uses
   `κ=M=64`, i.e. 6.4 ns — roughly a 3× hardware margin. This matters for hardware cost:
   the maximum-delay specification of a TTD can be designed as `M sin(Δθ/2)/W` instead of
   copying `M/W`.
2. Fig. 7 shows that with `max-iter=3` the 10th percentile already reaches 0.93 for
   behavior 1 and 0.98 for behavior 2. The paper's 10 iterations are conservative; 3–5 are
   enough, which matters for real-time beam tracking.

**To make the reproduction stricter**:
- The random distribution behind Fig. 7 should be confirmed with the authors (it is the
  only figure whose absolute values depend on an assumption).
- The PC initialization of Fig. 10 is not given in the paper; this implementation selects
  by `Fobj` and is slightly better than the paper. To match the published readings exactly,
  revert to the residual criterion with an SVD-based initial value.

---

*Report generated 2026-11-14; the simulation code and this report were produced with
DeepSeek-V4.1-Flash, without access to the authors' original code. All numbers can be
regenerated with `matlab/run_all.m`; the random seed is fixed (`rng(1)`).*
