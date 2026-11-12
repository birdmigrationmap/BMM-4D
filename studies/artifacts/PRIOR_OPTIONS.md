# Prior and regularization terms: names and options (dev, 5 November 2026)

Research note; nothing in the toolbox or the manuscript changes. Starting point: the correlated prior of fix 3 ("GP 100 km", `prior_correlation_m`) removes the coverage edges, but it needs 5–8 times more PCG iterations and does not predict held-out radars better (README, "Smoothness term and held-out radars"). The questions:
- What does "prior" mean in the paper, and what does it mean in dev?
- Can the vertical profile, the smoothness term and the prior term become a single term?

**In short**
- **What "prior" means.** In the paper, "prior" is a *field*, the mean: $\mathbf E_{VP}\mathbf p_t$, or $\boldsymbol\rho_t^{prior}$ in the mosaic. The dev "GP prior" is not a new field. It keeps that mean and replaces the identity weight of the $\lambda_{VP}$ term with a correlated precision. So it changes the *covariance* of the deviation from the mean, not the mean.
- **The paper prior is white.** As a Gaussian prior, the paper's $\lambda_{VP}$ and smoothness terms together have a correlation length below one grid cell: 0.06–0.33 at 4 km, and 0 at 20 km. So a correction from the data stops where the beams stop. This is the source of the edges.
- **Three terms are one square.** The white prior, the dev GP and the smoothness term are the three terms of one square: $(\kappa^2-\Delta)^2=\kappa^4-2\kappa^2\Delta+\Delta^2$. This is the SPDE form of a Matérn field with $\nu=1/2$ in 3D, which is the smoothness the Appendix C fits find (median $\nu=0.46$). One term, $\lambda_b^2\lVert\mathbf K(\boldsymbol\rho_t-\boldsymbol\mu_t)\rVert^2$, with a horizontal and a vertical range, can replace four settings: $\lambda_{VP}$, $\lambda_{s,3D}$, the horizontal/vertical coefficient and the GP length.
- **The 7× cost comes from the preconditioner.** In a toy problem with the diagonal (Jacobi) preconditioner, PCG needs:
  - 14 iterations for the production prior;
  - 119 for the GP 100 km;
  - 256 to 1289 for one Matérn term at 20 to 50 km.

  With a preconditioner that keeps the prior block whole, all of them need 10–15.
- **Recommendation.** Option B: one Matérn prior covariance around the radar-profile field, tested first in the dev solver with a prior-aware preconditioner. If B does not predict held-out radars better, try option C, where the profile field is also unknown. Do the renaming of option A in any case.

## 1. The terms as the paper writes them

| Term | Paper (Methods II-B, Appendix D) | Role | Code (paper / v2) |
| --- | --- | --- | --- |
| Data misfit | $\lVert\mathbf y_t-\mathbf F_t\boldsymbol\rho_t\rVert_2^2$, divided by the observation count | likelihood | v2 lean: optional weights $1/\sigma_i^2$ from the estimator's sd, scaled to the mean variance |
| Vertical-profile regularization | $\lambda_{VP}^2\lVert\boldsymbol\rho_t-\mathbf E_{VP}\mathbf p_t\rVert^2$ for one radar; $\lambda_{VP}^2\lVert\boldsymbol\rho_t-\boldsymbol\rho_t^{prior}\rVert^2$ for the mosaic | pulls every cell toward the mean field, with the same weight everywhere | `prior` |
| Vertical profile (the mean) | $\mathbf p_t=\arg\min\lVert\mathbf y_t-\mathbf F_t\mathbf E_{VP}\mathbf p\rVert^2+\lambda_{s,1D}^2\lVert\mathbf L_{1D}\mathbf p\rVert^2$, fitted per radar in the mosaic and interpolated by IDW (3 radars, power 2.25): the "radar-density prior" | mean field | `lambda_profile`; v2 interpolates with Gaussian weights of 150 km |
| Smoothness | $\lambda_{s,3D}^2\lVert\mathbf L_{3D}\boldsymbol\rho_t\rVert^2$, a 7-point Laplacian of the *field*, with land–water links cut | "secondary role": grid-scale artefacts | `smoothness`; horizontal and vertical link weights $(1+c)/6$ and $(1-2c)/6$, with $c$ = `coeff_horver` |
| Advection | $\lambda_a^2\sum_t\sum_{j\to k}(\rho_{t+1,k}-\rho_{t,j})^2$ | links the time steps | `advection`; v2 relaxes links where the mean changes (exponent 2, floor 5 birds km⁻³) |

**Scaling.** The paper code and v2 scale the terms the same way. Each regularization block is scaled to the Frobenius norm of the data block: the $\lambda_{VP}$ block by $1/\sqrt n$, the smoothness block by $1/\lVert\mathbf L_{3D}\rVert_F$. So the $\lambda$ values are dimensionless strengths, not variances.

**Values.**

| Case | $\lambda_{VP}$ | $\lambda_{s,3D}$ |
| --- | ---: | ---: |
| Paper, static single radar | 0.01 | 0.03 |
| Paper, mosaic | 0.1 | 0.1 |
| Paper, Midwest | 0.1 | 0.3 |
| v2 default | 0.3 | 0.3 |

All use $\lambda_a=1$ and $c=0$, except Key West ($c=-0.5$).

**How the manuscript uses "prior":**
- "prior structure", in the objective;
- "horizontally uniform prior", for $\mathbf E_{VP}\mathbf p_t$;
- "prior field $\boldsymbol\rho_t^{prior}$ … the radar-density prior", in the mosaic;
- "initialized from the prior field", in the numerical implementation;
- "VP anchoring weight", "density-prior anchoring weight $\lambda_{VP}$" and "Regional prior and smoothness weights", in Appendix D.

So "prior" names the mean field, and $\lambda_{VP}$ is the weight that anchors the field to it.

**Appendix C** fits an anisotropic Matérn (k-Bessel) covariance to the transformed residual $Z=T_t\{y-p_t(z),z\}$. It is used only to simulate the synthetic truths; the inversion never uses it. The 50 fitted cases (read from `cov.mat` for this note):

| Parameter | q25 | Median | q75 | Note |
| --- | ---: | ---: | ---: | --- |
| $\nu$ | 0.29 | 0.46 | 0.64 | close to 1/2 (exponential) |
| $\delta_z$ | 144 m | 188 m | 240 m | |
| $\delta_r$ | 13.9 km | 19.1 km | 20 km | 24 of 50 cases at the 20 km upper bound; fitted on 10–40 km ranges only |
| $\sigma_Z^2$ / sample variance | | 0.68 | | the remaining third is a nugget (spatially uncorrelated) |

**Appendix D** covers:
- the numerical strategy: LSQR for the static case; for the joint systems, matrix-free PCG with a Jacobi preconditioner; projected FISTA for non-negativity; single-precision CG for the full-US run;
- the static grid of $\lambda_{VP}$ × $\lambda_{s,3D}$;
- the advection weight;
- the IDW calibration of the radar-density prior;
- the regional grid of $\lambda_{VP}$ × $\lambda_{s,3D}$ (minimax regret);
- the regional advection weight.

**Possible inconsistencies, for discussion (nothing corrected):**
1. The manuscript does not mention the horizontal/vertical coefficient $c$; it only says "built from horizontal and vertical cell connections". Key West used $c=-0.5$, which makes vertical links 4 times stronger than horizontal ones; Midwest used $c=0$.
2. The manuscript does not report the Key West weights, and the two code versions disagree. The paper script (S09_C) uses 0.1 / 0.1, $\lambda_a=1$, shoreline 0.3. The v2 `keywest` case, labelled "paper Key West published weights", uses 0.003 / 1 / 3. One of them did not produce Fig. 6.
3. The Midwest and Key West runs smooth the profiles with $\lambda_{s,1D}=10$ (S11_B, S09_C). The manuscript gives 0.1 (synthetic mosaic) and does not say that the real cases differ. In v2, a value of 10 was found to flatten the low-level peak (e895f31).
4. "Prior" names the mean field, the anchoring term and, in dev, a covariance.

## 2. Names

**What each "prior" is**
- **Manuscript "prior":** the mean field $\boldsymbol\mu_t$. It is $\mathbf E_{VP}\mathbf p_t$ for one radar and $\boldsymbol\rho_t^{prior}$ for the mosaic. The $\lambda_{VP}$ term anchors the field to this mean with an identity weight.
- **Dev "GP prior":** the same mean, but the identity is replaced by $\mathbf Q=c_Q(\mathbf I+(\ell/\Delta x)^2\mathbf L_h)$. This is a covariance of the deviation, not a prior field. $\ell=100$ km is its scale parameter, not its correlation range: the correlation is 0.35 at 20 km and 0.08 at 100 km (§3).

**Proposed names.** In the Bayesian reading, every quadratic penalty is a Gaussian prior:

| Object | Proposed name | Symbol | Replaces |
| --- | --- | --- | --- |
| Field the inversion is pulled toward | prior mean (the *profile field*) | $\boldsymbol\mu_t$ | vertical-profile prior, radar-density prior, $\boldsymbol\rho_t^{prior}$ |
| Deviation from it | anomaly | $\boldsymbol\rho_t-\boldsymbol\mu_t$ | — |
| Its weight matrix | prior precision, with covariance $\mathbf B=\mathbf Q^{-1}$: *white* (paper) or *Matérn* (dev) | $\mathbf Q$ | "GP prior", `prior_correlation_m` |
| The term | prior term, of strength $\lambda_b$ | $\lambda_b^2(\boldsymbol\rho_t-\boldsymbol\mu_t)^\top\mathbf Q(\boldsymbol\rho_t-\boldsymbol\mu_t)$ | the $\lambda_{VP}$ term, and the smoothness term (§3) |
| Correlation scales | horizontal and vertical range | $\ell_r$, $\ell_z$ | GP length, `coeff_horver` |
| Advection term | transport (model-error) term | $\lambda_a$ | unchanged |
| Data term | likelihood, with observation errors | $\mathbf R$ | misfit |

Drop "GP prior". Every option in this note is a Gaussian-process prior; what differs between them is the covariance. Suggested v2 option names (code only, not maths): `prior_covariance` ("white" or "matern"), `prior_range_h_m`, `prior_range_v_m`.

## 3. What each term is, and where the terms overlap

**A penalty is a Gaussian prior.** Minimizing $\mathcal J$ gives the maximum a posteriori field of the model
- $\mathbf y_t\mid\boldsymbol\rho_t\sim\mathcal N(\mathbf F_t\boldsymbol\rho_t,\mathbf R)$,
- $\boldsymbol\rho\sim\mathcal N(\mathbf m,\mathbf Q_{tot}^{-1})$,

where $\mathbf Q_{tot}$ collects the regularization blocks. Each $\lambda_k^2$ is the ratio of an observation-error variance to a prior variance (Kaipio & Somersalo 2005).

**The paper's static prior.** Its prior part is $\lambda_{VP}^2\lVert\boldsymbol\rho-\boldsymbol\mu\rVert^2+\lambda_{s,3D}^2\lVert\mathbf L_{3D}\boldsymbol\rho\rVert^2$. This is a Gaussian Markov random field with:
- precision $\mathbf Q_{paper}=\lambda_{VP}^2\mathbf I+\lambda_{s,3D}^2\mathbf L_{3D}^\top\mathbf L_{3D}$;
- mean $\lambda_{VP}^2\mathbf Q_{paper}^{-1}\boldsymbol\mu$.

Because the smoothness acts on $\boldsymbol\rho$ and not on $\boldsymbol\rho-\boldsymbol\mu$, it also penalizes the curvature of the profile itself (only mildly at the paper weights).

In Fourier space the precision is $\lambda_{VP}^2+\lambda_{s,3D}^2\,l(\mathbf k)^2$, with $l(\mathbf k)\propto\lvert\mathbf k\rvert^2$ at low wavenumbers: a white "nugget" plus a biharmonic term. With the solver's scaling, its correlation length is below one cell for every paper and v2 setting.

**Implied correlations** (`ar_prior_spectrum`; periodic lattice, 4 km × 100 m; on the 2 km paper grid, halve the distances):

| Prior | r at 4 km | r at 20 km | r at 52 km | r at 100 km | Vertical 1/e | Prior-block condition number |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| paper static 0.01 / 0.03; Midwest 0.1 / 0.3 | 0.33 | 0.00 | 0.00 | 0.00 | 95 m | 32 |
| paper mosaic 0.1 / 0.1; v2 0.3 / 0.3 | 0.14 | 0.00 | 0.00 | 0.00 | 73 m | 4.4 |
| paper Key West 0.1 / 0.1, $c=-0.5$ | 0.06 | 0.00 | 0.00 | 0.00 | 85 m | 4.2 |
| dev GP 100 km + smoothness 0.3 | 0.57 | 0.18 | 0.07 | 0.03 | 86 m | 7.7·10³ |
| dev GP 100 km, no smoothness | 0.68 | 0.35 | 0.18 | 0.09 | 63 m | 5.0·10³ |
| one Matérn term ($\alpha=2$), 20 km × 190 m | 0.86 | 0.38 | 0.07 | 0.01 | 188 m | 4.6·10⁴ |
| one Matérn term ($\alpha=2$), 50 km × 190 m | 0.96 | 0.73 | 0.38 | 0.14 | 192 m | 1.6·10⁶ |

What the table shows:
- **Paper and v2 production.** The prior is white at the reconstruction scale, both horizontally and vertically (less than one layer). A correction made by the data stops where the beams stop, so the discs and arcs of fix 3 follow from the prior itself, not from an error.
- **Dev GP.** Its correlation is peaked at short lags and has a heavy tail. It is horizontal only, so it is vertically white (63 m) and still needs the smoothness term in the vertical: without it, AR_H found profiles 3.6 times rougher. In 2D, $\mathbf I+\ell^2\mathbf L_h$ is the $\alpha=1$ SPDE, a Matérn field with $\nu=0$, whose pointwise variance grows logarithmically as the grid is refined. This is why 93683f7 needs a numerical normalization.
- **One Matérn term.** It reproduces the form of Appendix C exactly: an exponential, with correlation 1/e at $\ell_r$ and $\ell_z$.

**The SPDE link.** A Matérn field with smoothness $\nu$ in $d$ dimensions solves $(\kappa^2-\Delta)^{\alpha/2}x=\mathcal W$, with $\alpha=\nu+d/2$ (Whittle 1954; Lindgren, Rue & Lindström 2011). On a grid, its precision is the sparse matrix $\tau^2(\kappa^2\mathbf I+\mathbf L)^\alpha$ for integer $\alpha$. In 3D with $\nu=1/2$ (Appendix C median: 0.46), $\alpha=2$, and

$$
(\kappa^2\mathbf I+\mathbf L)^2=\kappa^4\mathbf I+2\kappa^2\mathbf L+\mathbf L^2 .
$$

These are the three dev terms:
- the white prior, $\kappa^4\mathbf I$;
- the dev GP, $2\kappa^2\mathbf L$ (the $\alpha=1$ term);
- the smoothness, $\mathbf L^2=\mathbf L^\top\mathbf L$.

With anisotropy, written in stretched coordinates, the single term is

$$
\lambda_b^2\,\lVert\mathbf K(\boldsymbol\rho_t-\boldsymbol\mu_t)\rVert_2^2,\qquad
\mathbf K=\mathbf I+\frac{\ell_r^2}{\Delta x^2}\mathbf L_h+\frac{\ell_z^2}{\Delta z^2}\mathbf L_v,
$$

where $\mathbf L_h$ and $\mathbf L_v$ are the horizontal and vertical graph Laplacians of the active cells, with unit link weights and the land–water links cut. Its covariance is $\sigma^2\exp(-d)$, with $d=\sqrt{(\Delta r/\ell_r)^2+(\Delta z/\ell_z)^2}$. This is the Appendix C model with $\nu=1/2$, $\delta_r=\ell_r$ and $\delta_z=\ell_z$.

The current terms differ from this square in three ways:
1. the smoothness acts on $\boldsymbol\rho$, not on the anomaly;
2. the three weights are set independently, so their sum is not a Matérn covariance;
3. the dev GP and the smoothness have different horizontal/vertical structures (the GP is horizontal only; the smoothness depends on $c$).

**Other overlaps**
- **$c$ is a ratio of ranges.** It sets vertical against horizontal smoothing; in $\mathbf K$ it becomes $\ell_z/\ell_r$, in physical units.
- **Vertical correlation is set twice.** The profile fit penalizes first differences in height (a random walk in height, an intrinsic GMRF), and $\ell_z$ again encodes vertical correlation. That is two parameters for one property, at two stages.
- **The profile field is a plug-in mean.** It is estimated from the same data in a separate step and then treated as known (empirical Bayes). So the data enter twice, and leaving a radar out requires rebuilding the mean (the held-out leak of fix 2b).
- **The interpolation is fixed-weight kriging.** The IDW (paper) or Gaussian-weight (v2) interpolation of the radar profiles is a fixed-weight smoother of a coarse field.
- **Advection is a separate process term.** It is a Gaussian on the transport residual $\boldsymbol\rho_{t+1}-\mathbf A_t\boldsymbol\rho_t$: the model-error term of weak-constraint 4D-Var (Trémolet 2006). It couples the time steps, so the space–time prior precision is $\mathrm{blockdiag}_t(\lambda_b^2\mathbf Q)+\lambda_a^2\mathbf D_A^\top\mathbf D_A$. It describes a process, not a spatial covariance, and it stays a separate term in every option below.

## 4. Options

### A. Rename and document (method unchanged)
- **What:**
  - use the §2 names;
  - give the Bayesian reading of each term in Methods;
  - say that the smoothness and $\lambda_{VP}$ terms together form a short-range Gaussian prior;
  - report $c$, the Key West weights and the $\lambda_{s,1D}$ of the real cases (flags of §1).
- **Pros:** no rerun; removes the ambiguity of "prior"; the paper stays v1.
- **Cons:** keeps four unitless strengths; does not fix the edges.
- **Cost and risk:** none. Wording that adds claims still needs approval.

### B. One Matérn prior covariance around the profile field

$$
\mathcal J=\sum_t\lVert\mathbf R_t^{-1/2}(\mathbf y_t-\mathbf F_t\boldsymbol\rho_t)\rVert^2
+\lambda_b^2\sum_t\lVert\mathbf K(\boldsymbol\rho_t-\boldsymbol\mu_t)\rVert^2
+\lambda_a^2\sum_t\sum_{j\to k}(\rho_{t+1,k}-\rho_{t,j})^2
$$

- **Replaces:** the $\lambda_{VP}$ term, the smoothness term, $c$ and the GP length. Four settings become three parameters: $\lambda_b$, $\ell_r$ and $\ell_z$.
- **$\alpha=2$ is not a free choice:** it follows from $\nu\approx1/2$ in 3D.
- **Values:**
  - $\ell_z\approx190$ m, from Appendix C (interquartile range 144–240 m).
  - $\ell_r$ is at least about 20 km in Appendix C, but that fit hit its bound and used 10–40 km ranges only. Fix 3 favoured 50–150 km for the edges. Set $\ell_r$ from cross-radar statistics (§7, stage 0) and from held-out radars.
  - $\lambda_b$ keeps the current marginal-variance normalization, so its value has the same meaning as $\lambda_{VP}$.
- **Pros:**
  - One term with one sparse square-root operator, the same structure as the current smoothness block.
  - Parameters with units and a data source. Appendix C then links the simulator and the inversion.
  - It acts on the anomaly, so the profile shape is not flattened.
  - The cut coastline works as a barrier model (Bakka et al. 2019).
  - The variance is finite in 3D, so no grid-dependent normalization is needed.
  - It is a proper Gaussian prior: with the estimator's $\mathbf R$, the same terms give the posterior uncertainty.
- **Cons:**
  - With Jacobi, PCG is much slower. In the toy: 256 iterations at 20 km and 1289 at 50 km, against 14 for production and 119 for the dev GP. B therefore needs a prior-aware preconditioner (§5); this is the main engineering risk.
  - The mean stays two-stage, with its double use of the data and the held-out leak.
  - The dev GP did not predict held-out radars better, and B changes the same thing (the covariance). Expect cleaner edges and fewer parameters, not a better held-out score.
  - Inverse crime in the synthetic benchmarks if $\ell_r$ and $\ell_z$ come from the same per-case fits that simulate the truths.
- **Cost:** per iteration, about the same as production (two products with $\mathbf K$ per time step). The total cost depends on the preconditioner.
- **Variant B′, advect the anomaly:** use $\lambda_a^2\sum(\delta_{t+1,k}-\delta_{t,j})^2$ with $\boldsymbol\delta_t=\boldsymbol\rho_t-\boldsymbol\mu_t$.
  - The mean then carries its own change (take-off after dusk, landing before dawn), which would replace the relaxation exponent and floor (two fixed numbers).
  - It changes the physical meaning of the advection term, since the mean is not transported. Not for the first test.

### C. One hierarchical model: the profile field is unknown too

$$
\boldsymbol\rho_t=\mathbf W\mathbf c_t+\boldsymbol\delta_t
$$

- **The parts:**
  - $\mathbf c_t$ is a coarse field: either the radar profiles, with $\mathbf W$ the interpolation weights, or the nodes of a coarse grid (for example every 8 cells), with $\mathbf W$ bilinear. Its prior $\lVert\mathbf K_c\mathbf c_t\rVert^2$ has a long range (the 150 km now in the Gaussian weights), plus the vertical and time random walks of the profile fit.
  - $\boldsymbol\delta_t$ is the fine anomaly: white, like the paper's $\lambda_{VP}$ term, or a short-range Matérn as in B.
  - Both are estimated in one solve with the data and advection terms.
- **Replaces:** the separate profile fit, the IDW or Gaussian weights and their calibration (Appendix D), $\lambda_{VP}$, the smoothness and the GP. The paper's two-stage prior is the plug-in special case of C: fit $\mathbf c$ first, then hold it fixed. This is a multiresolution model (LatticeKrig; Nychka et al. 2015).
- **Pros:**
  - Everything is in one model, and the data are used once.
  - A held-out test only needs zero weights: no rebuilt mean, so no leak by construction.
  - The kriging weights adapt to the radar layout.
  - Advection also carries the coarse level.
  - Cheap to precondition: the coarse $\mathbf K_c$ is small enough to factorize exactly, and the white fine part works with Jacobi.
  - The option most likely to change held-out prediction, because the held-out error is set by the mean between radars.
- **Cons:**
  - The largest change, in both the method and the paper.
  - Near the radars, how the coarse and fine parts split depends on their variances.
  - A per-radar calibration offset still ends up in the field. The sandbox's "per-radar level" belongs in the observation model; no option here addresses it.
  - New solver code (augmented unknowns).
- **Cost:** per iteration, about production plus the coarse products. Iterations close to production with the block preconditioner (expected, not tested).

### D. Variance on relative or log density (a modifier of B or C)
- **D1, relative anomaly.** Weight the anomaly by $\mathbf S_t^{-1}$, with $\mathbf S_t=\mathrm{diag}(\boldsymbol\mu_t+\sigma_0)$: $\lambda_b^2\lVert\mathbf K\mathbf S_t^{-1}(\boldsymbol\rho_t-\boldsymbol\mu_t)\rVert^2$.
  - Still quadratic, so still one PCG solve.
  - Anomalies scale with the local level (multiplicative errors), as the level-dependent spread of Appendix C suggests and as the sandbox concluded (2 November).
  - $\sigma_0$ is the observation-error floor that v2 already has, so no new constant.
  - Cheap.
- **D2, log density.** $\boldsymbol\rho_t=\exp(\mathbf u_t)$, with prior $\lVert\mathbf K(\mathbf u_t-\log\boldsymbol\mu_t)\rVert^2$.
  - Non-negative by construction (no projected FISTA), and closest to the transformed space of Appendix C.
  - Non-linear: Gauss–Newton or L-BFGS, with 5–10 PCG solves.
  - Problems with very low densities aloft, and a bias in the mean ($e^{\sigma^2/2}$).
  - The highest cost and risk.

### Summary

| | Prior terms | Parameters | Removes | Toy PCG iterations (Jacobi / prior-aware) | Code | Risk | Manuscript |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Paper / v2 production | $\lambda_{VP}$, smoothness | $\lambda_{VP}$, $\lambda_{s,3D}$, $c$ (+ $\lambda_{s,1D}$, interpolation) | — | 14 / 15 | — | — | — |
| Dev fix 3 | $\lambda_{VP}$ with GP precision, smoothness | + $\ell$ | — | 119 / 12 | exists | 7× cost; no held-out gain | — |
| A | unchanged | unchanged | the ambiguity | 14 | none | none | wording |
| B | one Matérn term | $\lambda_b$, $\ell_r$, $\ell_z$ | smoothness, $c$, GP | 256–1289 / 10 | solver option and preconditioner | the preconditioner | Methods equation; Appendix D rerun |
| C | coarse and fine field, one solve | coarse range and sd, fine sd, $\ell_z$ | profile step, interpolation, smoothness, $c$, GP | not tested | augmented solver | the largest | Methods rewrite; Appendices C and D |
| D1 / D2 | + relative / log | + $\sigma_0$ (exists) | — / FISTA | as B / 5–10 solves | small / large | low / high | a sentence / a section |

### Recommendation
Test **B, with D1 as a variant**, in the dev solver first.
- It is the elegant form asked for: the three terms are one square, and the two ranges have units and a data source.
- It keeps the paper's structure (a mean built from radar profiles), so the manuscript change is one equation and a recalibration.

Run it only if the preconditioner works (§5), and judge it on held-out radars (§7). If B removes the edges but does not predict held-out radars better than the white prior, the lever is the mean, not the covariance: go to **C**. Do **A** in the manuscript in any case; the names and the flags of §1 do not depend on the method.

## 5. Why the GP costs 7×, and the fix
- **The preconditioner cannot see the prior's couplings.** The joint system is solved by PCG with a diagonal (Jacobi) preconditioner, which rescales each cell but cannot undo couplings between cells. The slowest modes are long-wavelength anomalies in cells that no beam reaches, which only the prior holds.
- **Their precision is far below the diagonal.** It is $\lambda_b^2c_Q$, against a diagonal of $\lambda_b^2c_Q(1+4\ell^2/\Delta x^2)$. The prior block's condition number is $1+8\ell^2/\Delta x^2$ for the dev GP (5·10³ at 100 km on 4 km cells), and its square for $\alpha=2$. CG iterations grow with its square root.

**Toy problem** (`ar_prior_pcg_toy`): 160 × 160 cells at 4 km; four radars 200 km apart, each seeing half the cells within 120 km through a 3 × 3 footprint; the solver's scaling; tolerance 10⁻⁴.

| Prior | Jacobi | Prior-aware |
| --- | ---: | ---: |
| white + smoothness 0.3 (v2 production) | 14 | 15 |
| dev GP 100 km + smoothness 0.3 | 119 | 12 |
| one Matérn term, 20 km | 256 | 10 |
| one Matérn term, 50 km | 1289 | 10 |

The toy reproduces the ratio seen on ar_plains: 119 / 14 = 8.5 in the toy, 281 / 38 = 7.4 there.

**The prior-aware preconditioner** keeps the prior block whole: per time step, $\mathbf M_t=\lambda_b^2\mathbf K^\top\mathbf K+\mathrm{diag}(\text{data}+\text{advection})$. The toy factorizes it exactly. At the size of us_east (millions of cells per time step) it needs an approximate inverse. Candidates:
1. **Two solves with $\mathbf K$ per application:** $\mathbf M_t^{-1}\approx\mathbf K^{-1}\mathbf D_t^{-1}\mathbf K^{-\top}$. $\mathbf K$ is a sparse 7-point M-matrix, so multigrid works, or an inner PCG with incomplete Cholesky.
2. **The same, with $\mathbf K^{-1}$ by a discrete cosine transform on the bounding box.** It is exact on a rectangular grid and approximate with the active mask and the coastline cuts, which is enough for a preconditioner.
3. **The control-variable transform of variational assimilation:** $\boldsymbol\rho_t=\boldsymbol\mu_t+\mathbf K^{-1}\mathbf v_t$. The prior becomes $\lambda_b^2\lVert\mathbf v_t\rVert^2$ and the Hessian is close to the identity. Weaver & Courtier (2001) build correlation operators the same way, from a diffusion equation.
4. **Option C's coarse field:** it is small enough to factorize exactly.

Advection couples the time steps. All four are per-time preconditioners, as the current Jacobi one is.

## 6. What option B would change (proposal only)

**Manuscript.** This needs explicit approval, because it changes the mathematical meaning of the regularization.
- **Methods, static reconstruction:**
  - replace the $\lambda_{VP}$ and $\lambda_{s,3D}$ terms by the prior term $\lambda_b^2\lVert\mathbf K(\boldsymbol\rho_t-\boldsymbol\mu_t)\rVert^2$, with $\boldsymbol\mu_t=\mathbf E_{VP}\mathbf p_t$;
  - replace the smoothness paragraph by one paragraph on the SPDE link and the ranges ($\alpha=2$ from Appendix C), and the coastline as a barrier.
- **Multi-radar reconstruction:** the same term, with $\boldsymbol\mu_t=\boldsymbol\rho_t^{prior}$.
- **Numerical implementation and Appendix D:**
  - describe the preconditioner;
  - recalibrate $\lambda_b$ × $\ell_r$ in place of $\lambda_{VP}$ × $\lambda_{s,3D}$, and rerun the static and regional sensitivity figures;
  - take $\ell_z$ from Appendix C.
- **Appendix C:** add one sentence saying that the inversion's prior covariance has the simulator's form, and how the inversion's ranges were set (avoiding the inverse crime).
- **The §1 flags,** whichever option is chosen.
- **Timing:** every number in the paper changes, so this belongs to the v2 rerun of the manuscript, not to the v1 text.

**v2 code.** Dev copy first; the toolbox only after approval.
- **solve_advection_lean:**
  - apply the prior term as $\mathbf K(\boldsymbol\rho_t-\boldsymbol\mu_t)$ (square-root form, on the active cells, on the terrain-following graph when the AGL or sigma links are used);
  - remove the smoothness term and `coeff_horver` from the lean path;
  - replace `prior_correlation_m`, `prior_correlation_vertical_m` and `smoothness` with `prior_covariance`, `prior_range_h_m` and `prior_range_v_m`;
  - add a preconditioner option.
- **Other files:**
  - the inversion struct in `config.m`, and the pass-through in `invert.m`;
  - a MIGRATION.md entry;
  - the "paper" solver path stays unchanged (white prior + smoothness, checked by the paper regression test);
  - `build_density_prior_fast` stays unchanged (B keeps the mean).

## 7. Later compute plan (not run)

**Stage 0: ranges from the data (minutes, no inversion)**
- **$\ell_z$ and $\nu$:** from Appendix C (done above).
- **$\ell_r$:** use the radar profiles (`info.radar_profiles` of `build_density_prior_fast`) of the shared night set. That set is 17 Dark Ecology nights from 2021–2025, both seasons, plus the 21–22 October 2026 reference: 107 us_east radars, every 30 min, cleaned with production settings and `"shape", false`. The list is in `studies/nights/figures/NS_A_nights.csv` (the rows whose uses include "prior"), and the plan in `studies/nights/PLAN.md`. Start with ar_plains_bb and ar_midwest, which are cached on the laptop.
  - take the anomaly of each radar profile from the network mean at the same height (as a log ratio);
  - correlate it between radar pairs against their distance (pairs 100–800 km apart), and fit an exponential per night and per season;
  - add the within-radar variogram of observation anomalies out to 120 km, which extends Appendix C beyond 40 km.

  This sets $\ell_r$ without using the simulator, so there is no inverse crime. Eighteen nights also show whether one range holds across nights and seasons. If the observations have to be built for every night, this is a cluster job rather than a matter of minutes (`hpc/README.md`).

**Stage 1 (about 1 h; dev solver; cached cases; production settings, with AGL links and prior as in the full-night run)**
- **Cases:**
  - 1a. ar_plains_bb (12 radars, 7 frames; tests the edges): every variant. Held out: KLNX, KUEX, KGLD, KOAX, KTWX (as AR_H).
  - 1b. ar_appalachians_bb (16 radars; close radars, terrain): V0, V1 and the best B. Held out: KMRX, KGSP, KFCX, KRLX, KJKL, KHTX (as AR_G).
- **Variants:**
  - V0: white + smoothness 0.3 (production);
  - V1: dev GP 100 km + smoothness 0.3 (the current fix);
  - V2: option B with $\ell_r$ = 25, 50 or 100 km, $\ell_z$ = 190 m, the same $\lambda_b$ normalization and the prior-aware preconditioner; plus one Jacobi run at 25 km, for the cost;
  - V3: the best V2 plus D1.
- **Held-out runs:** the prior is rebuilt without the held-out radar (`ar_holdout_score`); everything else as in AR_H.
- **Budget:** on ar_plains, one solve took 22 s in production and 193 s with the GP. Case 1a is 6 variants × 6 solves = 36 solves; case 1b is 3 × 7 = 21. At 60 s or less per solve with a working preconditioner, that is about 1 h. If the preconditioner does not work, run only 1a with Jacobi at 25 km, then stop.

**Metrics** (per variant and case):
- held-out relative L2 error and median absolute log ratio at 30–120 km: the mean over held-out radars, and per radar;
- edge score (`ar_edge_score`);
- misfit, on all rows and on rows at 30–120 km;
- roughness along z and along x;
- PCG iterations, wall time and peak memory;
- negative share before FISTA.

**Go / no-go for stage 2** (fixed before the runs; comparisons only, no new thresholds):
- **GO** if one B variant, on both cases:
  1. has held-out relative error and median absolute log ratio no worse than V0 on the mean, and better or equal at half or more of the held-out radars;
  2. has an edge score no worse than V1;
  3. takes no more wall time than V1.
- **If 1 and 2 hold but not 3:** work on the preconditioner before any larger run.
- **If 2 holds but 1 fails for every B variant:** keep V1 as the dev option, do the renaming (A), and test C, because the mean is the lever.
- **If 2 fails:** B's ranges are too short; stop B.

**Stage 2 (if GO; a few hours)**
- ar_midwest (28 radars), held out KDMX, KDVN, KOAX, KEAX, KARX (as AR_I).
- Full us_east nights with the winning variant and the current fix, on the cluster (`hpc/submit.sh`, about 8 threads and 36 GB per solve; the GP run took 64–100 min and 352 iterations). Use the night-set nights that are dry and strong over the Plains radars, where the edges show: 2021-05-12, 2023-05-22, 2024-04-14 and 2024-05-10. Run at most about 4 large jobs at once on the shared node, and ask the user before submitting. Compare the edge score, the held-out radars and the maps with the current fix (the video on the Mac: Linux has no MP4 export).
- The regional synthetic of Appendix D (5 cases): NRMSE against the truth, with $\ell_r$ from stage 0, not from the simulator fits.

## References
- Bakka, H., Vanhatalo, J., Illian, J., Simpson, D., Rue, H. (2019). Non-stationary Gaussian models with physical barriers. *Spatial Statistics* 29, 268–288.
- Kaipio, J., Somersalo, E. (2005). *Statistical and Computational Inverse Problems*. Springer.
- Lindgren, F., Rue, H., Lindström, J. (2011). An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. *J. R. Stat. Soc. B* 73, 423–498.
- Nychka, D., Bandyopadhyay, S., Hammerling, D., Lindgren, F., Sain, S. (2015). A multiresolution Gaussian process model for the analysis of large spatial datasets. *J. Comput. Graph. Stat.* 24, 579–599.
- Rue, H., Held, L. (2005). *Gaussian Markov Random Fields: Theory and Applications*. Chapman & Hall/CRC.
- Trémolet, Y. (2006). Accounting for an imperfect model in 4D-Var. *Q. J. R. Meteorol. Soc.* 132, 2483–2504.
- Weaver, A., Courtier, P. (2001). Correlation modelling on the sphere using a generalized diffusion equation. *Q. J. R. Meteorol. Soc.* 127, 1815–1846.
- Whittle, P. (1954). On stationary processes in the plane. *Biometrika* 41, 434–449.

## Files
- `functions/ar_prior_spectrum.m`: the correlation table of §3. Takes seconds and needs no data.
- `functions/ar_prior_pcg_toy.m`: the PCG table of §5. Takes seconds and needs no data.
- The Appendix C fits are read from BMM-within `data/derived/simulator/cov.mat` (variable `cov.parm`: sill, vertical range, horizontal range, $\nu$).
