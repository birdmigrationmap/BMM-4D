# Near radar: dropping the steep cuts close to each radar

Status 12 Nov 2026. Method and decision rule: `PLAN.md`. The production rule is drafted on `v3-dev` (`bmm4d.forward.drop_near_radar_cuts`, `c.near_radar`) and is a v3 candidate.

## Validation set and runs

The validation set is the 18 nights of the shared night set (`studies/nights`): 17 Dark Ecology nights from 2021–2025 in six migration phases, plus 21–22 Oct 2026. The 107 us_east radars are scanned every 30 min, cleaned with the production settings without the shape correction.

The runs:
- **Diagnostic:** every volume, 1,812 radar-nights.
- **Gap filling:** 36 regional solves (Midwest and East, each night), one per rule.

| Job | What | Time | Peak memory |
|---|---|---:|---:|
| 58740 | `nr_diagnose_nights`, all 18 nights | 30 min | 9.8 GB |
| 58741 | `nr_gap_fill` Midwest, 18 nights × 4 rules | 12 h 09 min | 56 GB |
| 58742 | `nr_gap_fill` East, 18 nights × 4 rules | 4 h 11 min | not read |
| 58867 | `NR_E_multi_night.m`: collect, analyse, plot | 4 min | |

The products are in `products/near_radar/NR_E/` on the unit share. The tables and the figure are in `figures/NR_E_*`.

Run 1 (NR_B–NR_D locally, including the synthetic gap-fill test and the NR_C clutter flag) has not been run. Steps 2 and 5c of the decision rule are therefore open.

## Diagnostic: the deficit by cut type

Each value is observed density divided by the density predicted by a beam-weighted profile fitted to the split cuts at 20–60 km, pooled over all volumes. Source: `NR_E_by_type.csv`.

| Range (km) | Split | Batch | Staggered | Doppler |
|---|---:|---:|---:|---:|
| 3–5 | 0.70 | 0.51 | 0.41 | 0.79 |
| 5–7.5 | 0.83 | 0.63 | 0.57 | 0.98 |
| 7.5–10 | 0.91 | 0.71 | 0.73 | 1.06 |
| 10–12.5 | 0.99 | 0.77 | 0.82 | 1.10 |
| 12.5–15 | 1.02 | 0.80 | 0.90 | 1.07 |
| 15–20 | 1.04 | 0.82 | 0.93 | 1.07 |
| 20–25 | 1.04 | 0.85 | 0.97 | 1.09 |
| 25–30 | 1.02 | 0.86 | 0.98 | 1.05 |
| 30–40 | 1.01 | 0.87 | 0.98 | 0.99 |
| 40–50 | 0.98 | 0.89 | 1.00 | 0.92 |

What the table shows:
- The batch cuts read low at every range. They also carry a far offset of about 0.87–0.89 against the split cuts at 30–50 km.
- Relative to each radar-night's own batch level at 25–50 km (`NR_E_deficit_spread.csv`), the median deficit is:
  - 0.42 at 3–5 km;
  - 0.66 at 7.5–10 km;
  - 0.85 at 12.5–15 km;
  - 0.90 at 15–20 km;
  - 0.95 at 20–25 km;
  - 0.98 at 25–30 km;
  - 0.99 at 30–40 km.

  The spread over radar-nights is wide: at 10–12.5 km the 10–90% range is 0.40–1.22.
- The split cuts read low only within ~7.5 km. They stay at every range in any rule.
- **Adaptive R_min per radar-night:** median 15 km, 10–90% 0–25 km. It exceeds 20 km in 31% of radar-nights. It does not follow track speed (Spearman 0.05) or density (−0.07). By season, phase and VCP it varies little (medians 12.5–15 km, and 20 km for VCP 32).
- **Radar or night** (two-way model of the log deficit at 5–15 km): radar explains a sum of squares of 73.2 (df 103), night 24.6 (df 17), with a residual of 302.4. Night has the larger effect per degree of freedom (F = 8.1 vs 4.0), but radar explains three times more of the variance. The deficit is not mostly a night effect.

## Gap filling: four rules on 36 night-regions

All four rules drop near-radar entries of the batch cuts, then rebuild the prior. The split rows within 20 km of every second radar are held out as an independent check.

| Rule | Column 0–10 km / 60–100 km | Column 10–20 km | Held-out split obs/pred | Entries dropped | Misfit 30–120 km |
|---|---|---|---|---|---|
| All rows | 0.77 (0.70–0.89) | 0.88 (0.81–0.99) | 1.14 (0.99–1.40) | 0 | 0.318 |
| Batch < 20 km (fixed) | 0.95 (0.87–1.01) | 0.97 (0.91–1.03) | 1.07 (0.93–1.34) | 2.9% | 0.322 |
| Batch < adaptive R_min | 0.93 (0.87–0.99) | 0.95 (0.90–1.00) | 1.09 (0.93–1.34) | 2.3% | 0.319 |
| Whole batch cuts | 0.99 (0.94–1.06) | 1.02 (0.97–1.08) | 1.05 (0.93–1.31) | 22% | 0.306* |

How to read the table:
- Values are medians over the 36 night-regions (10–90% in brackets). Source: `NR_E_decision.csv`; the per night-region values are in `NR_E_summary.csv`.
- \* Whole batch cuts removes batch entries at 30–120 km as well, so its misfit is over a different set of entries and is not comparable.
- Paired differences between rules, mean ± se over the 36 night-regions:

  | Pair | Column 0–10 km | Column 10–20 km | Held-out split | Misfit 30–120 km |
  |---|---|---|---|---|
  | Fixed − all rows | +0.164 ± 0.010 | +0.077 ± 0.007 | −0.059 ± 0.006 | +1.7% median (10–90%: 0.4–3.8%) |
  | Whole cut set − fixed | +0.045 ± 0.004 | +0.052 ± 0.005 | | |
  | Adaptive − fixed | −0.019 ± 0.004 | | | |

- **By region**, median column at 0–10 km for all rows / fixed / whole cut set:
  - Midwest: 0.82 / 0.97 / 1.01;
  - East: 0.75 / 0.92 / 0.97.

## Decision rule applied (PLAN.md)

1. **Which cuts.**
   - **Batch:** dropped.
   - **Staggered:** provisionally dropped too. Pooled, they read 0.57–0.93 at 5–15 km, close to batch. In 5 of the 7 radar-nights with staggered cuts, the mean deficit at 5–15 km is below 0.95. All 7 come from one night (21 Oct 2026), and the per radar-night se was not written out, so the 2-se condition is not checked. To confirm on more recent nights.
   - **Doppler:** stays. It reads 0.98–1.10 at 5–30 km.
2. **Clutter flag:** not assessed (NR_C has not been run).
3. **Fixed R\*: not settled.** Only 20 km and the whole cut were solved. Going from 20 km to the whole cut still raises the column at 0–20 km by about 10 se (0.045–0.052, se 0.004–0.005), so by the rule R\* is above 20 km. The diagnostic points to 30 km: the median deficit is 0.95 at 20–25 km and 0.98–0.99 beyond 25 km. A gap-fill run with batch < 30 km (and 25 km) is needed to fix R\*.
4. **Fixed, not adaptive.** Adaptive fails the first condition. Its worst 10% column at 0–10 km (mean of the 4 lowest night-regions) is 0.854, against 0.864 for fixed; its median is also lower (0.93 vs 0.95). It does drop fewer entries (2.3% vs 2.9%). The deficit varies more by radar than by night, so a static per-radar range table would be the alternative to adaptive. It is not needed while fixed does better.
5. **Gap filling.**
   - (a) **Misfit at 30–120 km:** fails narrowly with batch < 20 km. Its median relative change is +1.7% (10–90%: 0.4–3.8%) on the same entries, against the 1% tolerance.
   - (b) **Held-out split ratio:** passes. It moves toward 1 (1.14 → 1.07, and 1.05 for the whole cut).
   - (c) **Synthetic test:** not run.

## Outcome

- **The rule works:** dropping the batch cuts near the radar lifts the near-radar column from 0.77 to 0.95 of the far level (0–10 km), and the gap is filled from the split cuts, the neighbouring radars and advection. Dropping the whole batch cuts reaches 0.99, but removes 22% of all entries.
- **Candidate for v3:** drop batch (and probably staggered) cuts within a fixed R\*, keep the split and Doppler cuts at every range, and keep all biological echo in the cache.
- **R\* is not fixed yet:** 20 km is not enough by the decision rule. 30 km is the likely value.

Open before production:
- a gap-fill run at R\* = 25 and 30 km, Midwest and East. One rule per night-region instead of four, so roughly 3 h Midwest and 1–2 h East on the cluster.
- whether the misfit increase (+1.7% at 20 km) shrinks at 30 km;
- run 1: the synthetic gap-fill test and the clutter flag;
- staggered cuts on more nights with staggered VCPs;
- the far offset of batch against split (about 0.87 at 30–50 km). It is outside this rule's range, but it biases the upper layers.
