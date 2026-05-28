# GLM and data processing — worked example

A fully worked, runnable R example of the **count-data GLM workflow** used in
ecological video / survey studies: read raw effort + sightings → clean →
compute MaxN per deployment → join keeping true zeros → choose a Poisson or
negative-binomial GLM → test against the null model with a likelihood-ratio
χ² test → report rate ratios. Heavily commented for teaching.

The data are **fully synthetic** (Mesoamerican Barrier Reef video surveys of
Nassau grouper). The species and study area are illustrative — the workflow
is the point and is directly transferable to your own count data.

## How to use

```bash
git clone git@github.com:sharklover6000/GLM-and-data-processing-example-script-.git
cd GLM-and-data-processing-example-script-
Rscript example_GLM_workflow.R
```

Requires R with `dplyr` and `MASS` (`install.packages(c("dplyr","MASS"))`).
Output (a figure) is written to `outputs/`.

## Files

| File | What it is |
|---|---|
| `example_GLM_workflow.R` | The main teaching script — clean → join → model → report |
| `R/simulate_data.R` | How the example `data/*.csv` files were generated (for transparency / reproducibility) |
| `data/effort.csv` | One row per video deployment (some are failed and must be excluded) |
| `data/sightings.csv` | One row per observation; many species |
| `data/species_lookup.csv` | Species code → scientific + common name |
| `outputs/` | Figures and captured run output |

## The workflow this script teaches

1. **Read** the raw files (`read.csv`) and inspect them (`str`, `head`, `table`).
2. **Clean the effort table** — drop **failed deployments** (no usable video):
   these are *missing data*, not zeros.
3. **Set factor reference levels** so model coefficients read naturally.
4. **Compute the focal-species MaxN** (the max count of the species at each
   deployment) by filtering, grouping, and summarising the sightings file.
5. **Join MaxN back onto the full effort list** with `left_join`, filling
   absent species with `MaxN = 0`. *True zeros must be kept; this is the single
   most common mistake when working with sighting data.*
6. **Look at the response** — mean, variance, table of counts.
7. **Choose the right model.** Counts ≠ continuous-normal data, so a standard
   `lm` is wrong. Use a **GLM with a Poisson or negative-binomial error
   distribution and a log link**.
8. **Fit a Poisson GLM** and **check for overdispersion** (Pearson dispersion
   statistic). If overdispersed, refit as **negative binomial** (`MASS::glm.nb`)
   and compare AIC.
9. **Test the full model against the null model** with a **likelihood-ratio
   χ² test** (`anova(null, full, test = "Chisq")`). Report χ², df, and p.
10. **Test each predictor** with `drop1(model, test = "Chisq")`.
11. **Report effect sizes** as **rate ratios** (`exp(coef)`) with confidence
    intervals (`exp(confint(model))`) — these are interpretable as
    "x times more individuals".
12. **Plot a figure** that shows the headline pattern.

## Result you should see

When you run `Rscript example_GLM_workflow.R` you should see (numbers will
match exactly because the simulation is seeded):

- **Final model:** Poisson GLM (Pearson dispersion ≈ 1.2 → Poisson is
  adequate; negative binomial gives no improvement).
- **Full vs null model:** χ² ≈ 77.3, df = 6, **p ≈ 1.3 × 10⁻¹⁴**.
- **Per-predictor LRTs:** `reef_zone` ***, `depth_m` ***, `season` ***
  significant; **`current_cms` is non-significant** — the deliberate
  red herring that teaches what a null result looks like.
- A figure of mean MaxN by reef zone × season in `outputs/`.

## License

MIT — see `LICENSE`.
