################################################################################
##  WORKED EXAMPLE: count-data GLM in R for ecological survey data
##  ---------------------------------------------------------------------------
##  Scenario (synthetic):
##    Underwater-video surveys for Nassau grouper (Epinephelus striatus,
##    species code "Epistr") on the Mesoamerican Barrier Reef. Each video
##    "deployment" produced a list of fish observations; we want to know how
##    MaxN of the focal species varies with:
##
##        reef_zone    (4 levels: Forereef, Lagoon, Backreef, Channel)
##        depth_m
##        current_cms  (current speed, cm/s)
##        season       (Dry / Wet)
##
##  The scenario is illustrative -- the point of this script is the WORKFLOW
##  for clean -> compute response -> join with effort -> choose model ->
##  test against the null model -> report. You can re-use the same pattern
##  on your own count data with different predictors.
##
##  Run from the repo root:   Rscript example_GLM_workflow.R
################################################################################

# 1. Set-up -------------------------------------------------------------------
library(dplyr)   # data wrangling
library(MASS)    # glm.nb()   (load AFTER dplyr; both define select())

# 2. Read the raw data --------------------------------------------------------
effort <- read.csv("data/effort.csv",    stringsAsFactors = FALSE)
sights <- read.csv("data/sightings.csv", stringsAsFactors = FALSE)
lookup <- read.csv("data/species_lookup.csv", stringsAsFactors = FALSE)

cat("\n--- raw data ---\n")
cat("effort rows  :", nrow(effort),  "\n")
cat("sighting rows:", nrow(sights),  "\n")
str(effort)

# 3. Clean the effort table ---------------------------------------------------
#
# 3a. DROP FAILED DEPLOYMENTS.
# A failed deployment (no usable video) is MISSING DATA, not a zero. If you
# treat it as MaxN = 0 you are pretending the species was absent when really
# the camera never recorded -- a serious bias. Identify failures from the
# notes column.
cat("\nNotes present on the effort sheet:\n")
print(table(effort$notes[effort$notes != ""]))

fail_pattern <- "flooded|battery|tipped|corrupted|fogged"
effort_ok <- effort[!grepl(fail_pattern, effort$notes, ignore.case = TRUE), ]
cat("Valid (non-failed) deployments:", nrow(effort_ok), "\n")

# 3b. Make factor variables and SET MEANINGFUL REFERENCE LEVELS.
# By default R orders factor levels alphabetically; that controls which level
# becomes the "intercept" in a model. Here we make Forereef and Dry the
# references so coefficients read as "compared to a Forereef site in the
# Dry season".
effort_ok$reef_zone <- factor(effort_ok$reef_zone,
                              levels = c("Forereef","Channel","Lagoon","Backreef"))
effort_ok$season    <- factor(effort_ok$season, levels = c("Dry","Wet"))

# (If your own data only has a DATE column and no season, you can derive a
# 2-level factor from the month, e.g.:
#   season <- factor(ifelse(format(as.Date(date), "%m") %in% c("12","01","02","03"),
#                           "Dry", "Wet")) )

# 4. Compute focal-species MaxN per deployment -------------------------------
#
# MaxN = maximum count of the focal species seen in any single observation
# row at a deployment. It is the standard relative-abundance metric for
# baited/un-baited video. Compute it by:
#   (i)  filter sightings to the focal species,
#   (ii) group by deployment and take max(n).
focal <- "Epistr"
focal_maxn <- sights %>%
  filter(species == focal) %>%
  group_by(location_code) %>%
  summarise(MaxN = max(n), .groups = "drop")

# 5. JOIN back to the full effort list ---------------------------------------
#
# *** This is the step students most often get wrong. ***
# Deployments where the focal species was NEVER seen are TRUE ZEROS -- they
# must stay in the dataset. They are not in the sightings file, so a naive
# summary of the sightings file silently drops them. Use left_join onto the
# full effort list and fill any NA with 0.
dat <- effort_ok %>%
  left_join(focal_maxn, by = "location_code") %>%
  mutate(MaxN = ifelse(is.na(MaxN), 0L, as.integer(MaxN)))

cat("\nMaxN distribution (zeros are real absences, kept by the join):\n")
print(table(dat$MaxN))

# 6. Explore the response ----------------------------------------------------
cat(sprintf("Mean MaxN = %.3f   Variance = %.3f   (var > mean -> check Poisson)\n",
            mean(dat$MaxN), var(dat$MaxN)))
print(dat %>% group_by(reef_zone, season) %>%
        summarise(mean_MaxN = round(mean(MaxN), 2), n = n(), .groups = "drop"))

# 7. Choose an appropriate model --------------------------------------------
#
# The response is COUNTS (0, 1, 2, ...). It cannot be negative, cannot be a
# fraction, and the variance generally grows with the mean. A standard
# LINEAR REGRESSION (lm) assumes a continuous, normally-distributed response,
# which is wrong here.
#
# The natural choice is a GENERALISED LINEAR MODEL (GLM) with a
#   - Poisson error distribution            (when var(Y) ~ mean(Y))
#   - or NEGATIVE BINOMIAL distribution      (when var(Y) > mean(Y))
# and a LOG link, so that the predicted mean stays positive.

# 8. Fit a Poisson GLM -------------------------------------------------------
m_pois <- glm(MaxN ~ reef_zone + depth_m + current_cms + season,
              data = dat, family = poisson(link = "log"))
cat("\n--- Poisson summary ---\n"); print(summary(m_pois))

# 8a. CHECK FOR OVERDISPERSION.
# Poisson assumes variance == mean. If the variance is much bigger
# (overdispersion), standard errors are too small and p-values too small.
# The Pearson dispersion statistic should be ~1.
disp <- sum(residuals(m_pois, type = "pearson")^2) / df.residual(m_pois)
cat(sprintf("Pearson dispersion = %.2f   (~1 ok, >~1.5 suggests NB)\n", disp))

# 8b. Refit as negative binomial and compare AIC.
m_nb <- glm.nb(MaxN ~ reef_zone + depth_m + current_cms + season, data = dat)
cat(sprintf("AIC: Poisson = %.1f   Negative binomial = %.1f\n",
            AIC(m_pois), AIC(m_nb)))

# Choose the simpler model unless NB clearly improves fit.
final <- if (disp > 1.5 && AIC(m_nb) + 2 < AIC(m_pois)) m_nb else m_pois
cat("FINAL MODEL: ",
    ifelse(identical(final, m_nb), "negative binomial GLM", "Poisson GLM"), "\n")

# 9. Test the FULL MODEL against the NULL MODEL ------------------------------
#
# The null model has only an intercept -- "no predictors at all". To find out
# whether our predictors jointly explain anything we use a LIKELIHOOD-RATIO
# TEST, which asymptotically follows a CHI-SQUARED distribution with df equal
# to the difference in the number of parameters between the two models.
m_null <- update(final, . ~ 1)   # null model in the SAME family as `final`
chi2 <- 2 * (logLik(final)[1] - logLik(m_null)[1])
df_d <- attr(logLik(final), "df") - attr(logLik(m_null), "df")
pval <- pchisq(chi2, df_d, lower.tail = FALSE)
cat(sprintf("\nFull vs null:  chi^2 = %.2f,  df = %d,  p = %.3g\n",
            chi2, df_d, pval))
cat("\n--- anova(null, full) ---\n")
print(anova(m_null, final, test = "Chisq"))

# 10. PER-PREDICTOR likelihood-ratio tests -----------------------------------
#
# drop1() drops each term in turn and reports the LR test for that drop.
# For GLMs this is generally preferred to the Wald z-tests in summary().
cat("\n--- drop1(): which predictors matter ---\n")
print(drop1(final, test = "Chisq"))

# 11. EFFECT SIZES: rate ratios = exp(beta) ----------------------------------
#
# Coefficients on a log link are MULTIPLICATIVE on the count. exp(beta) is
# the "rate ratio": the factor by which the mean MaxN multiplies for a
# one-unit increase in the predictor (or vs the reference category).
cat("\n--- Rate ratios exp(beta) with 95% CI ---\n")
rr <- exp(cbind(RateRatio = coef(final), suppressMessages(confint(final))))
print(round(rr, 3))

# 12. A figure for reporting ------------------------------------------------
if (!dir.exists("outputs")) dir.create("outputs")
png("outputs/MaxN_by_zone_and_season.png",
    width = 1200, height = 700, res = 150)
par(mar = c(7, 4, 3, 1))
tab <- with(dat, tapply(MaxN, list(reef_zone, season), mean))
barplot(t(tab), beside = TRUE, las = 2, col = c("grey70", "grey30"),
        ylab = "Mean MaxN",
        main = sprintf("Mean MaxN of %s by reef zone and season", focal))
legend("topright", legend = colnames(tab),
       fill = c("grey70", "grey30"), bty = "n")
dev.off()
cat("\nFigure saved -> outputs/MaxN_by_zone_and_season.png\n")

# 13. How to REPORT the result in a paper ------------------------------------
#
# Example sentence you could put in a Results section:
#
#   "The full Poisson GLM explained significantly more variation in Nassau
#    grouper MaxN than a null intercept-only model
#    (likelihood-ratio chi^2 = X.XX, df = D, p = P.PPP). MaxN differed
#    significantly between reef zones (chi^2 = X.X, df = 3, p = X.XXX) and
#    seasons (chi^2 = X.X, df = 1, p = X.XXX); current speed had no
#    detectable effect on MaxN (chi^2 = X.X, df = 1, p = X.XXX). MaxN was
#    approximately X.X times higher in the wet season than the dry season
#    (rate ratio X.X, 95% CI ..)."
#
# Always quote: the test statistic, df, p, the direction and size of the
# effect, and whether you used a Poisson or negative-binomial GLM.
cat("\n================ END OF WORKED EXAMPLE ================\n")
