################################################################################
##  EXAMPLE GLM SCRIPT: analysing count-data using a GLM for BRUV data in Cancun, Mexico
##  ---------------------------------------------------------------------------
##  Scenario:
##    BRUV survey for Nassau grouper (Epinephelus striatus,
##    species code "Epistr"). Each deployment
##    produced a list of fish sighitngs. We want to know how
##    MaxN of the focal species varies with:
##
##        reef_zone    (4 levels: Forereef, Lagoon, Backreef, Channel)
##        depth_m
##        current_cms  (current speed, cm/s)
##        season       (Dry / Wet)
##
##  The the point of this script is to:
##  clean data -> compute MaxN -> join with effort -> select model ->
##  test against the null model -> report. You can re-use the same workflow
##  on your own count data with different predictors.
##
##  
################################################################################

# 1. Set-up -------------------------------------------------------------------
library(dplyr)   # to manipulate data structure
library(MASS)    # glm.nb() for negative binomial (overdispersion)

# 2. Read the raw data --------------------------------------------------------
effort <- read.csv("data/effort.csv",    stringsAsFactors = FALSE)
sights <- read.csv("data/sightings.csv", stringsAsFactors = FALSE)
lookup <- read.csv("data/species_lookup.csv", stringsAsFactors = FALSE)

# Inspect data
nrow(effort)
nrow(sights)
table(species)


# 3. Clean the effort table ---------------------------------------------------
# DROP FAILED DEPLOYMENTS.
# A failed deployment (no usable video) is MISSING DATA, not a zero. If you
# treat it as MaxN = 0 you are treating the species as absent when actually
# the camera never recorded. Identify failures from the
# notes column.
table(effort$notes[effort$notes != ""])
table

# 3b. Make factors from characters, necessary to add as predictors
effort_ok$reef_zone <- factor(effort_ok$reef_zone,
                              levels = c("Forereef","Channel","Lagoon","Backreef"))
effort_ok$season    <- factor(effort_ok$season, levels = c("Dry","Wet"))


# 4. Compute MaxN per deployment of species of interest -------------------------------
# MaxN = maximum count of a species seen across a single deployment 
# Compute it by:
#   (i)  filter sightings to the species of interest,
#   (ii) group by deployment and take max(n).
focal <- "Epistr"
focal_maxn <- sights %>%
  filter(species == focal) %>%
  group_by(location_code) %>%
  summarise(MaxN = max(n), .groups = "drop")

# 5. JOIN back to the full effort list ---------------------------------------
#
# REMEMBER
# Deployments where the focal species was not seen are TRUE ZEROS: they
# must stay in the dataset. Use left_join onto the
# full effort list and fill any NA with 0.
dat <- effort_ok %>%
  left_join(focal_maxn, by = "location_code") %>%
  mutate(MaxN = ifelse(is.na(MaxN), 0L, as.integer(MaxN)))

#Inspect
table(dat$MaxN)

# 7. Choose an appropriate model --------------------------------------------
#
# The response is COUNTS (0, 1, 2, ...). It cannot be negative, cannot be a
# fraction, and the variance generally grows with the mean. A standard
# LINEAR REGRESSION (lm) assumes a continuous, normally-distributed response,
# which is wrong here.
#
# The first choice is a GENERALISED LINEAR MODEL (GLM) with a
#   Poisson error distribution            (when var(Y) ~ mean(Y))
#   or NEGATIVE BINOMIAL distribution      (when var(Y) > mean(Y)) (If overdispersion)
# and a LOG link function.

# 8. Fit a Poisson GLM -------------------------------------------------------
m_pois <- glm(MaxN ~ reef_zone + depth_m + current_cms + season,
              data = dat, family = poisson(link = "log"))
summary(m_pois)

# CHECK FOR OVERDISPERSION.
#Use dharma resuduals to test for overdispersion (easier to report than other methods)
import(DHARMa)
library(DHARMa)
output <- simulateResiduals(model = model, plot = TRUE)
output

# Refit as negative binomial and compare AIC.
# Choose the simpler model unless NB clearly improves fit.
m_nb <- glm.nb(MaxN ~ reef_zone + depth_m + current_cms + season, data = dat)
cat(sprintf("AIC: Poisson = %.1f   Negative binomial = %.1f\n",
            AIC(m_pois), AIC(m_nb)))


# 9. Test the FULL MODEL against the NULL MODEL ------------------------------
#
# The null model has only an intercept -- "no predictors at all". To find out
# whether our predictors jointly explain anything we use a LIKELIHOOD-RATIO
# TEST, which asymptotically follows a CHI-SQUARED distribution with df equal
# to the difference in the number of parameters between the two models.

m_null <- glm(MaxN ~ reef_zone + depth_m + current_cms + season,
              data = dat, family = poisson(link = "log"))
anova(m_full, m_null)

#REPORTING
# Always quote: the test statistic, df, p, the direction and size of the
# effect, and whether you used a Poisson or negative-binomial GLM.

