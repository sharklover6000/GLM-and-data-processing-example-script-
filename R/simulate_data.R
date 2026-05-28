################################################################################
##  Generate the example dataset used by example_GLM_workflow.R
##
##  Scenario: a synthetic (made-up) underwater-video survey for Nassau grouper
##  (Epinephelus striatus, code "Epistr") on the Mesoamerican Barrier Reef.
##  Real-world locations, but the data are fully simulated for teaching.
##
##  Run from the repo root:   Rscript R/simulate_data.R
################################################################################

set.seed(20260528)

n_sites <- 60
zones   <- c("Forereef", "Lagoon", "Backreef", "Channel")

sites <- data.frame(
  site_id   = sprintf("S%02d", seq_len(n_sites)),
  reef_zone = sample(zones, n_sites, replace = TRUE),
  depth_m   = round(runif(n_sites, 5, 22), 1)
)

seasons <- c("Dry", "Wet")
dep <- merge(sites, data.frame(season = seasons), by = NULL)
dep$location_code <- sprintf("BLZ_2025_%03d", seq_len(nrow(dep)))

## date column: Dry = early March, Wet = mid September
dep$date <- ifelse(
  dep$season == "Dry",
  as.character(as.Date("2025-03-01") + sample(0:14, nrow(dep), replace = TRUE)),
  as.character(as.Date("2025-09-15") + sample(0:14, nrow(dep), replace = TRUE))
)

## current speed (cm/s): varies by zone with noise
zone_current <- c(Forereef = 35, Lagoon = 15, Backreef = 20, Channel = 55)
dep$current_cms <- round(zone_current[dep$reef_zone] + rnorm(nrow(dep), 0, 10), 1)

## ----- TRUE data-generating model for the focal species ----------------------
## log(mu) = b0
##         + bZone[reef_zone]                    (ref = Forereef)
##         + bSeasonWet * (season == "Wet")
##         + bDepth5    * (depth_m - 12) / 5     (per +5 m)
##         + bCurrent   * (current_cms - 30)/10  (TRUE = 0  -> red herring)
beta0       <- log(0.7)
beta_zone   <- c(Forereef = 0, Channel = 0.55, Lagoon = -0.95, Backreef = -1.50)
beta_seasW  <- 0.55
beta_depth5 <- 0.35
beta_curr   <- 0.00            # not a real driver
theta_nb    <- 4               # mild overdispersion

mu <- exp(
  beta0
  + beta_zone[dep$reef_zone]
  + ifelse(dep$season == "Wet", beta_seasW, 0)
  + beta_depth5 * (dep$depth_m - 12) / 5
  + beta_curr   * (dep$current_cms - 30) / 10
)
maxn_focal <- rnbinom(nrow(dep), size = theta_nb, mu = mu)

## ----- inject some failed deployments (no usable footage) --------------------
fail_idx <- sample(seq_len(nrow(dep)), 5)
dep$notes <- ""
dep$notes[fail_idx] <- sample(c("Camera flooded", "Battery died early",
                                "Frame tipped on retrieval", "GoPro corrupted file",
                                "Lens fogged - unscoreable"), 5)
maxn_focal[fail_idx] <- NA      # failed -> no MaxN, no sightings

## ----- build EFFORT.csv (no MaxN; no Round/Season column? keep season since the date is the same form as date in exam — let's keep season explicit too, since season is a known field in real BRUV datasets) -----
effort <- dep[, c("location_code", "date", "season", "reef_zone",
                  "depth_m", "current_cms", "notes")]
write.csv(effort, "data/effort.csv", row.names = FALSE)

## ----- build SIGHTINGS.csv ---------------------------------------------------
nonfocal <- c("Chrmul","Sphbar","Lutjoc","Halrad","Carrub","Acabah",
              "Scavet","Stetri","Mycti","Holads")

sights <- list()
for (i in seq_len(nrow(dep))) {
  if (is.na(maxn_focal[i])) next                      # failed deployment

  ## focal-species rows (0-3 observations whose max == MaxN)
  if (maxn_focal[i] > 0) {
    n_obs <- sample(1:3, 1)
    counts <- if (maxn_focal[i] == 1) rep(1L, n_obs)
              else sample.int(maxn_focal[i], n_obs, replace = TRUE)
    counts[1] <- maxn_focal[i]                        # ensure max == MaxN
    counts <- sample(counts)                          # randomise order in video
    for (c in counts) {
      sights[[length(sights) + 1]] <- data.frame(
        location_code = dep$location_code[i],
        time_s = sample(0:3600, 1),
        species = "Epistr",
        n = c, notes = "",
        stringsAsFactors = FALSE
      )
    }
  }

  ## non-focal community
  n_other <- rpois(1, lambda = 5)
  for (j in seq_len(n_other)) {
    sights[[length(sights) + 1]] <- data.frame(
      location_code = dep$location_code[i],
      time_s = sample(0:3600, 1),
      species = sample(nonfocal, 1),
      n = sample(1:4, 1), notes = "",
      stringsAsFactors = FALSE
    )
  }
}
sights <- do.call(rbind, sights)
sights <- sights[order(sights$location_code, sights$time_s), ]
write.csv(sights, "data/sightings.csv", row.names = FALSE)

## ----- species lookup --------------------------------------------------------
lookup <- data.frame(
  code = c("Epistr", nonfocal),
  scientific_name = c(
    "Epinephelus striatus",   "Chromis multilineata",  "Sphyraena barracuda",
    "Lutjanus jocu",          "Halichoeres radiatus",  "Caranx ruber",
    "Acanthurus bahianus",    "Scarus vetula",         "Stegastes triangulatus",
    "Mycteroperca tigris",    "Holocentrus adscensionis"
  ),
  common_name = c(
    "Nassau grouper",     "Brown chromis",      "Great barracuda",
    "Dog snapper",        "Puddingwife wrasse", "Bar jack",
    "Ocean surgeonfish",  "Queen parrotfish",   "Threespot damselfish",
    "Tiger grouper",      "Squirrelfish"
  ),
  is_focal = c(TRUE, rep(FALSE, length(nonfocal)))
)
write.csv(lookup, "data/species_lookup.csv", row.names = FALSE)

cat("Wrote data/effort.csv     (", nrow(effort), "rows;",
    sum(effort$notes != ""), "failed )\n")
cat("Wrote data/sightings.csv  (", nrow(sights), "rows;",
    length(unique(sights$species)), "species )\n")
cat("Wrote data/species_lookup.csv\n")
cat("\nTrue focal MaxN summary (incl. NA for failures):\n")
print(summary(maxn_focal))
