# Leaf temperature extraction from infrared snapshots
# Inputs: data/ir_example_snapshot.csv   (one radiometric snapshot)
#         data/ir_example_roi.csv        (Fiji ROI coordinates for that tray)
#
# Shown for a single snapshot from the 60 °C day. The full analysis applies the
# same extraction to every snapshot in each imaging window, producing
# data/fig3bcd_radiuses_leaves_combined_all.csv.

RADIUS <- 3     # pixels; ~2.9 mm at the imaging distance

# --- radiometric snapshot ----------------------------------------------------
# One temperature per detector pixel. Pixels outside the measurement range are
# recorded as -327.68.
snap <- read.csv("data/ir_example_snapshot.csv")
snap <- as.matrix(snap[-1])

# --- ROI coordinates ---------------------------------------------------------
# One point per tray cell, marked in Fiji on the TIFF converted from the same
# snapshot.
roi <- read.csv("data/ir_example_roi.csv")
names(roi)[1] <- "cell"
roi <- data.frame(cell = roi$cell, col = roi$X, row = roi$Y)

# --- pixels within a disc of RADIUS px around a point ------------------------
disc_pixels <- function(cx, cy, m, r) {
  xs <- max(1, floor(cx - r)):min(ncol(m), ceiling(cx + r))
  ys <- max(1, floor(cy - r)):min(nrow(m), ceiling(cy + r))
  g  <- expand.grid(col = xs, row = ys)
  g[sqrt((cx - (g$col + 0.5))^2 + (cy - (g$row + 0.5))^2) <= r, ]
}

# --- mean temperature per cell -----------------------------------------------
leaf_temp <- data.frame(
  record      = "Record_2025-04-16_10-38-00",
  cell        = roi$cell,
  area_mean_3 = mapply(function(cx, cy) {
    p <- disc_pixels(cx, cy, snap, RADIUS)
    mean(snap[cbind(p$row, p$col)], na.rm = TRUE)
  }, roi$col, roi$row)
)

write.csv(leaf_temp, "ir_extraction_example_temperatures.csv", row.names = FALSE)
print(head(leaf_temp))

# --- check: mask the frame to the sampled discs ------------------------------
# Written out as a PNG so ROI placement can be inspected against the image.
mask <- matrix(0, nrow(snap), ncol(snap))
for (i in seq_len(nrow(roi))) {
  p <- disc_pixels(roi$col[i], roi$row[i], snap, RADIUS)
  mask[cbind(p$row, p$col)] <- snap[cbind(p$row, p$col)]
}

ragg::agg_png("ir_extraction_example_visualization.png", width = 600, height = 600)
image(t(apply(mask, 2, rev)))
dev.off()