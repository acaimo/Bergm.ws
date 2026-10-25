.ws_pal <- list(
  surface  = "#fcfcfb",
  ink      = "#0b0b0b",
  ink2     = "#52514e",
  muted    = "#898781",
  grid     = "#e1e0d9",
  axis     = "#c3c2b7",
  neutral  = "#f0efec",
  pos      = "#2a78d6",
  neg      = "#e34948",
  critical = "#d03b3b",
  node     = "#52514e",
  node_cat = c("#1baf7a", "#eda100", "#6250d6"),
  node_other = "#b5b3ab",
  chains   = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
               "#e87ba4", "#008300", "#6250d6", "#e34948"),
  blue_ramp = c("#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec",
                "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab",
                "#184f95", "#104281", "#0d366b")
)

.ws_level_cols <- function(L, sign = c("pos", "neg")) {
  sign <- match.arg(sign)
  if (sign == "pos") {
    ramp <- colorRampPalette(c("#86b6ef", "#2a78d6", "#0d366b"), space = "Lab")
  } else {
    ramp <- colorRampPalette(c("#f0a39e", "#e34948", "#8e1f22"), space = "Lab")
  }
  ramp(max(L, 2))[seq_len(L)]
}

.ws_par <- function(...) {
  par(bg = .ws_pal$surface, fg = .ws_pal$axis, col.axis = .ws_pal$muted,
      col.lab = .ws_pal$ink2, col.main = .ws_pal$ink, font.main = 1,
      family = "sans", las = 1, tcl = -0.25, mgp = c(2.2, 0.45, 0), ...)
}
