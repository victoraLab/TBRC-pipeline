# Victora-lab standard ggplot2 themes.
#
# Canonical definitions live in the lab Obsidian vault at
# "Science and Programming/Code Snippets/ggplot2 themes.md". They are vendored
# here because the QC figures ship inside the archive users download, so they
# must render on the cluster without depending on anything outside this repo.
#
# These are partial theme() objects: add them AFTER a base theme.

suppressPackageStartupMessages(library(ggplot2))

# Helvetica is not guaranteed on a headless cluster node. Fall back to the
# device default rather than emitting a font warning per figure.
lab_font <- local({
  available <- tryCatch(
    "Helvetica" %in% systemfonts::system_fonts()$family,
    error = function(...) TRUE
  )
  if (isTRUE(available)) "Helvetica" else ""
})

gvic_theme <- theme(
  axis.text        = element_text(size = 16, family = lab_font, colour = "black"),
  axis.title       = element_text(size = 16, family = lab_font, colour = "black"),
  legend.title     = element_text(size = 16, colour = "black"),
  legend.text      = element_text(size = 16, colour = "black"),
  plot.title       = element_text(size = 16, family = lab_font),
  panel.grid       = element_line(colour = "grey92"),
  panel.background = element_rect(fill = "white", colour = NA)
)

vicLabTheme <- theme(
  axis.text        = element_text(size = 8, family = lab_font, colour = "black",
                                  angle = 90, hjust = 1, vjust = 0.5),
  axis.title       = element_text(size = 16, family = lab_font, colour = "black"),
  legend.title     = element_text(size = 16, colour = "black"),
  legend.text      = element_text(size = 16, colour = "black"),
  plot.title       = element_text(size = 16, family = lab_font),
  panel.grid       = element_line(colour = "grey92"),
  panel.background = element_rect(fill = "white", colour = NA)
)

gvic_theme_violin <- gvic_theme + theme(
  axis.text.x = element_text(angle = 45, hjust = 1)
)

# Plate maps are drawn by ggplate, which sets its own geometry. Only the text and
# background elements from gvic_theme apply cleanly on top.
gvic_theme_plate <- theme(
  plot.title       = element_text(size = 16, family = lab_font, colour = "black"),
  legend.title     = element_text(size = 12, colour = "black"),
  legend.text      = element_text(size = 12, colour = "black"),
  panel.background = element_rect(fill = "white", colour = NA)
)
