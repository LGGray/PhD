ggplot(data = ora.xcape,
       aes(axis1 = celltype, axis2 = ID, y = -log10(qvalue))) +
  geom_alluvium(aes(fill = ID),
                curve_type='arctangent',
                width = 1/18) +
  geom_stratum(width = 1/18) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)), size=3,
            hjust="inward") +
  scale_x_discrete(limits = c("Celltype", "ID"),
                   expand = c(0.05, 0.05)) +
  scale_fill_viridis_d() +
  theme_void() +
  theme(legend.position = 'none') 
