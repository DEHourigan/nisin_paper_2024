####################################
# Load the desired packages
####################################
.libPaths()
.libPaths( c( "/data/san/data0/users/david/rstudio/packages" , .libPaths() ) )
newlib <- "/data/san/data0/users/david/rstudio/packages"
.libPaths()
library(countrycode)
library(viridis)
library(ggsci)
library("scales")
library(data.table)
library(dplyr)
library(stringr)
library(janitor)
library(tibble)
library(lubridate)
library(countrycode)
library(ggplot2)
library(data.table)
library(tidyr)
library(ggnewscale)
pal_npg("nrc")(10)
library(networkD3)
library(plotly)
library(dplyr)
library(networkD3)
library(plotly)
library(ggsankey)
library(ape)

####################################
# Load the data
####################################
nisin_tsv = fread("/data/san/data1/users/david/mining_2023/nisins/tables/nisin_tsv.csv")
core_peptides = fread("/data/san/data1/users/david/mining_2023/nisins/tables/core_peptides.csv")
####################################
# ggtree of ALL organisms
####################################
library(ggtree)
tree = read.tree("../data/taxonomy/gtdbtk/tree/infer/intermediate_results/gtdbtk.unrooted.tree")
tree$tip.label
# rename tip labels by only taking first 15 characters
#tree$tip.label = sub("(.{15}).*", "\\1", tree$tip.label)
tree = ape::root.phylo(tree, 
 outgroup = "MGYG000121622")

tax_file <- fread('../data/taxonomy/gtdbtk/classify/gtdbtk.bac120.summary.tsv', 
  header = TRUE, sep="\t" ) %>%
  separate(classification, c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") 

####################################
############ TREE
####################################
tree_plot = ggtree(tree)
tree_plot$data
tree_plot_2 = tree_plot %<+% tax_file + 
	geom_tippoint(pch=15, aes(col=family), size=2) +
	#geom_hilight(mapping=aes(fill = as.factor(genus)), alpha = 0.5) +
	scale_color_manual(values = c(pal_npg("nrc")(10),pal_npg("nrc", alpha = 0.6)(10),pal_npg("nrc", alpha = 0.2)(10) )) +
	theme(axis.text.x = element_blank(), 
    axis.ticks = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.1, "cm")) +
	guides(col=guide_legend(ncol =1))

ggsave(tree_plot_2, 
	filename = "../figures/figure2_tree.jpg", 
	width = 5, 
	height = 5, 
	units = "in", 
	dpi = 600)



####################################
#  SANKEY plot
####################################

protein_dist_df = fread("/data/san/data1/users/david/mining_2023/nisins/tables/protein_dist_df.tsv")
unique(protein_dist_df$genus)
sankey_df = protein_dist_df %>%
  filter(n > 10) %>%
  filter(!is.na(genus)) %>% 
  filter(genus != "") %>%
  select(genus, species, protein_acc) %>%
#  distinct() %>%
  make_long(genus, protein_acc)


sankey_df <- sankey_df %>% 
  mutate(next_node = case_when(next_node == "WP_001058290.1" ~ "Suicin",
    next_node == "WP_015425978.1" ~ "Nisin-Z",
    next_node == "WP_081539212.1" ~ "Nisin-G-like",
    next_node == "WP_082026128.1" ~ "Geobacillin-X",
    next_node == "WP_032726862.1" ~ "Subtillin-like",
    next_node == "WP_079267797.1" ~ "novel",
    next_node == "WP_032726862.1" ~ "Entianin-like",
    next_node == "WP_003220055.1" ~ "Subtilin",
    next_node == "WP_210888916.1" ~ "Geobacillin-like",
    next_node == "WP_014570405.1" ~ "Nisin-A",
    next_node == "WP_073259420.1" ~ "Nisin-A-like",
    next_node == "WP_171299090.1" ~ "Nisin-A-like-2",
                               TRUE ~ next_node)) %>%
  mutate(node = case_when(node == "WP_001058290.1" ~ "Suicin",
    node == "WP_015425978.1" ~ "Nisin-Z",
    node == "WP_081539212.1" ~ "Nisin-G-like",
    node == "WP_082026128.1" ~ "Geobacillin-like",
    node == "WP_032726862.1" ~ "Subtillin-like",
    node == "WP_079267797.1" ~ "novel",
    node == "WP_032726862.1" ~ "Entianin-like",
    node == "WP_003220055.1" ~ "Subtilin",
    node == "WP_210888916.1" ~ "Geobacillin-like",
    node == "WP_014570405.1" ~ "Nisin-A",
    node == "WP_073259420.1" ~ "Nisin-A-like",
    node == "WP_171299090.1" ~ "Nisin-A-like",
                          TRUE ~ node))

# Convert 'node' and 'next_node' to factors
sankey_df$node <- factor(sankey_df$node)
sankey_df$next_node <- factor(sankey_df$next_node)

# Order the levels based on their frequency in descending order
sankey_df$node <- factor(sankey_df$node, levels = names(sort(table(sankey_df$node), decreasing = FALSE)))
sankey_df$next_node <- factor(sankey_df$next_node, levels = names(sort(table(sankey_df$next_node), decreasing = FALSE)))


core_peptide_sankey_plot = ggplot(sankey_df, 
  aes(x = x, 
    next_x = next_x, 
    node = node, 
    next_node = next_node, 
    fill = factor(node), 
    label = node)) +
  geom_sankey(flow.alpha = 0.7,
    color = "black",
    node.color = "black") +
  geom_sankey_label(size = 5, color = "black", fill="white", hjust = -0.2) +
  scale_fill_manual(values = c(
    pal_npg("nrc", alpha=1)(10),
    pal_npg("nrc", alpha = 0.8)(10),
    pal_npg("nrc", alpha = 0.6)(10),
    pal_npg("nrc", alpha = 0.2)(10),
    pal_npg("nrc", alpha = 0.1)(10))) +
  theme_sankey(base_size = 25) +
  theme(legend.position = "none",
    axis.text.x = element_blank(),
    plot.title = element_text(hjust = .5)) +
  xlab("") 

ggsave(core_peptide_sankey_plot, 
	file="../figures/figure_3_core_peptide_sankey_plot_greater_than_10.png", 
	width = 14, 
	height = 14)


################################
#  ST of potentially pathogenic strains
################################
st_df = fread("../data/tables/supplementary/supplementary_table_S5.csv", fill=TRUE)
sd_df = st_df %>% select(V1,V2,V3) %>%
  dplyr::rename(id = V1, sp = V2, st = V3) %>%
  group_by(sp, st) %>%
  summarise(n = n()) %>%
  filter(n >= 3 ) %>%
  filter(sp != "-") %>%
  filter(sp %in% c("bsubtilis", "sagalactiae", "ssuis", "saureus")) %>%
  mutate(sp = case_when(
    sp == "bsubtilis" ~ "B. subtilis",
    sp == "sagalactiae" ~ "S. agalactiae",
    sp == "ssuis" ~ "S. suis",
    sp == "saureus" ~ "S. aureus",
    TRUE ~ sp
  )) 


st_plot = ggplot(sd_df, aes(x = st, y = n)) +
  geom_bar(stat = "identity", width = 0.7, aes(fill=sp)) +
  #geom_violin(alpha = 0.5, fill="#21918c") +
  facet_wrap(~ sp, scales = "free_x", ncol = 4) +
  theme_bw() +
  labs(x = "Sequence Type (ST)", y = "Count", title = "") +
   theme(
    plot.title = element_blank(),  # Adjust size as needed
    axis.title.x = element_text(size = 8),  # Adjust size as needed
    axis.title.y = element_blank(),  # Adjust size as needed
    axis.text.x = element_text(size = 5, angle = 25),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 5)   # Adjust size as needed for y-axis labels
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(pal_npg("nrc")(10)))

ggsave(st_plot, 
  file = "../figures/supplementary_st_plot.png",
  width = 14, height = 4, units = "cm")