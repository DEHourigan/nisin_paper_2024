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
library(rhmmer)
pal_npg("nrc")(10)

####################################
# Set WD
####################################
setwd("/data/san/data2/users/david/nisin/code")

####################################
# Load the data
####################################

nisin_tsv = fread("../data/tables/nisin_tsv.csv")
core_peptides = fread("../data/tables/core_peptides.csv")











################################
# ggtree of ALL organisms
################################

library(ggtree)
tree = read.tree("../data/tree/gtdbtk.unrooted.tree")
tree$tip.label
tree = ape::root.phylo(tree, 
 outgroup = "MGYG000121622")

tax_file <- fread('../data/tree/gtdbtk.bac120.summary.tsv', 
  header = TRUE, sep="\t" ) %>%
  filter(user_genome != "GCA_014488735.1_ASM1448873v1_genomic") %>% # this is the phage genome
  separate(classification, c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") 



################################
# taxa plots
################################

fig2_tax_df = tax_file %>% select(user_genome, family, order, genus, species)

fig2_tax_long <- fig2_tax_df %>%
  pivot_longer(cols = c(family, order, genus, species), 
               names_to = "rank", 
               values_to = "taxon")

# Count and filter for the top 10 in each rank
taxon_counts_top10 <- fig2_tax_long %>%
  group_by(rank, taxon) %>%
  summarize(count = n(), .groups = 'drop') %>%
  group_by(rank) %>%
  slice_max(order_by = count, n = 10)  # select the top 10 taxa by count within each rank

top_species = taxon_counts_top10 %>% filter(rank == "species")
top_genus = taxon_counts_top10 %>% filter(rank == "genus")
top_family = taxon_counts_top10 %>% filter(rank == "family")
top_ordder = taxon_counts_top10 %>% filter(rank == "order")

taxa_plot_top10_species = ggplot(top_species, aes(x = reorder(taxon, count), y = count, fill = taxon)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +  # Flip coordinates for better horizontal layout
  facet_wrap(~ rank, scales = "free_x", ncol = 4) +
  labs(x = NULL, y = NULL, title = "Top 10 Counts of Taxonomic Ranks") +
  theme_bw() +
   theme(
    plot.title = element_blank(),  # Adjust size as needed
    axis.title.x = element_text(size = 8),  # Adjust size as needed
    axis.title.y = element_blank(),  # Adjust size as needed
    axis.text.x = element_text(size = 5, angle = 25),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 7),
    plot.margin = margin(0, 0.02, 0, 0.4, "cm")   # Adjust size as needed for y-axis labels
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(pal_npg("nrc")(10), pal_npg("nrc", alpha = 0.6)(10), pal_npg("nrc", alpha = 0.2)(10)))

taxa_plot_top10_genus = ggplot(top_genus, aes(x = reorder(taxon, count), y = count, fill = taxon)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +  # Flip coordinates for better horizontal layout
  facet_wrap(~ rank, scales = "free_x", ncol = 4) +
  labs(x = NULL, y = NULL, title = "Top 10 Counts of Taxonomic Ranks") +
  theme_bw() +
   theme(
    plot.title = element_blank(),  # Adjust size as needed
    axis.title.x = element_text(size = 8),  # Adjust size as needed
    axis.title.y = element_blank(),  # Adjust size as needed
    axis.text.x = element_text(size = 5, angle = 25),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 7),
    plot.margin = margin(0, 0.02, 0, 0.4, "cm")   # Adjust size as needed for y-axis labels
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(pal_npg("nrc")(10), pal_npg("nrc", alpha = 0.6)(10), pal_npg("nrc", alpha = 0.2)(10)))

taxa_plot_top10_family = ggplot(top_family, aes(x = reorder(taxon, count), y = count, fill = taxon)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +  # Flip coordinates for better horizontal layout
  facet_wrap(~ rank, scales = "free_x", ncol = 4) +
  labs(x = NULL, y = NULL, title = "Top 10 Counts of Taxonomic Ranks") +
  theme_bw() +
   theme(
    plot.title = element_blank(),  # Adjust size as needed
    axis.title.x = element_text(size = 8),  # Adjust size as needed
    axis.title.y = element_blank(),  # Adjust size as needed
    axis.text.x = element_text(size = 5, angle = 25),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 7),
    plot.margin = margin(0, 0.02, 0, 0.4, "cm")   # Adjust size as needed for y-axis labels
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(pal_npg("nrc")(10), pal_npg("nrc", alpha = 0.6)(10), pal_npg("nrc", alpha = 0.2)(10)))

taxa_plot_top10_order = ggplot(top_ordder, aes(x = reorder(taxon, count), y = count, fill = taxon)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +  # Flip coordinates for better horizontal layout
  facet_wrap(~ rank, scales = "free_x", ncol = 4) +
  labs(x = NULL, y = NULL, title = "Top 10 Counts of Taxonomic Ranks") +
  theme_bw() +
   theme(
    plot.title = element_blank(),  # Adjust size as needed
    axis.title.x = element_text(size = 8),  # Adjust size as needed
    axis.title.y = element_blank(),  # Adjust size as needed
    axis.text.x = element_text(size = 5, angle = 25),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 7),
    plot.margin = margin(0, 0.02, 0, 0.4, "cm")   # Adjust size as needed for y-axis labels
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(pal_npg("nrc")(10), pal_npg("nrc", alpha = 0.6)(10), pal_npg("nrc", alpha = 0.2)(10)))


# ggsave(taxa_plot_top10_species, 
#   file = "../figures/taxa_plot_top10_species.png",
#   width = 6, height = 5, units = "cm")
# ggsave(taxa_plot_top10_genus, 
#   file = "../figures/taxa_plot_top10_genus.png",
#   width = 6, height = 5, units = "cm")
# ggsave(taxa_plot_top10_family, 
#   file = "../figures/taxa_plot_top10_family.png",
#   width = 6, height = 5, units = "cm")
# ggsave(taxa_plot_top10_order,
#   file = "../figures/taxa_plot_top10_order.png",
#   width = 6, height = 5, units = "cm")

library(gridExtra)

figure2_taxa_plot = grid.arrange(taxa_plot_top10_species, taxa_plot_top10_genus, taxa_plot_top10_family, taxa_plot_top10_order, 
  ncol = 4, nrow = 1, widths = c(1.2, 1, 1, 1)
)

ggsave(figure2_taxa_plot,
  file = "../figures/figure2_taxa_plot.png",
  width = 22, height = 6, units = "cm", dpi=600)

n_distinct(fig2_tax_df$family)
n_distinct(fig2_tax_df$genus)
n_distinct(fig2_tax_df$species)



################################
# ggtree of ALL organisms
################################

library(ape)
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

################################
#  plot a subset of the tree and heatmap of nisin genes they have
################################

tax_file <- fread('../data/tree/gtdbtk.bac120.summary.tsv', 
  header = TRUE, sep="\t" ) %>%
  separate(classification, c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  select(user_genome, family, genus, species) 
small_annot_df = fread("../data/tables/annot_df1.txt", header = TRUE, sep = ",") %>%
select(-prediction) %>% distinct()
small_annot_df_tax = left_join(small_annot_df, tax_file, by = c("filename" = "user_genome")) 


small_annot_df_tax %>%
  group_by(genus, filename) %>%
  summarise(type_count = n_distinct(type)) %>%
  ungroup() %>%
  arrange(genus, desc(type_count)) %>%
  group_by(genus) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(filename = paste0(filename, ".fna")) %>%
  select(filename) %>%
  fwrite("../data/tables/reduced_tree.txt", sep = "," , col.names = FALSE, row.names = FALSE)

################################
# ggtree of producing ogranisms
################################
tree = read.tree("../data/tree/gtdbtk.unrooted_subset.tree")
tree_tip_labels_old  <- tree$tip.label
tree$tip.label  <- ifelse(substr(tree$tip.label, 1, 1) == "G", sub("^(G[^_]+_[^_]+)_.*$", "\\1", tree$tip.label), tree$tip.label)

# df from hungate genomes 2023
tax_file <- fread('../data/tree/gtdbtk.bac120.summary.tsv', 
  header = TRUE, sep="\t" ) %>%
  separate(classification, c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") 

# taxonomy for subset
small_tax = tax_file %>% 
  filter(user_genome %in% tree_tip_labels_old) %>%
  select( user_genome, family, genus, species) %>%
  distinct() %>%
  mutate(assembly = if_else(str_starts(user_genome, "G"), 
                               str_replace(user_genome, "(^G[^_]+_[^_]+)_.*", "\\1"), 
                               user_genome)) %>%
  select(assembly, genus)


# read in nisin tsv 
tsv_df = fread("../data/tables/nisin_tsv.csv")

# hmmer
hmmer_wgs_dir ="../data/tables/hmm"
hmmer_wgs_files = list.files(hmmer_wgs_dir, pattern = "*hmm.out", recursive = FALSE, full.names = TRUE)
hmmer_wgs_files = grep("checker", hmmer_wgs_files, invert = TRUE, value = TRUE)
hmmer_wgs_df = do.call(rbind, lapply(hmmer_wgs_files, read_tblout)) %>%
  janitor::clean_names() %>%
  dplyr::rename(protein_acc = domain_name,
  type = query_name,
  evalue = sequence_evalue) %>%
  filter(case_when(
    type == "nisB" ~ evalue < 1e-02,
    type == "nisC" ~ evalue < 1e-02,
    type == "nisins" ~ evalue < 10,
    type == "nisT" ~ evalue < 1e-15,
    type == "nisF" ~ evalue < 1e-20,
    type == "nisI" ~ evalue < 1e-5,
    type == "nisE" ~ evalue < 1e-20,
    type == "nisG" ~ evalue < 1e-20,
    type == "nisR" ~ evalue < 1e-100,
    type == "nisk" ~ evalue < 1e-100,
    type == "nisP" ~ evalue < 1e-10,
    TRUE ~ FALSE)) %>%
   distinct(protein_acc, type, .keep_all = TRUE) %>%
   mutate(protein_acc = if_else(str_detect(protein_acc, "Bakta"), 
                               str_extract(protein_acc, "(?<=\\|Bakta\\|)[^|]+$"), 
                               protein_acc))

hmmer_wgs_master = hmmer_wgs_df %>%
  select("protein_acc","type","evalue")
assemly_protein_df = tsv_df %>%
  select(assembly, protein_acc) %>% 
  distinct()
small_names_df = tsv_df %>%
	filter(assembly %in% tree$tip.label) %>%
  left_join(., hmmer_wgs_master, by = "protein_acc")

wide_df <- small_names_df %>%
  pivot_wider(
    id_cols = c(assembly), 
    names_from = type, 
    values_from = type,
    values_fill = list(type = NA), # Fill absent types with NA
    values_fn = list(type = function(x) unique(x)) # Insert the type name
  ) %>%
  filter(assembly %in% tree$tip.label)


wide_df %>%
  filter(assembly == "GCF_000407445.1") # moraviensis

wide_df = wide_df %>% # these have cores predicted by rodeo that are not in the TSV files
  mutate(
  nisins = case_when(
    assembly == "GCF_003391095.1" ~ "nisins", 
    assembly == "GCF_008369905.1" ~ "nisins",
    assembly == "GCF_029395415.1" ~ "nisins",
    assembly == "GCF_015711475.1" ~ "nisins",
    TRUE ~ nisins
  ))

############  annot df
annot_df = column_to_rownames(wide_df , var = "assembly")  %>%
  as.data.frame() %>%
  select(nisins, nisB, nisC, nisI)

annot_df$nisins <- gsub("nisins", "core", annot_df$nisins)
annot_df$nisins <- as.factor(annot_df$nisins) # safe because all have a core peptide
annot_df$nisB <- as.factor(annot_df$nisB)
annot_df$nisC <- as.factor(annot_df$nisC)
# annot_df$nisR <- as.factor(annot_df$nisR)
# annot_df$nisk <- as.factor(annot_df$nisk)
# annot_df$nisF <- as.factor(annot_df$nisF)
# annot_df$nisE <- as.factor(annot_df$nisE)
# annot_df$nisG <- as.factor(annot_df$nisG)
# annot_df$nisT <- as.factor(annot_df$nisT)
annot_df$nisI <- as.factor(annot_df$nisI)


############ TREE
# ape root tree to Bacteroides
tree = ape::root.phylo(tree, 
 outgroup = "MGYG000121622")

tree_plot = ggtree(tree)



small_tax$genus <- gsub("g__Phocaeicola","g__Phocaeicola*",small_tax$genus)
small_tax$genus <- gsub("g__Bacteroides","g__Bacteroides*",small_tax$genus)

tree_plot_2 = tree_plot %<+% small_tax + 
  geom_tiplab(aes(label=genus), offset=0.5, hjust=0, size=3.2, align=TRUE) +
  scale_color_manual(values = carto_pal("ag_GrnYl", n = n_distinct(small_tax$genus))) 


# Create the heatmap
p1 <- gheatmap(tree_plot_2, annot_df, 
      legend_title = "Presence/Absence", 
			color = "black",
      colnames_angle=-45,
      colnames_offset_y = -0.5, 
			font.size = 2.7,
      colnames_offset_x = 0.03,
      offset=1.04,
      width=0.3,) +
      scale_fill_manual(
    values = c(pal_npg("nrc")(10), pal_npg("nrc", alpha = 0.4)(4)),
    breaks = unique(annot_df$values),  # Set unique values from annot_df as breaks
    na.value = "white"  # Set NA values to white
  ) +
  theme(legend.position = "bottom",
    axis.title.x = element_text(size = 10),  
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_blank(),  
    axis.text.y = element_blank(),  
    legend.title = element_text(size = 8),  
    legend.text = element_text(size = 8),  
    legend.key.size = unit(0.5, 'cm'),
    legend.key = element_rect(colour = "black",  size = 0.5)) +
    coord_cartesian(clip = "off")  # Adjust ratio as needed
    
ggsave(p1,
  file = "../figures/ggtree_heatmap.png", 
  width = 25, height = 25, limitsize = FALSE,dpi=600, units="cm")
 





# Create a distinct dataset with assembly and type
distinct_genomes <- hmm_table_genomes %>%
  	select(assembly, type) %>%
  	distinct()

tsv_df %>%
	distinct(assembly) %>%
	filter(!assembly %in% distinct_genomes$assembly)

phage_id="GCA_014488735.1"

distinct_genomes %>%
	filter(assembly == phage_id)



# Count total number of assemblies
total_assemblies <- n_distinct(tsv_df$assembly)

# Calculate the number of assemblies with nisB, nisC, nisI, and both
summary_table <- distinct_genomes %>%
	# mutate(assembly = gsub("GCA_", "GCF_",assembly)) %>%
	group_by(assembly) %>%
	summarize(
		nisB = any(type == "nisB"),
		nisC = any(type == "nisC"),
		nisI = any(type == "nisI")
	) %>%
	summarize(
		nisB_count = sum(nisB),
		nisC_count = sum(nisC),
		nisI_count = sum(nisI),
		both_count = sum(nisB & nisC & nisI),
		nisI_percentage = (nisI_count / total_assemblies) * 100
	) %>%
	mutate(
		nisB_percentage = (nisB_count / total_assemblies) * 100,
		nisC_percentage = (nisC_count / total_assemblies) * 100,
		nisI_percentage = (nisI_count / total_assemblies) * 100,
		both_percentage = (both_count / total_assemblies) * 100,
	) %>%
	select(nisB_count, nisB_percentage, nisC_count, nisC_percentage, nisI_count, nisI_percentage, both_count, both_percentage) %>%
  mutate(total_count = total_assemblies)


lanb <- hmmer_wgs_master[hmmer_wgs_master$type == "nisB",]$protein_acc
writeLines(lanb, "/data/san/data2/users/david/nisin/data/lanBCI_analysis/lanb.txt")

lanc <- hmmer_wgs_master[hmmer_wgs_master$type == "nisC",]$protein_acc
writeLines(lanc, "/data/san/data2/users/david/nisin/data/lanBCI_analysis/lanc.txt")

lani <- hmmer_wgs_master[hmmer_wgs_master$type == "nisI",]$protein_acc
writeLines(lani, "/data/san/data2/users/david/nisin/data/lanBCI_analysis/lani.txt")



# Load necessary libraries
library(data.table)
library(dplyr)
library(tidyr)

# Read the data without row names
core_pid <- fread("/data/san/data2/users/david/nisin/data/core_peptides/pid/total_unique_core.pid", header = TRUE, data.table = FALSE)

# Ensure column names are unique
colnames(core_pid) <- make.unique(colnames(core_pid))

# Convert the data to long format
core_pid_long <- core_pid %>%
  pivot_longer(cols = -V1, names_to = "V2", values_to = "V3")

# Perform the filtering operation
# WP_195320622.1 = ce02
query_Core = c(
	"NZ_WHVU01000008.1_17530_17713", # nisin J
	"WP_014570405.1", # 
	"WP_117640235.1", # nisin O
	"WP_008881441.1",
	"WP_003220055.1", # subtilin
	"WP_195320622.1")
filtered_core_pid <- core_pid_long %>%
  filter(V1 %in% query_Core & V2 %in% query_Core) %>%
  filter(V1 != V2) %>%
  arrange(V3)

# Print the filtered data
print(filtered_core_pid, n = Inf)

x = core_pid_long %>%
  filter(V1 == "NZ_WHVU01000008.1_17530_17713")
View(x)
