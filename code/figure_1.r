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
library(cowplot)
pal_npg("nrc")(10)
setwd("/data/san/data2/users/david/nisin/code")

nisin_meta <- fread("../data/tables/supplementary/supplementary_table_S1.csv") 

nisin_countres_all = c(nisin_meta$country_value)
world_coordinates <- map_data("world") 

country_counts <- table(nisin_countres_all)
country_counts_df <- as.data.frame(country_counts)
colnames(country_counts_df) <- c("region", "count")
world_coordinates <- left_join(world_coordinates, country_counts_df, by = "region")

# Create bins for the counts
world_coordinates$count_bin <- cut(world_coordinates$count, 
                   breaks = c(0, 10, 50, Inf), 
                   labels = c("0-10", "10-50", ">50"), 
                   include.lowest = TRUE)

# Update the world plot to use a manual fill based on the bins
world_plot <- ggplot() + 
  geom_map(data = world_coordinates, map = world_coordinates, color = "black", 
       aes(long, lat, map_id = region, fill = count_bin)) + 
  scale_fill_manual(values = c(
  "0-10" = "#b8c5e0",
  "10-50" = "#5f7cb9",
  ">50" = "#3c5488"
  ), na.value = "grey88") + 
  theme_minimal() +
  labs(fill = "Count", font = "Arial") + 
  theme(
  legend.position = "bottom",
  legend.margin = margin(-5, 15, 15, 15),
  legend.text = element_text(size = 10),
  legend.box.margin = margin(0, 0, 0, 0),
  plot.background = element_rect(fill = "white"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank()
  )

ggsave(world_plot, file = "../figures/world_map_colored.png", 
  width = 18, height = 10, units = "cm", dpi = 600)

# SET PIE COLOURS
pie_colors = c(pal_npg("nrc", alpha = 0.6)(10), pal_npg("nrc")(10), pal_npg("nrc", alpha = 0.2)(10))


#################################
# bar charts
#################################
isolation_source_plot <- nisin_meta %>%
  count(isolation_source) %>%
  filter(!is.na(isolation_source) & isolation_source != "na") %>%
  arrange(desc(n)) %>%
  slice_max(n, n = 10) %>%
  ggplot(aes(x = reorder(isolation_source, n), y = n, fill = isolation_source)) +
  geom_bar(stat = "identity", color="black") +
  coord_flip(clip = "off") +
  labs(x = NULL, y = NULL, title = "Source") +
  theme_classic() +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.ticks = element_blank(),
        plot.title = element_text(face = "bold")) +
  geom_text(aes(label = n), hjust = -0.3, size = 4) +
  scale_fill_npg() +
  expand_limits(y = 350) 

# Create the host attribute value plot
host_attribute_value_plot <- nisin_meta %>%
  filter(host_attribute_value != "") %>%
  count(host_attribute_value) %>%
  filter(!is.na(host_attribute_value) & host_attribute_value != "na" & n >= 4) %>%
  arrange(desc(n)) %>%
  slice_max(n, n = 10) %>%
  ggplot(aes(x = reorder(host_attribute_value, n), y = n, fill = host_attribute_value)) +
  geom_bar(stat = "identity", color="black") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = "Host") +
  theme_classic() +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.ticks = element_blank(),
        plot.title = element_text(face = "bold")) +
  geom_text(aes(label = n), hjust = -0.3, size = 4) +
  scale_fill_npg() +
  expand_limits(y = 550)

# Create the country value plot
country_value_plot <- nisin_meta %>%
  count(country_value) %>%
  filter(!is.na(country_value) & country_value != "na") %>%
  arrange(desc(n)) %>%
  slice_max(n, n = 10) %>%
  ggplot(aes(x = reorder(country_value, n), y = n, fill = country_value)) +
  geom_bar(stat = "identity", color="black") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = "Country") +
  theme_classic() +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 10),
        axis.ticks = element_blank(),
        plot.title = element_text(face = "bold")) +
  geom_text(aes(label = n), hjust = -0.3, size = 4) +
  scale_fill_npg() +
  expand_limits(y = 450)

# Combine the plots side by side using cowplot
combined_plot <- plot_grid(
  isolation_source_plot, 
  host_attribute_value_plot, 
  country_value_plot, ncol = 3, align = "h")

# Save the combined plot
ggsave(combined_plot, 
  filename = "../figures/figure_1_bcd.png", 
  width = 18, height = 6, units = "cm", dpi = 600)



#################################
# get meta data for assemblies
#################################
library(lubridate)
epidermin = fread("../data/tables/supplementary/supplementary_epidermin.csv", 
  header=FALSE) 

meta_df = fread("../data/tables/assemblies_meta.tsv") %>%
  clean_names() %>%
  mutate(assembly_accession = sub("GCF_", "GCA_", assembly_accession)) %>%
  filter(!assembly_accession %in% epidermin$V3)

meta_df$seq_release_format <- as.Date(meta_df$assembly_bio_sample_submission_date, format="%Y%M")
meta_df$year <- floor_date(meta_df$seq_release_format, "year")

monthly_counts <- meta_df %>%
  group_by(year, assembly_accession) %>% 
  summarise(total_genomes = n_distinct(assembly_accession)) %>%
  ungroup() %>%
  mutate(total_genomes = cumsum(total_genomes))  %>%
  filter(!assembly_accession %in% epidermin$V3)

# Calculate the average number of genomes per year in total
overall_avg <- meta_df %>%
  select(assembly_accession, year) %>%
  distinct() %>%
  group_by(year) %>%
  summarise(n = n()) %>%
  summarise(overall_avg = mean(n, na.rm = TRUE))

# Calculate the average number of genomes per year for the last 5 years
last_5_years_avg <- meta_df %>%
  select(assembly_accession, year) %>%
  distinct() %>%
  group_by(year) %>%
  summarise(n = n()) %>%
  filter(year >= as.Date("2018-01-01")) %>%
  summarise(last_5_years_avg = mean(n, na.rm = TRUE))


overall_avg
last_5_years_avg

date_plot = ggplot(monthly_counts) +
  geom_line(aes(x = year, y = total_genomes)) +
  theme_bw() +
  labs(x = "", y = "", title = "Release Date") +
   theme(
    axis.title.x = element_text(size = 7),  
    axis.title.y = element_blank() ,  
    axis.text.x = element_text(size = 7),   
    axis.text.y = element_text(size = 7) ,
    plot.title = element_text(face = "bold",size = 10)   
   ) +
  theme(legend.position = "none")


#################################
# plot the counts of assembly levels
#################################
assembly_level_counts = meta_df %>% # replace Chromosome with Complete Genome
  mutate(assembly_level = ifelse(assembly_level == "Chromosome", "Complete", assembly_level)) %>%
  mutate(assembly_level = ifelse(assembly_level == "Complete Genome", "Complete", assembly_level)) %>%
  group_by(assembly_level) %>%
  summarise(count = n_distinct(assembly_accession)) %>%
  ungroup() %>%
  mutate(assembly_level = factor(assembly_level, levels = c("Complete", "Chromosome", "Scaffold", "Contig")))

assembly_level_plot = ggplot(assembly_level_counts) +
  geom_bar(aes(x = assembly_level, y = count), stat = "identity", 
    fill = "#2A788EFF",
    color = "black") +
  theme_bw() +
  labs(x = "", y = "Count", title = "Assembly level") +
   theme(
    plot.title = element_text(size = 10, face="bold"),
    axis.title.x = element_text(size = 7),  
    axis.title.y = element_text(size= 7) ,  
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7)
  )

#################################
# genome size vs protein coding sequence count
#################################
genome_size_vs_protein_coding_gene_count_plot = meta_df %>%
  filter(annotation_count_gene_protein_coding > 500) %>%
  ggplot() +
  geom_jitter(aes(x = assembly_stats_total_sequence_length, 
    y = annotation_count_gene_protein_coding, 
    fill= "#414487FF"), 
    size = 1, 
    color = "black", pch=21) +
  scale_fill_identity() +
  theme_bw() +
  labs(x = "Genome Size (bp)", y = "CDS Count", title = "Genome size vs CDS") +
   theme(
    plot.title = element_text(size = 10, face="bold"),
    axis.title.x = element_text(size = 7),  
    axis.title.y = element_text(size=7) ,  
    axis.text.x = element_text(size = 4),   
    axis.text.y = element_text(size = 7)
  ) +
  theme(legend.position = "none") 

####################################
# checkm2
####################################
checkm2_df=fread("../data/tables/quality_report.tsv") %>%
  clean_names()

comp_vs_cont = checkm2_df %>% 
  ggplot() +
  geom_jitter(aes(x = completeness, 
    y = contamination, 
    fill= "#414487FF"), 
    size = 1, 
    #alpha=0.5,
    color = "black", pch=21) +
  scale_fill_identity() +
  theme_bw() +
  labs(x = "Completeness (%)", y = "Contamination (%)", title = "Completeness") +
   theme(
    plot.title = element_text(size = 10, face="bold"),  
    axis.title.x = element_text(size = 7),  
    axis.title.y = element_text(size=7) ,  
    axis.text.x = element_text(size = 7),   
    axis.text.y = element_text(size = 7)    
  ) +
  theme(legend.position = "none")


# combine the 4 of these and plot
combined_plot = plot_grid(
  date_plot, 
  assembly_level_plot, 
  genome_size_vs_protein_coding_gene_count_plot, 
  comp_vs_cont, ncol = 4)

ggsave(combined_plot, 
  filename = "../figures/figure_1_e.png", 
  width = 18, height = 5, units = "cm", dpi = 600)



####################################
# descriptive statistics on the metadata
####################################
nisin_meta %>%
  filter(isolation_source == 'na' & host_attribute_value == 'na') %>%
  count(assembly)

# count all sources
nisin_meta %>%
  count(isolation_source) %>%
  arrange(desc(n))



table_species_per_source = nisin_meta %>%
  count(species, isolation_source) %>%
  arrange(desc(n))

nisin_meta %>%
  filter(grepl("diseased pig_unknown|systemic|infectioncerebrospinal|blood|brain|organ", isolation_source, ignore.case = TRUE)) %>%
  count(isolation_source,host_attribute_value,species) %>%
  filter(grepl("suis", species, ignore.case = TRUE)) %>%
  mutate( total = sum(n))

# count genera per source
table_host = nisin_meta %>%
  count(host_attribute_value) %>%
  arrange(desc(n))

table_species_per_host = nisin_meta %>%
  count(species, host_attribute_value) %>%
  arrange(desc(n))

# count species of s. suis
table_species_suis = nisin_meta %>%
  count(species, host_attribute_value) %>%
  filter(grepl("suis", species, ignore.case = TRUE)) %>%
  arrange(desc(n)) %>%
  mutate( total = sum(n)) %>%
  mutate(pig_total = sum(n))

# species table
table_species = nisin_meta %>%
  count(species) %>%
  arrange(desc(n))

unique(nisin_meta$species)
unique(nisin_meta$family)
unique(nisin_meta$genus)
