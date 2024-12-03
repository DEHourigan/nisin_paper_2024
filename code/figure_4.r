####################################
# Load the desired packages
####################################
.libPaths()
.libPaths( c( "/data/san/data0/users/david/rstudio/packages" , .libPaths() ) )
newlib <- "/data/san/data0/users/david/rstudio/packages"
.libPaths()
library(thacklr)
library(gggenomes)
library(ggsci)
library(scales)
library(data.table)
library(dplyr)
library(stringr)
library(janitor)
library(tibble)
library(ggplot2)
library(tidyr)
library(ggnewscale)
library(msa)
library(gridExtra)
library(janitor)
library(rhmmer)


setwd("/data/san/data2/users/david/nisin/code")

tax_nuc = fread("../data/tables/supplementary/supplementary_table_S1.csv")
contigs_w_nisin = fread("../data/tables/supplementary/contigs.txt") 



platon_df = fread("../data/plasmids/all_plasmids.tsv") %>%
  janitor::clean_names() %>%
  mutate(nucleotide_acc = id) %>% select(nucleotide_acc) %>%
  filter(nucleotide_acc %in% contigs_w_nisin$nucleotide_acc)


tax_nuc$assembly %>% unique()

############################
# What is on MGEs - the plasmid content
############################
plasmids_all = fread("../data/tables/plasmids_all.tsv")
mlplasmids_df = fread("../data/tables/mlplasmids_df.tsv",
   sep="\t")
plasmids_all = left_join(plasmids_all, mlplasmids_df, by = "nucleotide_acc")

plasmids_10000to200000 = plasmids_all %>%
  filter(contig_length > 10000 & contig_length < 250000) %>%
  distinct()

  

plasmids_10000to200000_plot = ggplot(plasmids_10000to200000, 
    aes(x = "All plasimids", y = contig_length)) +
  geom_violin(alpha = 0.5, fill="#21918c") +
  geom_jitter(aes(
    shape = ifelse(nucleotide_acc %in% platon_df$nucleotide_acc, "star", "dot")),
    width = 0.3, size = 2, alpha = 0.5) +
  theme_bw() +
  labs(x = "", y = "Plasmid contig length", title = "Size") +
   theme(
    plot.title = element_text(size = 10),  # Adjust size as needed
    axis.title.x = element_text(size = 10),  # Adjust size as needed
    axis.title.y = element_text(size = 10),  # Adjust size as needed
    axis.text.x = element_text(size = 10),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 10)   # Adjust size as needed for y-axis labels
  ) +
  theme(legend.position = "none")

ggsave(plasmids_10000to200000_plot,
  file = "../figures/plasmids_10000to200000_plot.png",
  width = 6, height = 10, units = "cm")


plasmids_tax_plot = ggplot(plasmids_10000to200000, aes(x = contig_length, fill = genus)) +
  geom_histogram(binwidth = 20000, position = "stack", color="black") + # Set bin width to 20kb
  theme_bw() +
  labs(x = "Plasmid contig (in 20kb windows)", y = "Count", title = "Plasmid contig Distribution by Genus") +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10), 
    axis.text.y = element_text(size = 10)
  ) + # add text saying n=150
  geom_text(aes(x = 125000, y = 40, label = paste0("Count: ", n_distinct(plasmids_10000to200000$nucleotide_acc))), size = 4) +
  scale_fill_manual(values = c(
    pal_npg("nrc", alpha=1)(10),
    pal_npg("nrc", alpha = 0.4)(10))) +
  coord_cartesian(ylim = c(0, 10)) + # Adjust the y-axis limits as needed
  scale_y_continuous(breaks = seq(0, 10, 1), labels = function(x) as.integer(x)) # Display whole numbers on y-axis

ggsave(plasmids_tax_plot,
  file = "../figures/plasmids_taxonomy_plot.png",
  width = 20, height = 10, units = "cm")


plasmid_plot = grid.arrange(plasmids_10000to200000_plot, plasmids_tax_plot, 
  ncol = 2,
  widths = c(1, 5))

ggsave(plasmid_plot,
  file = "../figures/plasmid_double_plot.png", 
  width = 26, height = 10, units = "cm", dpi=600)


############################
# HGT of some nisin core peptides
############################

# note this part of the code requires all the annotated tables from bakta 
# and cannot be uploaded to github due to the size
hmm_dir = "../data/tables/hmm_bakta/"
hmm_files = list.files(hmm_dir, pattern = "*hmm.out", recursive = FALSE, full.names = TRUE)
hmm_df = do.call(rbind, lapply(hmm_files, read_tblout)) %>%
  janitor::clean_names() %>%
  dplyr::rename(protein_acc = domain_name,
  type = query_name,
  evalue = sequence_evalue)

hmm_core_BCIT = hmm_df %>%
  select("protein_acc","type","evalue") %>%
  filter(evalue < 0.0005) %>%
  filter(type %in% c("nisins", "nisB", "nisC", "nisI", "nisT")) %>%
  distinct()


fwrite(hmm_core_BCIT, "../tables/hmm_core_BCIT.tsv", sep="\t")
unique(hmm_df$type) # contig must contain core, nisB or nisC

############################
# 5 core peptides that are present in 5 seperate species
############################
# read in 6 genome GCA_014488735.1, 
hgt_assemblies = c("GCA_014488735.1","GCF_000644315.1", "GCF_024580395.1",
  "GCF_002077325.1","GCF_009771385.1","GCF_018424365.1")

tsv_dir = "../data/bakta_out"
# read in hgt.tsv files and bind
hgt_files = list.files(tsv_dir, full.names = TRUE, recursive = TRUE, pattern = ".tsv")
matching_files = grep(paste(hgt_assemblies, collapse = "|"), hgt_files, value = TRUE)
matching_files = grep("hypothetical", matching_files, invert = TRUE, value = TRUE)


contigs = c("NZ_KK234519.1", "NZ_JANILA010000008.1", "NZ_MJDR01000005.1", 
  "NZ_WNIR01000003.1", "NZ_VZNI01000014.1","MN270279.1")

hgt_df = do.call(rbind, lapply(matching_files, fread)) %>%
  janitor::clean_names() %>%
  mutate(protein_acc = locus_tag) %>%
  dplyr::rename(
    feat_id = locus_tag,
    gene_type = type,
    end = stop,
    gene_id = gene) %>%
  left_join(., hmm_core_BCIT, by = "protein_acc") %>%
  dplyr::rename(
    hmm_type = type,
    seq_id = number_sequence_id) %>%
  mutate(type = "CDS") %>%
  filter(seq_id %in% contigs) %>%
  mutate(bin_id = case_when(
    seq_id == "NZ_KK234519.1" ~ "S. aureus",
    seq_id == "NZ_JANILA010000008.1" ~ "S. hyicus", 
    seq_id == "NZ_MJDR01000005.1" ~ "E. hirae",
    seq_id == "NZ_WNIR01000003.1" ~ "S. agalactiae",
    seq_id == "NZ_VZNI01000014.1" ~ "S. suis",
    seq_id == "MN270279.1" ~ "phage-SsuZKB4",
    TRUE ~ seq_id)) %>%
  filter(seq_id %in% contigs) %>%
  filter(!type %in% c("nisI_and_PF18217")) %>%
  mutate(hmm_type = case_when(
    hmm_type == "nisins" ~ "core",
    hmm_type == "nisB" ~ "lanB",
    hmm_type == "nisI" ~ "lanI",
    hmm_type == "nisC" ~ "lanC",
    hmm_type == "nisT" ~ "lanT",
    TRUE ~ hmm_type))

hgt_df$strand <- case_when(
  hgt_df$strand == "+" ~ "+",
  hgt_df$strand == "-" ~ "-",
  hgt_df$strand == "?" ~ NA,
  TRUE ~ hgt_df$strand
)


hgt_df <- hgt_df %>%
  mutate(hmm_type = case_when(
    product == "Conjugative transposon protein" ~ "Conjugative transposon protein",
    product == "Leader peptide-processing serine protease" ~ "lanP",
    product == "Lantibiotic protection ABC transporter ATP-binding protein" ~ "lanFEG",
    product == "Lantibiotic protection ABC transporter permease protein" ~ "lanFEG",
    product == "ABC transporter permease" ~ "lanFEG",
    product == "Plasmid mobilization relaxosome protein MobC" ~ "MobC",
    product == "Two-component system, response regulator, lantibiotic associated" ~ "lanRK",
    product == "histidine kinase" ~ "lanRK",
    TRUE ~ hmm_type))


fwrite(hgt_df, file="../data/tables/gggenomes_table.tsv",
  sep="\t")
  
# add links and GC content
blast <- read_links("../data/tables/ava_nuc.o6")
gc <- fread("../data/tables/gc_out.tsv")

blast <- blast %>% #filter(seq_id %in% c("NZ_VZNI01000014.1","MN270279.1") & seq_id2 %in% c("NZ_VZNI01000014.1","MN270279.1")) %>% 
  filter(seq_id2 != "NZ_VZNI01000014.1" & seq_id != "MN270279.1")

################################
# HGT protein level
################################
prot_blast = read_sublinks("../data/tables/gggenomes_proteins.o6") %>%
  filter(pident > 30, length > 100)


hgt_plot2 = gggenomes(genes = hgt_df) %>%
  add_sublinks(links = prot_blast) %>%
  pick(2,6,5,1,4,3) %>%
  add_feats(gc) %>%
  flip(3,5,6) %>%
  focus(hmm_type == "core" | gene_id == "nsuT", 
    .expand = c(1.1e4, 1.5e4),
    .overhang = c("trim")) +
  geom_seq() +
  geom_link(aes(fill=pident, color=pident),
          offset = 0.01) +
  scale_color_gradient("%ID",low="grey88", high="grey50") +
  scale_fill_gradient("%ID",low="grey88", high="grey50") +
  ggnewscale::new_scale("fill") +
  geom_wiggle(aes(z=score), fill="#2367db", alpha=.5, offset=-0.3, height=.3) +
  geom_gene(aes(fill=hmm_type), size=8) +
  scale_fill_manual("gene",values = c(
    pal_npg("nrc", alpha=1)(10),
    pal_npg("nrc", alpha = 0.4)(10)), 
      na.value = "grey88") +
  geom_bin_label(expand_left = 0.3, size=6) +
  theme_gggenomes_clean(base_size = 20)

ggsave(hgt_plot2,
  file = "../figures/hgt_plot2.png", 
  width = 13,
  dpi = 600,
  height = 8)



################################
# gggenomes of all ICE
################################
ice_to_query <- c(
  "NZ_ATWZ01000005.1",
  "NZ_ALUP01000012.1",
  "NZ_AEYY01000049.1",
  "NC_015600.1",
  "NZ_BCCR01000016.1",
  "NZ_BCFB01000042.1",
  "NZ_CEDG01000006.1",
  "NZ_CKGJ01000002.1",
  "NZ_CP116958.1",
  "NZ_MJDR01000005.1"
)
fwrite(as.data.frame(ice_to_query), "../data/tables/ice_to_query.tsv", sep="\t",col.names = FALSE)

ice_bakta_df = fread("../data/tables/bakta_df.tsv", sep="\t") %>%
  filter(seq_id %in% ice_to_query) %>%
  janitor::clean_names() %>%
  mutate(protein_acc = locus_tag,
    feat_id = locus_tag) %>%
  dplyr::rename(
    gene_type = type,
    # end = stop,
    gene_id = gene) %>%
  left_join(., hmm_core_BCIT, by = "protein_acc") %>%
  dplyr::rename(
    hmm_type = type
   ) %>%
  mutate(type = "CDS") %>%
  filter(!type %in% c("nisI_and_PF18217")) %>%
  mutate(hmm_type = case_when(
    hmm_type == "nisins" ~ "core",
    hmm_type == "nisB" ~ "lanB",
    hmm_type == "nisI" ~ "lanI",
    hmm_type == "nisC" ~ "lanC",
    hmm_type == "nisT" ~ "lanT",
    TRUE ~ hmm_type))

ice_bakta_df$strand <- case_when(
  ice_bakta_df$strand == "+" ~ "+",
  ice_bakta_df$strand == "-" ~ "-",
  ice_bakta_df$strand == "?" ~ NA,
  TRUE ~ ice_bakta_df$strand
)


ice_plot1 = gggenomes(genes = ice_bakta_df) %>%
  focus(gene_id %in% c("nsuT"), .expand = c(2e4, 2e4)) +
  geom_seq() +
  geom_gene(aes(fill = hmm_type), size = 6) +
  scale_fill_manual(values = c(
    pal_npg("nrc", alpha = 1)(10), 
    pal_npg("nrc", alpha = 0.4)(10)
  )) +
  geom_bin_label(expand_left = 0.3, size = 4) +
  ggtitle("") +
  theme_gggenomes_clean(base_size = 12)  # Explicitly defining base_size


ggsave(ice_plot1,
  file = "../figures/ice_plot1.png", 
  width = 12,
  dpi=600,
  height = 4)

################################
# plot gene content of ICE
################################
n_ice = 145 # number of ice predicted with nisin core peptides

file_path <- "../data/ice/protein_content_analysis/in"
file_pattern <- "*.faa"
file_list <- list.files(path = file_path, pattern = file_pattern, full.names = TRUE)
all_data <- data.frame()
ice_contigs = gsub("_I.E.*", "", basename(file_list))

extract_data_from_file <- function(file) {
  lines <- readLines(file)
  lines <- grep(">", lines, value = TRUE)
  data <- data.frame(
    ICE = sapply(lines, function(x) strsplit(basename(file), "_")[[1]][1]),
    contig = sapply(lines, function(x) sub(":>.*", "", basename(file))),
    protein = sapply(lines, function(x) sub(".*:>(.*?) .*", "\\1", x)),
    ID = sapply(lines, function(x) sub(".* (TMPID_[0-9]+) .*", "\\1", x)),
    product = sapply(lines, function(x) sub(".* TMPID_[0-9]+ (.*)", "\\1", x))
  )
  return(data)
}

for (file in file_list) {
  file_data <- extract_data_from_file(file)
  all_data <- bind_rows(all_data, file_data)
}

ice_content <- all_data %>% select(product, contig)

# Count the number of unique contigs for each product
product_contig_count <- ice_content %>%
  distinct(product, contig) %>%
  group_by(product) %>%
  summarise(contig_count = n())

# Calculate the total number of unique contigs, this should equal the number of ICEs
total_contigs <- n_distinct(ice_content$contig)
total_contigs == n_ice

# Calculate the percentage of contigs for each product
product_contig_percentage <- product_contig_count %>%
  mutate(percentage = (contig_count / total_contigs) * 100) %>%
  filter(percentage > 50) %>%
  filter(product != "hypothetical protein")


nisin_products = c("ABC transporter ATP-binding protein", 
      "NisI/SpaI family lantibiotic immunity lipoprotein",
      "gallidermin/nisin family lantibiotic",
      "Lantibiotic",
      "Lantibiotic ABC transporter ATP-binding protein",
      "Lantibiotic ABC transporter permease",
      "Lantibiotic immunity ABC transporter MutG family permease subunit",
      "lanthionine synthetase C family protein",
      "lantibiotic ABC transporter permease",
      "lantibiotic dehydratase")
MGE_products = c("")

ice_content_plot = ggplot(product_contig_percentage, 
  aes(x = percentage, y = reorder(product, percentage))) +
  geom_bar(stat = "identity", color = "black", size=0.25, aes(
    fill = ifelse(product == "tetracycline resistance ribosomal protection protein Tet(O)", "#E64B35B2", 
    ifelse(product %in% nisin_products, "#3C5488B2", "gray")))) +
  theme_bw(base_size = 6, base_family = "Arial") +
  labs(x = "% of ICE elements", y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.8, vjust = 0), 
    plot.title = element_text(size = 12),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "bottom") +
  scale_fill_identity()

ggsave(ice_content_plot,
  file = "../figures/ice_content_plot.png",
  width = 10, height = 8, units = "cm", dpi=600)  

fwrite(product_contig_count, "../data/tables/supplementary/supplementary_table_S8.csv", sep="\t")

################################
# plot gene content of plasmids
################################
plasmid_df = fread("../data/tables/plasmid_content_df.tsv", sep="\t") %>%
  select(seq_id, product) %>%
   mutate(product = case_when(
      product == "Nisin biosynthesis protein NisB" ~ "lantibiotic dehydratase",
      product == "Lantibiotic dehydratase"~ "lantibiotic dehydratase",
      product == "NsjB lantibiotic dehydratase" ~ "lantibiotic dehydratase",
      product == "Lantibiotic modifying enzyme" ~ "lantibiotic dehydratase",
      product == "Lant-dehydr-N domain-containing protein" ~ "lantibiotic dehydratase",
      product == "thiopeptide-type bacteriocin biosynthesis protein" ~ "lantibiotic dehydratase",
      product == "(pseudo) lantibiotic dehydratase" ~ "lantibiotic dehydratase",
      product == "lanthionine synthetase C family protein" ~ "lanthionine synthetase C family protein",
      product == "Lanthionine synthetase C-like protein" ~ "lanthionine synthetase C family protein",
      product == "Lantibiotic cyclase" ~ "lanthionine synthetase C family protein",
      product == "NsjC Lanthionine synthetase C-like protein" ~ "lanthionine synthetase C family protein",
      TRUE ~ product
    )
  )

plasmid_df %>%
  group_by(seq_id) %>%
  filter(grepl("dehy", product, ignore.case = TRUE)) 

seq_ids_without_dehydratase <- plasmid_df %>%
  group_by(seq_id) %>%
  filter(!any(grepl("dehydratase", product, ignore.case = TRUE))) %>%
  ungroup() %>%
  select(seq_id, product) %>%
  distinct()

seq_ids_without_dehydratase

# Check if ICE part of plasmids
ice_contigs %in% plasmid_df$seq_id
product_seqid_count <- plasmid_df %>%
  distinct(seq_id, product)	 %>%
  group_by(product) %>%
  summarise(seq_id_count = n())

total_seq_ids <- n_distinct(plasmid_df$seq_id)

# Calculate the percentage of seq_ids for each product and standardize product names
product_seqid_percentage <- product_seqid_count %>%
  mutate(percentage = (seq_id_count / total_seq_ids) * 100
  ) %>%
  filter(percentage > 50) %>%
  filter(product != "hypothetical protein") %>%
  arrange(desc(percentage)) 

product_seqid_percentage %>% filter(grepl("thio", product, ignore.case = TRUE))

product_seqid_count %>% 
  filter(grepl("Tet", product, ignore.case = TRUE))

product_seqid_percentage_80 <- product_seqid_count %>%
  mutate(percentage = (seq_id_count / total_seq_ids) * 100) %>%
  filter(percentage > 80) %>%
  filter(product != "hypothetical protein")


high_quality_plasmid_content_plot = ggplot(product_seqid_percentage, 
  aes(x = percentage, y = reorder(product, percentage))) +
  geom_bar(stat = "identity", color = "black", size=0.25, aes(
    fill = ifelse(product == "Tetracycline resistance protein", "#E64B35B2", 
    ifelse(product %in% nisin_products, "#3C5488B2", "gray")))) +
  theme_bw(base_size = 6, base_family = "Arial") +
  labs(x = "% of plasmids", y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.8, vjust = 0), 
    plot.title = element_text(size = 12),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "bottom") +
  scale_fill_identity()

ggsave(high_quality_plasmid_content_plot, 
  file = "../figures/high_quality_plasmid_content_plot_plot.png",
  width = 8, height = 6, units = "cm", dpi=600)  

fwrite(product_seqid_count, "../data/tables/supplementary/supplementary_table_S9.csv", sep="\t")


################################
# content of transposons
################################
tn_content <- fread("../data/tables/tn_content.tsv", sep="\t")

# Count the number of unique seq_id for each product
product_seqid_count_tn <- tn_content %>%
  distinct(nucleotide_acc, product) %>%
  group_by(product) %>%
  summarise(seq_id_count = n()) %>%
  arrange(desc(seq_id_count))

# Calculate the total number of unique seq_id
total_seq_ids_tn <- n_distinct(tn_content$assembly)
total_seq_ids_tn == 123

# Calculate the percentage of seq_id for each product
product_seqid_percentage_60_75 <- product_seqid_count_tn %>%
  mutate(percentage = (seq_id_count / total_seq_ids_tn) * 100) %>%
  filter(percentage > 60 & percentage < 100) %>%
  filter(product != "hypothetical protein") %>%
  arrange(desc(percentage))

product_seqid_count %>% 
  filter(grepl("Tet", product, ignore.case = TRUE))

product_seqid_percentage_80 <- product_seqid_count_tn %>%
  mutate(percentage = (seq_id_count / total_seq_ids_tn) * 100) %>%
  filter(percentage > 80) %>%
  filter(product != "hypothetical protein")

product_seqid_percentage <- product_seqid_count_tn %>%
  mutate(percentage = (seq_id_count / total_seq_ids_tn) * 100) %>%
  filter(product != "hypothetical protein")

################################
# plot of counts on MGE (Figure 4)
################################
mge_table_class = fread("../data/tables/supplementary/supplementary_table_S6.csv") 

n_distinct(mge_table_class$nucleotide_acc) / 915 * 100

count_mge_species = mge_table %>% select(species, nucleotide_acc, mge_class) %>%
  group_by(nucleotide_acc, species, mge_class) %>%
  distinct() %>%
  group_by(nucleotide_acc, species) %>%
  summarise(mge_class = paste(mge_class, collapse = ", "), .groups = 'drop') %>% # t
  group_by(species, mge_class) %>%
  summarise(count = n_distinct(nucleotide_acc)) %>%
  arrange(desc(count))

count_mge = mge_table %>% select(species, nucleotide_acc, mge_class) %>%
  group_by(nucleotide_acc, species, mge_class) %>%
  distinct() %>%
  group_by(nucleotide_acc, species) %>%
  summarise(mge_class = paste(mge_class, collapse = ", "), .groups = 'drop') %>% # t
  group_by(mge_class) %>%
  summarise(count = n_distinct(nucleotide_acc)) %>%
  arrange(desc(count)) %>%
  mutate(pecent = (count / 915) * 100)

mge_table_df_for_plot = mge_table_class %>%
  # left_join(., dplyr::select(mge_table, nucleotide_acc, genus), by = "nucleotide_acc") %>%
  select(nucleotide_acc, mge_class, genus) %>%
  group_by(nucleotide_acc, genus, mge_class) %>%
  distinct() %>%
  dplyr::slice(1) %>%
  group_by(nucleotide_acc, genus) %>%
  summarise(mge_class = paste(mge_class, collapse = ", "), .groups = 'drop') %>% # this stops double counting
  mutate(mge_class = case_when(
    mge_class == "ICE_element, transposon" ~ "ICE_element",
    TRUE ~ mge_class)) %>%
  group_by(genus, mge_class) %>%
  summarise(count = n_distinct(nucleotide_acc), .groups = 'drop') %>% 
  filter(count > 0) %>%
  mutate(subset = case_when(count > 50 ~ "yes", TRUE ~ "no")) %>%
  arrange(desc(count)) 


  
strep_df = mge_table_df_for_plot %>%
  filter(grepl("Streptococcus", genus, ignore.case = TRUE))

mge_table_plot = ggplot(strep_df, 
  aes(x = count, y = reorder(mge_class, -count, function(x) -sum(x)), fill = mge_class)) +
  geom_bar(stat = "identity", color="black") +
  theme_bw(base_size = 10) +
  labs(x = "g__Streptococcus MGE count", y = "", fill = "MGE Class") +  # Add fill legend name
  theme(axis.text.x = element_text(angle = 90, hjust = 0.8, vjust = 0), 
    plot.title = element_text(size = 12),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "none") +
  scale_fill_manual(values = c(pal_npg("nrc")(10) )) 

mge_table_plot2 = ggplot(mge_table_df_for_plot %>% filter(genus != "g__Streptococcus"), 
  aes(x = count, y = reorder(genus, -count, function(x) -sum(x)), fill = mge_class)) +
  geom_bar(stat = "identity", color="black") +
  theme_bw(base_size = 10) +
  labs(x = "Count", y = "Genus", fill = "MGE Class") +  # Add fill legend name
  theme(axis.text.x = element_text(angle = 90, hjust = 0.8, vjust = 0), 
    plot.title = element_text(size = 12),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "right") +
  scale_fill_manual(values = c(pal_npg("nrc")(10) )) 


ggsave(mge_table_plot, 
  file = "../figures/mge_table_plot_plot.png",
  width = 14, height = 6, units = "cm", dpi=600)

mge_table_plot3 = grid.arrange(mge_table_plot, mge_table_plot2, 
  ncol = 2,
  widths = c(4, 5))

ggsave(mge_table_plot3,
  file = "../figures/mge_table_3_plot.png", 
  width = 26, height = 8, units = "cm", dpi=600)
