########################################
## Packages - don't package shame people
########################################

.libPaths()
.libPaths( c( "/data/san/data0/users/david/rstudio/packages" , .libPaths() ) )
newlib <- "/data/san/data0/users/david/rstudio/packages"
.libPaths()
library(data.table)
library(ape)
library(dplyr)
library(gggenomes)
library(thacklr)
library(tools)
library(msa)
library(Biostrings)
library(tinytex)
library(remotes)
library(data.table)
library(dplyr)
library(stringr)
library(ggtree)
library(ggsci)

texlib <- "/data/san/data0/users/david/rstudio/tinytex"
setwd("/data/san/data2/users/david/nisin/code")


##############################
#### nisin phylogenetic tee
##############################
core <- readAAStringSet("../data/nisin_vp/nisin_corepeptides.faa", format = "fasta")
core_alignment = msaMuscle(core, order = "aligned")
msaPrettyPrint(core_alignment, file="../figures/core_peptides_of_known.tex", output="tex",
               showNames="left", showNumbering="none", showLogo="top",
               showConsensus="bottom",
               consensusThreshold= c(90,100),
               shadingMode = "identical",
               logoColors="rasmol",
               verbose=FALSE, askForOverwrite=FALSE, paperWidth = 7, paperHeight = 6, 
               alFile =  "../data/nisin_vp/core_peptides_of_known_aln.fasta")
    leader_texfile <- "../figures/core_peptides_of_known.tex"
tinytex::pdflatex(leader_texfile)

### RAXML was run on the alignment file
tree_core <- read.tree("../data/nisin_vp/core_peptides_of_known_aln.fasta.raxml.supportFBP")
tree_core <- ape::root.phylo(tree_core, 
     outgroup = "kunkecin-A")

gg_core <- ggtree(tree_core, size = 0.8) +
  geom_tiplab(size=5, align=FALSE, linesize=.5, offset = 0.02, aes(color=(label == "nisin-VP")), fontface="bold", fill="grey") +
  scale_color_manual(values=c("TRUE"="purple", "FALSE"="black"), guide=FALSE) + # Custom color scale, no legend
  theme_tree2() +
  ggtitle("") +
  coord_cartesian(clip="off") +
  expand_limits(x = 1.7, y = 0.1)

ggsave(gg_core, 
	file ="../figures/nisin_core_tree.png",
	height = 5, width = 5, dpi = 300)



##############################
# SYNTENY of operons
##############################	
ava = read_sublinks("../data/nisin_vp/synteny/region_proteins.blast", format = "blast") 


synteny_table = fread("../data/nisin_vp/synteny/vp_synteny_table.tsv") 




nisin_operons_plot = gggenomes(synteny_table) %>%
  add_sublinks(ava %>% filter(length > 200)) %>%
  pick_seqs(c(
    "Blautia producta blauticin",
    # "Blautia Obeum nisin-O1-4",
    "Velocimicrobium porci nisin-VP",
    "Streptococcus agalactiae nisin-P",
    "Lactococcus lactis subsp. lactis nisin-Z",
    "Lactococcus lactis subsp. lactis nisin-A",
    "Lactococcus lactis nisin-Q",
    "Streptococcus hyointestinalis nisin-H",
    "Streptococcus salivarius nisin-G",
    "Streptococcus uberis nisin-U",
    "Apilactobacillus kunkeei kunkecin-A"  
    )) +
  geom_seq() +
  geom_link(aes(fill=pident, color=pident), offset = 0.05) +
  scale_color_gradient(low="grey88", high="grey50") +
  scale_fill_gradient(low="grey88", high="grey50") +
  ggnewscale::new_scale("fill") +
  geom_gene(aes(fill = gene), size = 6, show.legend = NA) +
    scale_fill_manual("gene",
      values = c("grey88",pal_npg("nrc", alpha=1)(10)), 
      na.value = "grey88", 
      guide = "legend") +
  coord_cartesian(clip = 'off') +
  geom_bin_label(
    size = 4, 
    # fontface = "italic", 
    family = "Arial", 
    expand_left = 1) +
  theme(legend.position="right") +
  theme_gggenomes_clean(base_size = 14, base_family = "Arial") 

nisin_operons_plot_flip = nisin_operons_plot  %>% 
  flip_seqs((c("Velocimicrobium porci nisin-VP",unique_seq_ids_vector[c(8, 5,10)])))

ggsave(nisin_operons_plot_flip, 
  file="../figures/nisin_operons.png", 
  width = 22, 
  height = 16, 
  dpi=600,
  units="cm",
  limitsize = FALSE)
