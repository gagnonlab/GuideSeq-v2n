# Suppress package load messages
suppressMessages(library(tidyverse))
suppressMessages(library(BSgenome))

# Get path to input and output
args <- commandArgs(trailingOnly=TRUE)
input_file = args[2]
output_file = args[3]

# Load libraries
library(tidyverse)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38

# Load input
message(input_file)
CHANGEseq_matched = suppressMessages(read_tsv(input_file,col_names=TRUE))
names(CHANGEseq_matched)[names(CHANGEseq_matched) == "BED_Max.Position"] <- "end"
names(CHANGEseq_matched)[names(CHANGEseq_matched) == "BED_Min.Position"] <- "start"
names(CHANGEseq_matched)[names(CHANGEseq_matched) == "#BED_Chromosome"] <- "chr"
names(CHANGEseq_matched)[names(CHANGEseq_matched) == "total.sum"] <- "reads"
names(CHANGEseq_matched)[names(CHANGEseq_matched) == "Cell"] <- "sample"
names(CHANGEseq_matched)[names(CHANGEseq_matched) == "Site_SubstitutionsOnly.NumSubstitutions"] <- "distance"
#names(CHANGEseq_matched)[names(CHANGEseq_matched) == "Site_SubstitutionsOnly.Sequence"] <- "name"

# Get chromosome lengths into hg38_tbl
#   Note as.numeric() required to store numbers greater than 2E9
#   Filter out chromosome strings that contain '_'
#   Make a cumulative sum to get y-coordinate of Manhattan plot
#   Select all row except for chr_len (temporary variable)
hg38_tbl = tibble( chr = names(seqlengths(hg38)), chr_len = as.numeric(seqlengths(hg38)) ) %>%
  filter(!str_detect(chr,'_')) %>%
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len)

# Make an annotated dataset called circleseq_matched_manhattan_annotated
#   Join circleseq_matched table with hg38_tbl made above; 
#   Left join means that all rows in left table are joined with something in right
#   Arrange by sample, chromosome, start
#   Add a variable called BPcum that adds the position to the total calculated above
#   Select only sample, chr, start, tot, BPcum, and reads for final table
circleseq_matched_manhattan_annotated = CHANGEseq_matched %>% 
filter(chr!="chrM") %>%
left_join(hg38_tbl, by="chr") %>% 
mutate(sample = str_replace(sample, "CRL[0-9]{3}_", "")) %>%
arrange(sample, chr, start) %>%
mutate(BPcum=start+tot) %>%
select(sample, chr, start, tot, BPcum, reads, distance) %>%
group_by(sample) %>%
mutate(label_x=ifelse(distance==0, BPcum, NA), label_y=ifelse(distance==0, reads, NA), 
       label_distance=ifelse(distance==0, 0.1*max(reads), NA))

# Calculate positions to put labels for X-axis  
axisdf = hg38_tbl %>% 
  mutate(center=(tot + lead(tot))/2) %>%
  arrange(center)

# Calculate number of pages - Unused?
npages = ceiling(length(unique(circleseq_matched_manhattan_annotated$sample))/10)

# Generate plot
#   Puts a dark-red arrow over the on-target
p <- circleseq_matched_manhattan_annotated %>%
  ggplot(aes(x=BPcum, y=reads, color=chr)) +
  geom_segment(aes(x=BPcum, y=0, xend=BPcum, yend=reads, color=as.factor(chr)), alpha=0.9, size=0.5) +
  geom_segment(aes(x=label_x, y=reads+2*label_distance, xend=label_x, yend=reads+0.5*label_distance),
               color="darkred", linejoin="mitre", arrow=arrow(length=unit(0.2, "cm")), size=0.5) +
  scale_x_continuous(name="Chromosome", label=axisdf$chr, breaks=axisdf$center, limits=c(0, tail(hg38_tbl$tot,1))) +
  scale_y_continuous(name="GUIDE-seq V2n read count", expand=c(0,0)) +
  theme_classic(base_size=15) +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle=60, hjust=1.1, vjust=1.1, size=8))

# Save plot to output
suppressWarnings(ggsave(output_file, plot = p, width=11, height=4))

