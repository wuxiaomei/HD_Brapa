library(dplyr)

## do fisher.test and adjust p values
work_path <- "./"

##########################################
# Read in command line options
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
min_kmer <- args[2]
max_kmer <- args[3]
#file <- "reg_BraA01g023140.3.5C_targ"

# do fisher.test
print (file)
for (k in seq(min_kmer, max_kmer, 1))
#for (k in seq(6, 6, 1))
{
    print (k)
    count_table0 <- read.table(paste0(file, k), sep='\t', header=TRUE)
    count_table <- count_table0 %>% filter(Nmer_t > 4 & Nother_t > 4)
    count_table <- count_table %>% mutate(odds.ratio = NA, p.value.fish = NA, p.adjust.fish = NA)
    for (i in 1:nrow(count_table))
    {
        x<-matrix(c(
            count_table[i, 'Nmer_t'],
            count_table[i, 'Nmer_c'],
            count_table[i, 'Nother_t'],
            count_table[i, 'Nother_c']
        ), nrow=2)
        res <- fisher.test(x)
        count_table[i, 'odds.ratio'] <- res$estimate
        count_table[i, 'p.value.fish'] <- res$p.value
    }
    # adjust p-values
    p.adjusted.bh <- p.adjust(count_table$p.value.fish, method = "BH")
    count_table$p.adjust.fish <- p.adjusted.bh

    # choose motifs with adjusted p.value < 0.05 to file, 
    #  users can also reset the cutoff
    count_table %>% filter(p.adjusted.bh < 0.05) %>% 
        filter(odds.ratio > 1) %>%
        write.table(paste0(file, k, ".p"), quote=FALSE, sep="\t", row.names=FALSE)
} # k

