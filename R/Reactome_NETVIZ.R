#script to create term-term netowrk for the terms in Reactome filtered
#GMT file. Terms are connected by the significance of the sharing of genes
#between them (represented by the -log10 of the Fisher's exact test pval)

library(dplyr)
library(clusterProfiler)
library(stringr)
options(stringsAsFactors = F)

#load necessary objects
GMT <- read.gmt("data/ReactomePathwaysAll.gmt")
GMT$ont <- str_remove_all(string = GMT$ont,pattern = "REACTOME_")
GMT$ont <- str_replace_all(string = GMT$ont,pattern = "_",replacement = " ")
GMT$ont <- tolower(GMT$ont)

#create term network
term_network <- GMT
colnames(term_network) <- c("Source","Target")

term_nodes <- data.frame(Id=unique(c(term_network$Source,term_network$Target)),
                         Label=unique(c(term_network$Source,term_network$Target)))
term_nodes$Class <- c(rep("Term",length(unique(term_network$Source))),
                      rep("Gene",length(unique(term_network$Target))))

write.csv(term_network,file="data/term_network_edges.csv",quote = T,row.names = F)
write.csv(term_nodes,file="data/term_network_nodes.csv",quote = T,row.names = F)

all_terms <- term_nodes %>%
  filter(Class=="Term") %>%
  pull(Label) %>%
  unique()

term_term_network <- as.data.frame(t(combn(all_terms,2)))
colnames(term_term_network) <- c("Source","Target")

for (i in 1:length(term_term_network$Source)){
  term_1 <- term_term_network$Source[i]
  term_2 <- term_term_network$Target[i]
  genes_1 <- term_network %>%
    filter(Source==term_1) %>%
    pull(Target) %>%
    unique()
  genes_2 <- term_network %>%
    filter(Source==term_2) %>%
    pull(Target) %>%
    unique()
  comm <- intersect(genes_1,genes_2)
  pval <- round(phyper(q = length(comm)-1,
                       m = length(genes_1),
                       n = length(unique(term_network$Target))-length(genes_1),
                       k = length(genes_2),
                       lower.tail = F),
                digits = 10000)
  term_term_network$pval[i] <- pval
  term_term_network$comm[i] <- length(comm)
  term_term_network$logpval[i] <- -log(pval+0.0000000000000000001)
}

term_term_nodes <- data.frame(Id=all_terms_evolution,
                              Label=all_terms_evolution)

write.csv(term_term_network,file="data/term_term_network_edges.csv",row.names = F)
write.csv(term_term_nodes,file="data/term_term_network_nodes.csv",row.names = F)
