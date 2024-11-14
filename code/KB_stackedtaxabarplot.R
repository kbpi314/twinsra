# stacked taxabarplot from AC
library(tidyverse)
library(ggplot2)

#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")

library(qiime2R)

taxa_barplot<-function(features, metadata, category, normalize, ntoplot, sort_by="none"){
  q2r_palette<-c(
    "blue4",
    "olivedrab",
    "firebrick",
    "gold",
    "darkorchid",
    "steelblue2",
    "chartreuse1",
    "aquamarine",
    "yellow3",
    "coral",
    "grey"
  )
  
  
  if(missing(ntoplot) & nrow(features)>10){ntoplot=10} else if (missing(ntoplot)){ntoplot=nrow(features)}
  if(missing(normalize)){normalize<-"percent"}
  if(normalize=="percent"){features<-make_percent(features)} else if(normalize=="proportion"){features<-make_proportion(features)}
  if(missing(metadata)){metadata<-data.frame(SampleID=colnames(features))}
  if(!"SampleID" %in% colnames(metadata)){metadata <- metadata %>% rownames_to_column("SampleID")}
  if(!missing(category)){
    if(!category %in% colnames(metadata)){message(stop(category, " not found as column in metdata"))}
  }
  
  plotfeats<-names(sort(rowMeans(features), decreasing = TRUE)[1:ntoplot]) # extract the top N most abundant features on average
  
  suppressMessages(
    suppressWarnings(
      fplot<-
        features %>%
        as.data.frame() %>%
        rownames_to_column(var="Taxon") %>%
        gather(-Taxon, key="SampleID", value="Abundance") %>%
        mutate(Taxon=if_else(Taxon %in% plotfeats, Taxon, "Remainder")) %>%
        group_by(Taxon, SampleID) %>%
        summarize(Abundance=sum(Abundance)) %>%
        ungroup() %>%
        mutate(Taxon=factor(Taxon, levels=rev(c(plotfeats, "Remainder")))) %>%
        left_join(metadata)
    ))
  bplot<-
    ggplot(fplot, aes(x=reorder(SampleID, ave(setNames(Abundance, Taxon), SampleID, FUN = function(x) replace(x[sort_by], is.na(x[sort_by]), 0))), y=Abundance, fill=Taxon)) +
    geom_bar(stat="identity") +
    theme_q2r() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    coord_cartesian(expand=FALSE) +
    xlab("Sample") +
    ylab("Abundance")
  if(ntoplot<=10){bplot<-bplot+scale_fill_manual(values=rev(q2r_palette), name="")}
  if(!missing(category)){bplot<-bplot + facet_grid(~get(category), scales="free_x", space="free")}#, labeller = as_labeller(wrap_text))}
  return(bplot)
}

# read inputs
asv_table = read_qza(file='/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs00/metaphlan_taxa_table_abs.qza')$data
taxon_table = read_qza(file='/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/outputs/jobs03/taxonomy.qza')$data %>% parse_taxonomy() 
sum_taxa = summarize_taxa(features=asv_table, taxonomy=taxon_table)
sp_taxasums = sum_taxa$Species

# if you don't want to use qza stuff you can make your own summarized stacked barplot dataframe
sp_taxasums = read.csv(file='/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/df_q2R.tsv', 
                         sep='\t',
                         row.names='Species')

metadata=read_delim(file='/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/qiime_mapping_file_noctrl_no182183_no516517_ensemble.tsv',
                    delim='\t')#,row
metadata=metadata[-1,]

p <- taxa_barplot(sp_taxasums, 
                  metadata, 
                  category="Diagnosis", 
                  ntoplot=10, 
                  #sort_by="d__Bacteria; Actinobacteria; Actinobacteria; Bifidobacteriales; Bifidobacteriaceae; Bifidobacterium; Bifidobacterium_longum")
                  sort_by="Remainder",
                  normalize='none') + 
                ylab('Relative abundance (%)') + 
                #ggtitle('Taxa Bar Plot') +
                guides(fill = guide_legend(ncol = 1))# +
                #scale_fill_manual(values = col1) +
                #bkg
p

