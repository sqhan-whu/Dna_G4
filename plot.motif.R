library(gridExtra)
library(tidyverse)
piePlot <- function(count, categories) {
    dat <- data.frame(count = count, category = categories)
    dat$fraction <- dat$count / sum(dat$count)
    dat$ymax <- cumsum(dat$fraction)
    dat$ymin <- c(0, head(dat$ymax, n = -1))
    dat$label <- factor(paste(dat$category, dat$count), levels = paste(dat$category, dat$count))
    plot <-
        ggplot(dat, aes(
            fill = label, # fill by label not category
            ymax = ymax,
            ymin = ymin,
            xmin = 0,
            xmax = 1
        )) +
        geom_rect() +
        coord_polar(theta = "y") +
        theme(legend.position="top") + theme_void() # no need for labels anymore
    plot
}

generate_structure_plots <- function(dataframe_peaks,outuput_file,outuput_file_pie,outuput_file_fe, label_to_print,labels_legend){
  
  #pie chart
  pie_all_peaks <- ggplot(dataframe_peaks, aes(x = "", y = actual_count, fill = motif_name)) +
  geom_bar(width = 1, stat = "identity", 
           color = "white") +
  coord_polar("y", start = 0) + 
  scale_fill_manual(values=c("skyblue4","tomato3","darkseagreen3","mediumpurple1","lightskyblue","coral","azure3"),
                    label=labels_legend)+
  #geom_text(aes(label = actual_count), position = position_stack(vjust = 0.5)) +
  ggtitle(label_to_print) + 
  theme(legend.position ="rigth") +
  theme_void()
  
  bar_fold_enrichments <- ggplot(dataframe_peaks, 
       aes(x = motif_name, y = fold_change_over_random, fill = motif_name)) +
  geom_bar(stat = "identity", color = "white",
           width = 0.7) + 
  scale_fill_manual(values=c("skyblue4","tomato3","darkseagreen3","mediumpurple1","lightskyblue","coral","azure3"),
                    label=label_to_print)+
  geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.8) +
  ggtitle(paste0("ESC - fold_enrich")) + ylab('fold_enrichm') + 
  theme(legend.position = "none",
        text = element_text(size=12),
        #labels=labels_legend,
        axis.text.x = element_text(angle = 90))
  
  Final_plot <- grid.arrange(pie_all_peaks,bar_fold_enrichments)
  ggsave(file=outuput_file,Final_plot)
  ggsave(file=outuput_file_pie,pie_all_peaks)
  ggsave(file=outuput_file_fe,bar_fold_enrichments)
  rm(Final_plot)
}
#setwd('/Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/PDTX_only_dm6_inputSubtraction/DG4R/hg19_old_new_q005_fasta')
assign('hg19_all_peaks_motifs', get(load('2_peak.results_motifs.Rdata')))

hg19_all_peaks_motifs <- hg19_all_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))

labels_legend_hg19_all_peaks_motifs <- paste0(hg19_all_peaks_motifs$actual_count,' (',hg19_all_peaks_motifs$percentage,'%)')


generate_structure_plots(hg19_all_peaks_motifs,
                         'hg38_2_all_peaks_motifs.pdf',
                         'hg38_2_all_peaks_motifs_pie_counts.pdf',
                         'hg38_2_all_peaks_motifs_fold_enrich.pdf',
                         'hg38_2_all_peaks',labels_legend_hg19_all_peaks_motifs)
write.table(results_motifs,file='hg38_2.sorted.results_motifs.csv',sep = ",",row.names =F, col.names = T)
