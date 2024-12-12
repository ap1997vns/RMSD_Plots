library(ggplot2)
library(gridExtra)
library(plyr)
library(grid)

load_and_filter_data <- function(file_path, selected_proteins) {
  data <- read.table(file_path, header = FALSE, sep = ',')
  colnames(data) <- c("protein", "domain", "time", "val")
  data <- data[data$protein %in% selected_proteins & data$val <= 11, ]
  data <- data[data$domain == "Protein", ]  
  data <- data[c(TRUE, rep(FALSE, 19)), ] 
  return(data)
}

transform_data <- function(data) {
  transform(data, domain = factor(domain, 
                                  levels = c("Protein"), 
                                  labels = c("Full length Protein")))
}

data_files <- list(
  list(file = "./newT180I.txt", proteins = c("APOT180I", "LFT180I", "LBT180I", "APOEFT180I")),
  list(file = "./newS171L1.txt", proteins = c("APOS171L", "LFS171L", "LBS171L", "APOEFS171L")),
  list(file = "./newA97V.txt", proteins = c("APOA97V", "LFA97V", "LBA97V", "APOEFA97V")),
  list(file = "./newM188I.txt", proteins = c("APOM188I", "LFM188I", "LBM188I", "APOEFM188I")),
  list(file = "./newWT.txt", proteins = c("APOWT", "LFWT", "LBWT", "APOEFWT")),
  list(file = "./finalR102Q.txt", proteins = c("APOR102Q", "LFR102Q", "APOEFR102Q", "LBR102Q")),
  list(file = "./newR71C.txt", proteins = c("APOR71C", "LFR71C", "LBR71C", "APOEFR71C")),
  list(file = "./newE152A.txt", proteins = c("APOE152A", "LFE152A", "LBE152A", "APOEFE152A"))
)

all_data <- ldply(data_files, function(df) transform_data(load_and_filter_data(df$file, df$proteins)))
y_limits <- range(all_data$val)

create_plot <- function(data, protein_labels, color_values, y_limits) {
  ggplot(data, aes(x = time, y = val)) +
    geom_line(aes(colour = protein, linetype = protein), linewidth = 1.2) +
    ylab(expression(RMSD(ring(A)))) +
    xlab("Time (ns)") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 1.0, linetype = 1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 45, face = "bold")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = y_limits) +
    scale_linetype_manual(values = rep(c("solid", "dashed", "dotted"), length.out = length(protein_labels)),
                          name = "protein",
                          limits = protein_labels) +
    scale_colour_manual(values = color_values,
                        name = "protein",
                        limits = protein_labels)
}

plot_params <- list(
  list(data = subset(all_data, protein %in% c("APOWT", "LFWT", "LBWT", "APOEFWT")), protein_labels = c("APOWT", "LFWT", "LBWT", "APOEFWT"), colors = c("#FF9933", "#00CCCC", "#FF3333", "#3333FF")),
  list(data = subset(all_data, protein %in% c("APOA97V", "LFA97V", "LBA97V", "APOEFA97V")), protein_labels = c("APOA97V", "LFA97V", "LBA97V", "APOEFA97V"), colors = c("#FF9933", "#00CCCC", "#FF3333", "#3333FF")),
  list(data = subset(all_data, protein %in% c("APOS171L", "LFS171L", "LBS171L", "APOEFS171L")), protein_labels = c("APOS171L", "LFS171L", "LBS171L", "APOEFS171L"), colors = c("#FF9933", "#00CCCC", "#FF3333", "#3333FF")),
  list(data = subset(all_data, protein %in% c("APOM188I", "LFM188I", "LBM188I", "APOEFM188I")), protein_labels = c("APOM188I", "LFM188I", "LBM188I", "APOEFM188I"), colors = c("#FF9933", "#00CCCC", "#FF3333", "#3333FF")),
  list(data = subset(all_data, protein %in% c("APOR102Q", "LFR102Q", "LBR102Q", "APOEFR102Q")), protein_labels = c("APOR102Q", "LFR102Q", "LBR102Q", "APOEFR102Q"), colors = c("#FF9933", "#00CCCC", "#FF3333", "#3333FF")),
  list(data = subset(all_data, protein %in% c("APOT180I", "LFT180I", "LBT180I", "APOEFT180I")), protein_labels = c("APOT180I", "LFT180I", "LBT180I", "APOEFT180I"), colors = c("#FF9933", "#00CCCC", "#FF3333", "#3333FF")),
  list(data = subset(all_data, protein %in% c("APOR71C", "LFR71C", "LBR71C", "APOEFR71C")), protein_labels = c("APOR71C", "LFR71C", "LBR71C", "APOEFR71C"), colors = c("#FF9933", "#00CCCC", "#FF3333", "#3333FF")),
  list(data = subset(all_data, protein %in% c("APOE152A", "LFE152A", "LBE152A", "APOEFE152A")), protein_labels = c("APOE152A", "LFE152A", "LBE152A", "APOEFE152A"), colors = c("#FF9933", "#00CCCC", "#FF3333", "#3333FF"))
)

plots <- lapply(plot_params, function(params) create_plot(params$data, params$protein_labels, params$colors, y_limits))

ppi <- 300

save_combined_plot <- function(filename, plot_list, titles, nrow, width, height) {
  for (i in seq_along(plot_list)) {
    plot_list[[i]] <- plot_list[[i]] + ggtitle(titles[i])
  }
  png(filename, width = width * ppi, height = height * ppi, res = ppi)
  do.call(grid.arrange, c(plot_list, nrow = nrow))
  dev.off()
}

save_combined_plot("Combined_Plots.png", list(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]]), 
                   c("WT", "A97V", "S171L", "M188I", "R102Q", "T180I", "R71C", "E152A"), 
                   nrow = 3, width = 60, height = 27)
