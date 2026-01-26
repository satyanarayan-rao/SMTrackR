savePlotNanopore <- function (label = "smac_seq", span_left = 1000,
                              span_right = 1000, stride = 5, 
                              plot_title = "scerevisiae BY4741 strain",
                              x_label = "smac_seq (sacCer3 chrIII:114300-114600)",
                              target_dir = ""){
    data_to_plot = read.table (paste(target_dir, 
                                     paste0(label, ".nanopore_methylation.tsv"),
                                     sep = "/"),
                               sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                               row.names = 1)
    strand_counts_df = read.table(paste(target_dir, 
                                        paste0(label, ".strand_wise_count.tsv"),
                                        sep = "/"),
                                  sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    row.names(strand_counts_df) = strand_counts_df$V1
    mid_point_minus_strand = strand_counts_df["-", "V2"]/2
    mid_point_plus_strand = strand_counts_df["-", "V2"] + strand_counts_df["+", "V2"]/2
    
    x_width = dim(data_to_plot)[2]
    total_molecules = dim(data_to_plot)[1]
    main_dir <- target_dir
    if (!dir.exists(main_dir)){
        dir.create(main_dir, recursive = TRUE)
    }
    plot_dir <- paste(main_dir, "plots", sep = "/")
    if (!dir.exists(plot_dir)){
        dir.create(plot_dir, recursive = TRUE)
    }
    
    if (nrow(data_to_plot) > 15){
        
        pdf (file.path(plot_dir, paste0(label, ".heatmap.pdf")),
             height = 5, width = 10, bg = "white")
        dev.control("enable")
        
        par(mgp=c(1.5,0.25,0), cex = 1)
        image(1:ncol(data_to_plot), 1:nrow(data_to_plot), t(data_to_plot),
              axes = FALSE, useRaster = TRUE,
              oldstyle = FALSE, col = c("-1" = "#bdbdbd", "0" = "#FFFFFF",
                                        "1"  = "Red" , "2" = "#C7E9C0"),
              xlim = c(-30, x_width), xlab = paste0("           ", x_label),
              ylab = "", main = paste0("      ",
                                       plot_title))
        
        
        # Draw vertical red line +/- 15 bases from `0`
        left_line = span_left/stride
        right_line = x_width - (span_right/stride)
        segments(left_line , 0, left_line,
                 (total_molecules + 0.5), lwd = 0.5, col = "Blue")
        segments(right_line, 0, right_line,
                 (total_molecules + 0.5), lwd = 0.5, col = "Blue")
        segments(-30, strand_counts_df["-", "V2"],
                 x_width, strand_counts_df["-", "V2"], lwd = 0.5, col = "Blue")
        
        ticks_at <- seq(0, x_width, 100/stride)
        lab_at_ticks <- as.character((ticks_at - as.integer(x_width/2) )*stride)
        axis(1, ticks_at,
             lab_at_ticks,
             las = 2, cex.lab = 0.5, tck = -0.01, tick = FALSE)
        
        box(lwd=0.5)
        
        text (-15, mid_point_minus_strand, paste0(" - (n = ",
                                                  strand_counts_df["-", "V2"],
                                                  " )"), srt = 90, cex = 1)
        text (-15, mid_point_plus_strand, paste0(" + (n = ",
                                                 strand_counts_df["+", "V2"],
                                                 " )"), srt = 90, cex = 1)
        #saveToPDF(file.path(plot_dir, paste0(label, ".heatmap.pdf")),
        #          height = 5, width = 10)
        #saveToPNG(file.path(plot_dir, paste0(label, ".plot.png")),
        #          height = 5, width = 10, units = "in", res = 300)
        #saveToEPS(file.path(plot_dir, paste0(label, ".heatmap.eps")),
        #          height = 5, width = 10,
        #          horizontal = FALSE, paper = "special")
        dev.copy(postscript, file.path(plot_dir, paste0(label, ".heatmap.eps")),
                 height = 5, width = 10, 
                 horizontal = FALSE, paper = "special")
        dev.control("enable")
        
        dev.copy(png, file.path(plot_dir, paste0(label, ".plot.png")),
                 height = 5, width = 10, units = "in", res = 300)
        dev.off()
        dev.off()
        dev.off()
        cat (paste0("Heatmap is saved in the file ", 
                    file.path(plot_dir, paste0(label, ".heatmap.pdf"))))
        cat ("\n")
        
    }else {
        
        pdf (file.path(plot_dir, paste0(label, ".heatmap.pdf")),
             height = 6, width = 4.5, bg = "white")
        dev.control("enable")
        par(mgp=c(1.5,0.25,0), cex = 0.75)
        plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n',
             xaxt = 'n', yaxt = 'n');
        text(x = 0.5, y = 0.5, paste0("Only ", nrow (data_to_plot),
                                      " molecules, not enough for plotting"))
        
        dev.copy(postscript, file.path(plot_dir, paste0(label, ".heatmap.eps")),
                 height = 6, width = 4.5, 
                 horizontal = FALSE, paper = "special")
        dev.control("enable")
        
        dev.copy(png, file.path(plot_dir, paste0(label, ".plot.png")),
                 height = 6, width = 4.5, units = "in", res = 300)
        dev.off()
        dev.off()
        dev.off()
        cat (paste0("Heatmap is saved in the file ", 
                    file.path(plot_dir, paste0(label, ".heatmap.pdf"))))
        cat ("\n")
        
    }
    
}

