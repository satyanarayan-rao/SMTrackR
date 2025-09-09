savePlotSMForDSMF <- function (label = "peak229",
                           span_left = 150, span_right = 150, 
                           plot_title = "dmelanogaster S2 WT",
                           x_label = "peak229 (dm6 chr2L:480290-480320)"){
    dat_for_plot <- read.table(paste0(label, ".num.fp.tsv"),
                               sep = "\t", header = FALSE,
                               stringsAsFactors = FALSE,
                               row.names = 1, comment.char = "%")
    dat_for_plot <- dat_for_plot[!grepl("#3", row.names (dat_for_plot)), ]
    #reverse the dataframe: bring bottom row to the top
    #remove not valid states i.e. state 3 for now.

    dat_for_plot <- apply(dat_for_plot, 2, rev)

    main_dir <- file.path("./", "plots/")
    if (!dir.exists(main_dir)){
        dir.create(main_dir, recursive = TRUE)
    }
    plot_dir <- file.path("plots/")
    if (!dir.exists(plot_dir)){
        dir.create(plot_dir, recursive = TRUE)
    }


    all_valid_states <- data.frame(table(
        unlist(lapply (row.names(dat_for_plot),
                       function (x) {
                           as.integer(unlist(strsplit(x, split= "#"))[-1]) }))))
    names(all_valid_states) <- c("state", "count")
    state_info <- c("2" = "N", "1" = "T", "0" = "D")
    all_valid_states$state_label <- unlist(
        lapply (all_valid_states$state,
                function (x){state_info[as.character(x)]}))
    all_valid_states$line_loc <- rev(cumsum(rev(all_valid_states$count)))
    all_valid_states$prev_point <- c(all_valid_states$line_loc[-1], 0)
    all_valid_states$label_loc <- (all_valid_states$line_loc +
                                       all_valid_states$prev_point)/2
    all_valid_states$percentage <- all_valid_states$count /
        sum (all_valid_states$count)
    all_valid_states$total <- sum (all_valid_states$count)
    all_valid_states$annotation <- label
    write.table(all_valid_states,
                file.path(plot_dir, paste0(label, ".states.tsv")),
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    total_molecules <- nrow(dat_for_plot)

    jj <- as.matrix (dat_for_plot)
    x_width <- span_left + span_right + 1
    if (nrow(jj) > 15){

        pdf (file.path(plot_dir, paste0(label, ".plot.pdf")),
             height = 6, width = 4.5, bg = "white")
        dev.control("enable")
        par(mgp=c(1.5,0.25,0), cex = 0.75)
        image(1:ncol(jj), 1:nrow(jj), t(jj),  axes = FALSE, useRaster = TRUE,
              oldstyle = FALSE, col = c("-1" = "#bdbdbd", "0" = "#FFFFFF",
                                        "1"  = "Red" , "2" = "#2b8cbe"),
              xlim = c(-30, x_width), xlab = paste0("           ", x_label), 
              ylab = "", main = paste0("      ", plot_title))


        counter <- 1
        for (i in rev(seq (nrow(all_valid_states)) )){
            line_loc <- all_valid_states[i, "line_loc"]
            label_loc <- all_valid_states[i, "label_loc"]
            state_label <- all_valid_states[i, "state_label"]
            cnt <- all_valid_states[i, "count"]

            if (counter < nrow(all_valid_states)) {
                segments(-30, (line_loc + 0.5), x_width,
                         (line_loc + 0.5), lwd = 0.5)
            }
            if (counter == 1) {
                text(-15, label_loc, paste0 ("State ",  state_label,
                                             " ( n = ", cnt, ")"), srt = 90)
            }else{
                text(-15, label_loc, paste0 (state_label,
                                             " (", cnt, ")"), srt = 90)
            }
            counter <- counter + 1
        }



        # Draw vertical red line +/- 15 bases from `0`
        segments((span_left - 15) , 0, (span_left - 15),
                 (total_molecules + 0.5), lwd = 0.5, col = "Red")
        segments((span_left + 15), 0, (span_left + 15),
                 (total_molecules + 0.5), lwd = 0.5, col = "Red")
        ticks_at <- seq(0, span_left + span_right, 50)
        lab_at_ticks <- as.character(ticks_at - as.integer(
            (span_left + span_right)/2))
        axis(1, ticks_at,
             lab_at_ticks,
             srt = 2, cex.lab = 0.5, tck = -0.01, tick = FALSE)

        box(lwd=0.5)

        #saveToPDF(file.path(plot_dir, paste0(label, ".plot.pdf")),
        #          height = 6, width = 4.5)
        #saveToPNG(file.path(plot_dir, paste0(label, ".plot.png")),
        #          height = 6, width = 4.5, units = "in", res = 300)
        #saveToEPS(file.path(plot_dir, paste0(label, ".plot.eps")),
        #          height = 6, width = 4.5,
        #          horizontal = FALSE, paper = "special")
        dev.copy(postscript, file.path(plot_dir, paste0(label, ".plot.eps")),
                 height = 6, width = 4.5, 
                 horizontal = FALSE, paper = "special")
        dev.control("enable")

        dev.copy(png, file.path(plot_dir, paste0(label, ".plot.png")),
                 height = 6, width = 4.5, units = "in", res = 300)
        dev.off()
        dev.off()
        dev.off()
        
    }else {

        pdf (file.path(plot_dir, paste0(label, ".plot.pdf")),
             height = 6, width = 4.5, bg = "white")
        dev.control("enable")
        par(mgp=c(1.5,0.25,0), cex = 0.75)
        plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n',
             xaxt = 'n', yaxt = 'n');
        text(x = 0.5, y = 0.5, paste0("Only ", nrow (jj),
                                      " molecules, not enough for plotting"))

        dev.copy(postscript, file.path(plot_dir, paste0(label, ".plot.eps")),
                 height = 6, width = 4.5, 
                 horizontal = FALSE, paper = "special")
        dev.control("enable")

        dev.copy(png, file.path(plot_dir, paste0(label, ".plot.png")),
                 height = 6, width = 4.5, units = "in", res = 300)
        dev.off()
        dev.off()
        dev.off()
        
    }

}
