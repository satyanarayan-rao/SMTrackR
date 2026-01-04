#' Visualize Single Molecule Footprint Patterns with Ideogram and Gene Tracks
#'
#' @description Generates Gviz-compatible R code that puts SMF footprints below ideogram and genetracks. This call assumes that you have already ran `plotFootprints` function. 
#' Designed for Drosophila melanogaster analysis with expandable parameters for other organisms.
#'
#' @param organism Organism code (default: "dmelanogaster")
#' @param model Biological model/system (default: "S2" cells)
#' @param condition Experimental condition (default: "WT")
#' @param genome_assembly Genome version (default: "dm6")
#' @param type Data type ("dSMF" = dual enzyme Single Molecule Footprinting)
#' @param chromosome Chromosome ID (e.g., "chr2L")
#' @param start Genomic start position (numeric/character)
#' @param stop Genomic end position (numeric/character)
#' @param tr Track name for BigBed file resource (default: "fp_and_mvec")
#' @param label Plot title annotation
#' @param span_left Upstream window size from region (default: 150)
#' @param span_right Downstream window size from region (default: 150)
#' @param remove_dup Remove duplicate reads? (default: FALSE)
#' @param fp_cap Maximum footprint value for y-axis scaling (default: 50)
#' @param gviz_left left zoom (bases) from the start
#' @param gviz_right right zoom (bases) from the stop
#' @param target_dir destination directory (must be the same used in `plotFootprints` function)
#'
#' @return Saves a heatmap in pdf, png, and eps file formats.
#'
#' @details
#' This function generates a R code running which one can get ideogram, genetracks and SMF heatmap altogether in one plot.
#'
#' @examples
#' # Basic usage with default parameters
#' generateGvizCodeforSMF(organism = "mmusculus", model = "16cell", 
#'                          condition = "WT", genome_assembly = "mm10",
#'                          type = "SMF", chromosome = "chr1",start = 191718250,
#'                          stop = 191718280, tr = "16cell", label = "tss",
#'                          fp_cap = 50, remove_dup = F, gviz_left = 500,
#'                          gviz_right = 500, target_dir = "")
#'                          
#'
#' @export
generateGvizCodeforSMF <- function (organism = "mmusculus", model = "16cell", 
                       condition = "WT", genome_assembly = "mm10",
                       type = "SMF", chromosome = "chr1", start = 191718250,
                       stop = 191718280, tr = "16cell", label = "tss",
                       fp_cap = 50, remove_dup = F, gviz_left = 500, 
                       gviz_right = 500, target_dir = "") {
    # Function implementation
}

generateGvizCodeforSMF <- function (organism = "demlanogaster", model = "S2",
        condition = "WT", genome_assembly = "dm6",
        type = "dSMF", chromosome = "chr2L",
        start = 480290, stop = 480320,
        tr = "fp_and_mvec", label = "peak229",
        span_left = 150, span_right = 150, remove_dup = FALSE,
        fp_cap = 50, gviz_left = 1000, gviz_right = 1000, target_dir = ""){

    fp_file = paste (target_dir, paste0(label, ".num.fp.tsv"), sep = "/")
    occ_file = paste (target_dir, paste0("plots/", label, ".states.tsv"), 
                      sep = "/")

    fp_file_flag = file.exists(fp_file)
    occ_file_flag = file.exists(occ_file)
    if (fp_file_flag && occ_file_flag){

        # get number of valid states and their count
        occ_dt = read.table(occ_file, sep = "\t", header = T,
                            stringsAsFactors = F)
        num_states = dim(occ_dt)[1]
        total_molecules = sum(occ_dt$count)

        gviz_code_file = paste(target_dir, paste0(label, ".gviz.R"), sep = "/")
        libs_and_options = c("library(biomaRt)", "library(Gviz)",
                             "options(ucscChromosomeNames=TRUE)")
        cat (libs_and_options, sep = "\n", file = gviz_code_file)

        chrom = paste0("chrom = ", '"', chromosome, '"')
        genome = paste0("genome = ", '"', genome_assembly, '"')


        cat (paste0("locus = ", as.integer((start + stop)/2)),
             file = gviz_code_file, sep = "\n", append = TRUE)
        cat (paste0("start = ", start),
             file = gviz_code_file, sep = "\n", append = TRUE)
        cat (paste0("gviz_left = ", gviz_left),
             file = gviz_code_file, sep = "\n", append = TRUE)
        cat (paste0("gviz_right = ", gviz_right),
             file = gviz_code_file, sep = "\n", append = TRUE)
        cat (paste0("span_left = ", span_left),
             file = gviz_code_file, sep = "\n", append = TRUE)
        cat (paste0("span_right = ", span_right),
             file = gviz_code_file, sep = "\n", append = TRUE)
        cat ("starts <- seq (locus - span_left, locus + span_right, by = 1)",
             file = gviz_code_file, sep = "\n", append = TRUE)
        cat ("ends <- starts + 1",
             file = gviz_code_file, sep = "\n", append = TRUE)
        cat (chrom, genome, file = gviz_code_file, sep = "\n", append = TRUE)
        ideogram = paste0("itrack = IdeogramTrack(genome = ", '"',
                          genome_assembly,'"',
                          ", chromosome = ", '"',  chromosome, '"',")" )
        cat (ideogram, file = gviz_code_file, sep = "\n", append = TRUE)

        mart_str = paste0 ("mart <- useMart('ensembl', dataset = '",
        organism, "_gene_ensembl'",
        ", host = 'https://nov2020.archive.ensembl.org/')")

        cat (mart_str, file = gviz_code_file, sep = "\n", append = TRUE)

        cat (paste0("genes <- BiomartGeneRegionTrack(
        genome = genome,
        biomart = mart,
        chromosome = chrom,
        start = start - ", gviz_left,
        ", end = start + ", gviz_right,
        ", name = 'Gene Model')"),
             file = gviz_code_file, sep = "\n", append = TRUE)
        cat ("displayPars(genes) = list(transcriptAnnotation='symbol',
    showExonId=FALSE)",
             file = gviz_code_file, sep = "\n", append = TRUE)
        fp_file = paste0(label, ".num.fp.tsv")
        cat (paste0("real_data = read.table('", fp_file, "'",
                    ", sep = '\\t', header = F , stringsAsFactors = F,
                comment.char = '%')"),
             file = gviz_code_file, sep = "\n", append = TRUE)

        cat (paste("row.names(real_data) = real_data$V1",
                   "real_data$V1 = NULL",
                   "real_data = as.matrix(real_data)",
                   "df_to_list = list()", sep = "\n"),
             file = gviz_code_file, sep = "\n", append = TRUE)
        prev = 1
        cat ("axisTrack <- GenomeAxisTrack()",
             file = gviz_code_file, sep = "\n", append = TRUE)

        cat ("bstarts <- seq (locus - span_left - 20, locus + span_right, by = 1)
        bends <- bstarts + 1
        gr = GRanges(seqnames = chrom,
                     ranges = IRanges(start = bstarts, end = bends),
                     score = rep (-1, 321))
        red_line <- DataTrack(
            range = gr,
            type = 'heatmap',        # heatmap mode
            name = '',
            gradient = c('red', 'white', 'blue'),  # color scale
            ylim = c(-1,0),
            showAxis = FALSE,
            background.title = 'white')",
             file = gviz_code_file, sep = "\n", append = TRUE)

        # set track sizes :
        # ideogram : fixed,10% of the plot height
        # axisTrack: fixed, 5% of the plot height
        # geneTrack: fixed, 25% of the plot height
        # rest is SMF plots
        all_tracks = "c(axisTrack, genes"
        track_sizes = "c(0.1, 0.05, 0.25"
        line_counter = 0
        heatmap_height = 0.4 # 40% of the plot height
        per_line_width = round (heatmap_height/total_molecules, 6)
        for (i in seq(num_states)){
            state_label = occ_dt[i, "state_label"]
            mol_cnt = occ_dt[i, "count"]
            current = prev + mol_cnt - 1

            cat (paste0("df_to_list[[", "'", state_label,"'", "]] = ",
            "lapply(seq(", prev, ",", current, "), ", "function(x) {
        gr = GRanges(seqnames = chrom,
                     ranges = IRanges(start = starts, end = ends),
                     score = unname(real_data[x,]))
        heatmapTrack <- DataTrack(
            range = gr,
            type = 'heatmap',
            name = '',
            #gradient = c('grey', 'white', '#41b6c4'),  # color scale
            gradient = c('grey', 'white', '#1f78b4'),  # color scale
            ylim = c(-1, 1),
            showAxis = FALSE,
            background.title = 'white')
        return (heatmapTrack)})"),
                 file = gviz_code_file, sep = "\n", append = TRUE)
            prev = current + 1

            track_sizes = paste0(track_sizes, ", rep (",
                                 per_line_width, ", ", mol_cnt, ")")
            all_tracks = paste0(all_tracks, ",
                                df_to_list[['", state_label, "']]")
            if (line_counter < num_states - 1){ # if 3 states then counter will run for 0 and 1
                track_sizes = paste0(track_sizes, ", ", per_line_width/2)
                line_counter = line_counter + 1
                all_tracks = paste0(all_tracks, ", red_line")
            }
        }

        all_tracks = paste0(all_tracks, ")")
        track_sizes = paste0(track_sizes, ")")
        # highLightTracks

        cat (paste0("hl_tracks <- HighlightTrack (trackList = ", all_tracks,
                    ", genome = '", genome_assembly, "'",
                    ", chromosome = '", chromosome, "'",
                    ", col = 'red', fill = 'transparent', start = c (",
                    start, ")", ", width = ", stop - start, ", name = 'RoI' )")
             , file = gviz_code_file, sep = "\n", append = TRUE)
        cat ("\n", file = gviz_code_file, sep = " ", append = TRUE)

        cat ("displayPars(hl_tracks) = list (inBackground = FALSE, col = 'red')",
             file = gviz_code_file, sep = "\n", append = TRUE)
        cat ("\n", file = gviz_code_file, sep = " ", append = TRUE)
        cat (paste0("pdf('plots/", label, ".gviz.pdf'"), ",", " height = 10, width = 8)",
             file = gviz_code_file, sep = " ", append = TRUE)
        cat ("\n", file = gviz_code_file, sep = "\n", append = TRUE)

        cat (paste0("plotTracks(c(itrack, hl_tracks)",
            ", chromosome = chrom,
            from = locus - gviz_left,  # zoom start
            to = locus + gviz_right,    # zoom end
            background.title = 'white',
            showColorBar = F,
            showAxis = F, showLabel = F,
            sizes =", track_sizes, ")"), file = gviz_code_file,
             sep = "\n", append = TRUE)

        cat ("dev.off()",
             file = gviz_code_file, sep = "\n", append = TRUE)
        cat (paste0("A R code file ", gviz_code_file, 
                    " is written. Please run the command: Rscript ",
                    gviz_code_file, " to place the heatmap on gene track!\n"))

    }else{
        cat (paste0("Both files", fp_file, " ", occ_file, " ", 
                    "must exist to generate the code\n"))
        return ("Please call plotFootprints function before you run this function")
    }

}

