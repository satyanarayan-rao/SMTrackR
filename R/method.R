fp_info_for_complete_molecule  <- function (
        footprint_string = "F...FFFF...FFF.....FFFFF."){
    # use stringr library
    first_char <- substr(footprint_string,1,1)
    footprint_loci <- data.frame(str_locate_all(footprint_string, "\\.F"))
    st <- footprint_loci$start
    en <- footprint_loci$end
    st <- st + 1
    if (first_char == "F"){
        st <- c (1, st)
    }
    
    #print (st)
    total_length <- str_length(footprint_string)
    fp_len_vec <- rep(0, total_length)
    abs_fp_start_vec <- rep(0, total_length)
    just_length_vec <- c()
    for (k in st){
        cnt <- 0
        bp <- k
        for (i in seq (k,total_length)){
            char_i <- substr(footprint_string, i,i)
            if (char_i == "."){
                break
            }else {
                cnt <- cnt + 1
                abs_fp_start_vec [i] <- k
                bp <- i
            }
        }
        
        fp_len_vec[k:bp] <- cnt
        just_length_vec <- c(just_length_vec, cnt)
        abs_fp_start_vec[k:bp] <- k
    }
    return (list(fp_len_at_each_pos = fp_len_vec,
                 just_length_vec = just_length_vec,
                 fp_start_locus = st,
                 abs_fp_start_vec = abs_fp_start_vec))
}

locate_fp_start  <- function (footprint_string = "F...FFFF...FFF.....FFFFF."){
    # use stringr library
    first_char <- substr(footprint_string,1,1)
    footprint_loci <- data.frame(str_locate_all(footprint_string, "\\.F"))
    st <- footprint_loci$start
    en <- footprint_loci$end
    st <- st + 1
    if (first_char == "F"){
        st <- c (1, st)
    }
    return (st)
}


fp_info_for_sub_molecule <- function (
        query_footprint = "FF...FFF.....FF",
        query_start = 7,
        query_end = 21,
        complete_footprint = "F...FFFF...FFF.....FFFFF."){
    fp_start_loci_on_query <- locate_fp_start(
        footprint_string = query_footprint)
    fp_info_complete_fp <- fp_info_for_complete_molecule(complete_footprint)
    fp_len_at_each_pos <- fp_info_complete_fp$fp_len_at_each_pos
    abs_fp_start_at_each_pos <- fp_info_complete_fp$abs_fp_start_vec
    abs_fp_start_in_query <- abs_fp_start_at_each_pos[query_start:query_end]
    fp_len_query <- c()
    for (i in fp_start_loci_on_query){
        idx_on_long_molecule <- i + query_start - 1
        fp_len_query <- c (fp_len_query,
                           fp_len_at_each_pos[idx_on_long_molecule])
    }
    return (list (fp_len_query = fp_len_query,
                  fp_start_loci_on_query = fp_start_loci_on_query,
                  abs_fp_start_in_query = abs_fp_start_in_query))
}

extend_fp_from_center <- function (
        sub_moleculue = "FF...FFF.....FF",
        complete_molecule = "F...FFFF...FFF.....FFFFF.",
        abs_roi_start = 7, abs_roi_end = 21, abs_complete_molecule_start = 1,
        abs_complete_molecule_end = 25, lextend = 30,
        rextend = 30, strand = "+") {
    rel_roi_start <- abs_roi_start - abs_complete_molecule_start + 1
    rel_roi_end <- abs_roi_end - abs_complete_molecule_start + 1
    peak_center <- as.integer((rel_roi_start + rel_roi_end)/2)
    len_complete_molecule <- str_length(complete_molecule)
    m_to_add_in_left <- abs (min ( (peak_center - 1) - lextend, 0))
    m_to_add_in_right <- abs ( min (
        len_complete_molecule - (peak_center + rextend), 0 ) )
    extended_sub_molecule_left <- NULL
    extended_sub_molecule_right <- NULL
    if (m_to_add_in_left > 0){
        extended_sub_molecule_left <- 1
    }else {
        extended_sub_molecule_left <- peak_center - lextend
    }
    if (m_to_add_in_right > 0){
        extended_sub_molecule_right <- len_complete_molecule
    }else{
        extended_sub_molecule_right <- peak_center + rextend
    }
    extended_sub_molecule <- str_sub(complete_molecule,
                                     extended_sub_molecule_left,
                                     extended_sub_molecule_right)
    
    final_string <- paste0(paste0(rep("M", m_to_add_in_left), collapse = ""),
                           extended_sub_molecule,
                           paste0(rep("M", m_to_add_in_right), collapse = "" ),
                           collapse = "")
}
expand_cigar_string <- function (cigar_string = ".131F52.77F35.5"){
    idx_vec <- c()
    char_vec <- c()
    for(i in seq (1, nchar(cigar_string))) {
        char <- substr(cigar_string,i,i)
        if (!(char %in% as.integer(0:9))){
            idx_vec <- c(idx_vec, i)
            char_vec <- c(char_vec, char)
        }
    }
    num_rep <- c()
    rept_vec <- c ()
    length_cigar_string <- str_length(cigar_string)
    if (length(idx_vec) == 1) {
        if (length_cigar_string > 1){
            rept <- as.integer(substr(cigar_string, 2, length_cigar_string))
            expand_cigar <- paste0 (rep(char_vec[1], rept), collapse = "")
            return (expand_cigar)
        }else { return (cigar_string) }
    }else {
        for (j in 1:(length(idx_vec) - 1) ){
            idx_c <- idx_vec[j]
            idx_n <- idx_vec[j+1]
            if (idx_n - idx_c == 1){
                rept_vec <- c(rept_vec, "1")
            }else{
                num_rep <- c(num_rep, j)
                rept <- substr(cigar_string, idx_c + 1, idx_n - 1)
                rept_vec <- c (rept_vec, rept)
            }
        }
    }
    idx_l <- idx_vec[length(idx_vec)]
    if (idx_l < nchar(cigar_string)){
        rept_l <- substr(cigar_string, idx_l+1, nchar(cigar_string))
        rept_vec <- c(rept_vec, rept_l)
    }
    expand_cigar <- c()
    for (k in seq(1, length(rept_vec))) {
        k_rep <- as.integer(rept_vec[k])
        char_ex <- char_vec[k]
        expand_cigar <- c (expand_cigar, rep(char_ex, k_rep))
    }
    if (length(char_vec) > length(rept_vec)){
        expand_cigar <- c(expand_cigar, char_vec[length(char_vec)])
    }
    expand_cigar <- paste0(expand_cigar, collapse = "")
    return (expand_cigar)
}



add_footprints_on_edges <- function(methylation_string = "..Z..x..Z...",
                                    footprint_string = "...FFFFF....") {
    len_vec <- nchar(methylation_string)
    left_border <- 0
    right_border <- 0
    cnt <- 0
    for (i in seq_len(len_vec)) {
        char <- substr(methylation_string, i, i)
        if (char == "." || char == tolower(char)) {
            cnt <- cnt + 1
        } else if (char == toupper(char)) { break }
    }
    left_border <- cnt
    cnt <- 0
    for (i in seq(len_vec, 1, -1)) {
        char <- substr(methylation_string, i, i)
        if (char == "." || char == tolower(char)) {
            cnt <- cnt + 1
        } else if (char == toupper(char)) { break }
    }
    right_border <- len_vec - cnt + 1
    if (left_border == len_vec) {
        left_border <- 0
        right_border <- 0
    }
    rep_v <- "F"
    if (left_border == 0  && right_border == 0) {
        o_str <- paste0(rep(rep_v, len_vec), collapse = "")
    } else {
        if (left_border > 0 && right_border > 0) {
            o_str <- paste0(
                paste(rep(rep_v, left_border), collapse = ""),
                substr(footprint_string, (left_border + 1),
                       (right_border - 1) ),
                paste(rep(rep_v, (len_vec - right_border + 1) ),
                      collapse = ""))
        }else{
            if (left_border > 0 && right_border == 0){
                o_str <- paste0(
                    paste(rep(rep_v, left_border), collapse = ""),
                    substr(footprint_string, (left_border + 1) , len_vec),
                    collapse ="")
            }else{
                if (left_border == 0 && right_border > 0){
                    o_str <- paste0(substr(footprint_string,1,(right_border-1)),
                                    paste(rep(rep_v,
                                              (len_vec - right_border + 1))),
                                    collapse = "")
                }
            }
        }
    }
    return(o_str)
}

get_count_and_percentage_methylation <-function (methylation_vector = ""){
    total <- 0
    total_lower <- 0
    total_upper <- 0
    per_c <- 0
    total_length <- str_length(methylation_vector)
    for (idx_char in seq (1, total_length) ){
        char <- str_sub(methylation_vector, idx_char, idx_char)
        if (char == "." || char == "M"){ # M is just the fill-in (M: missing)
            next
        }else if (tolower(char) == char){
            total  <- total + 1
            total_lower <- total_lower + 1
            
        }else {
            total <- total + 1
            total_upper <- total_upper + 1
        }
    }
    if (total > 0){
        per_c <- as.character(round (total_upper/total, 3) * 100)
    }else {
        per_c <- "NA"
    }
    return (list (per_c = per_c, total_upper = total_upper,
                  total_lower = total_lower, total = total))
}

capital_percentage_and_stretch_of_unbound <- function (methylation_vector = ""){
    total_length <- str_length(methylation_vector)
    check_vector <- paste0(rep(".", total_length), collapse = "")
    total_capital <- 0
    total_small <- 0
    left_extreme <- 1
    right_extreme <- total_length
    first_found <- FALSE
    base_counter <- -1
    last_capital <- -1
    for (idx_char in seq (1, total_length)){
        char <- str_sub(methylation_vector, idx_char, idx_char)
        base_counter <- base_counter + 1
        
        if (char == "."){
            next
        } else if (char == tolower(char)) {
            total_small <- total_small + 1
        } else if (char == toupper(char)){
            total_capital <- total_capital + 1
            if (first_found == FALSE){
                left_extreme <- base_counter
                first_found = TRUE
            } else {
                last_capital <- base_counter
            }
        } else {
            next
        }
        
    }
    naked_dna_stretches = c()
    if (first_found == TRUE) {
        if (last_capital == -1){
            naked_dna_stretches <- c (naked_dna_stretches, left_extreme)
            naked_dna_stretches <- c (naked_dna_stretches, 0)
            naked_dna_stretches <- c (right_extreme - left_extreme)
        } else {
            naked_dna_stretches <- c (naked_dna_stretches, left_extreme)
            naked_dna_stretches <- c (naked_dna_stretches,
                                      last_capital - left_extreme) # check if -1 is required
            naked_dna_stretches <- c (naked_dna_stretches,
                                      right_extreme - last_capital) # check if -1 is required
        }
    }else {
        naked_dna_stretches = c(0,0,0)
        
    }
    total <- total_capital + total_small
    percentage_capital <- "NA"
    if (total > 0) {
        percentage_capital <- round (total_capital/total, 2)
    }
    return (list (total_length = total_length,
                  total = total,
                  percentage_capital = percentage_capital,
                  naked_dna_stretches = naked_dna_stretches))
}

get_percent_for_each_footprint <- function (footprint_vector = ""){
    flen_vec <- c()
    cnt <- 0
    f_switch <- FALSE
    total_length <- str_length(footprint_vector)
    for (idx_char in seq (1, total_length)){
        char <- str_sub(footprint_vector, idx_char, idx_char)
        if (char == "F"){
            cnt <- cnt + 1
            f_switch = TRUE
        }else{
            if (cnt != 0){
                per <- round ((cnt/total_length)*100, 2)
                cnt <- 0
                flen_vec <- c (flen_vec, per)
            }
        }
        
    }
    if (cnt != 0){
        per <- round ((cnt/total_length)*100, 2)
        flen_vec <- c (flen_vec, per)
    }
    return (flen_vec)
}

assign_binding_states <- function (overlap_string = "",
                                   lextend = 30, rextend = 30,
                                   fp_cap = 50){
    all_info <- unlist (str_split(overlap_string, pattern = "\t"))
    labels_for_info = c("chrQ", "abs_roi_start", "abs_roi_end",
                        "chrS", "abs_read_start", "abs_read_end",
                        "read_id", "score", "strand", "complete_fp_str",
                        "complete_m_str")
    chrQ <- all_info[1]
    abs_roi_start <- as.integer(all_info[2])
    abs_roi_end <- as.integer(all_info[3])
    chrS <- all_info[4]
    abs_read_start <- as.integer(all_info[5])
    abs_read_end <- as.integer(all_info[6])
    read_id <- all_info[7]
    score <- all_info[8]
    strand <- all_info[9]
    complete_fp_str <- all_info[10]
    complete_m_str <- all_info[11]
    rel_roi_start <- abs_roi_start - abs_read_start + 1
    rel_roi_end <- abs_roi_end - abs_read_start  + 1
    intersected_fp_str <- str_sub (complete_fp_str, rel_roi_start, rel_roi_end)
    intersected_m_str <- str_sub (complete_m_str, rel_roi_start, rel_roi_end)
    
    #stop("here")
    fp_info_query = fp_info_for_sub_molecule (
        query_footprint = intersected_fp_str,
        query_start = rel_roi_start,
        query_end = rel_roi_end,
        complete_footprint = complete_fp_str)
    total_length <- str_length(complete_fp_str)
    half_of_query <- as.integer ((rel_roi_end - (rel_roi_start - 1))/2)
    left_neighbor <- extend_fp_from_center(
        sub_moleculue = intersected_m_str,
        complete_molecule = complete_m_str,
        abs_roi_start = abs_roi_start,
        abs_roi_end = abs_roi_end,
        abs_complete_molecule_start = abs_read_start,
        abs_complete_molecule_end = abs_read_end,
        lextend = lextend,
        rextend = half_of_query, strand = strand)
    right_neighbor <- extend_fp_from_center(
        sub_moleculue = intersected_m_str,
        complete_molecule = complete_m_str,
        abs_roi_start = abs_roi_start,
        abs_roi_end = abs_roi_end,
        abs_complete_molecule_start = abs_read_start,
        abs_complete_molecule_end = abs_read_end,
        lextend = half_of_query,
        rextend = rextend, strand = strand)
    info_for_left <- get_count_and_percentage_methylation(left_neighbor)
    info_for_right <- get_count_and_percentage_methylation(right_neighbor)
    percentage_methylation <- "NA"
    per_capital_l <- info_for_left$per_c
    total_upper_l <- info_for_left$total_upper
    total_lower_l <- info_for_left$total_lower
    total_l <- info_for_left$total
    per_capital_r <- info_for_right$per_c
    total_upper_r <- info_for_right$total_upper
    total_lower_r <- info_for_right$total_lower
    total_r <- info_for_right$total
    
    if (per_capital_l !="NA" && per_capital_r !="NA"){
        if (as.numeric (per_capital_l) <= as.numeric(per_capital_r)) {
            percentage_methylation <- per_capital_l
        }else {
            percentage_methylation <- per_capital_r
        }
    }else if (per_capital_l == "NA" && per_capital_r != "NA") {
        percentage_methylation <- per_capital_r
    } else if (per_capital_l != "NA" && per_capital_r == "NA"){
        percentage_methylation <- per_capital_l
    } else{
        percentage_methylation <- "0"
    }
    percentage_methylation <- as.numeric(percentage_methylation)
    
    #   return (list (per_c = per_c, total_upper = total_upper,
    #                  total_lower = total_lower, total = total))
    edge_info <- capital_percentage_and_stretch_of_unbound(complete_m_str)
    complete_m_str_len <- edge_info$total_length
    total_letters_on_complete_m_str <- edge_info$total
    ratio_capital_letters <- edge_info$percentage_capital
    edges <- edge_info$naked_dna_stretches
    max_of_edges <- max (edges[1], edges[2])
    intersect_len <- str_length(intersected_fp_str)
    indiv_fp_percent <- get_percent_for_each_footprint(intersected_fp_str)
    total_percent_fp_in_intersect <- round ((
        str_count(intersected_fp_str,
                  "F")/intersect_len)*100, 3)
    labels_on_reads <- c()
    label <- "NA"
    labelled <- "NA"
    label_list <- list ("Naked_DNA" = 0, "TF" = 1,
                        "Nuc" = 2, "discard" = 3)
    if (length(indiv_fp_percent) == 0){
        labelled <- paste(overlap_string, intersected_fp_str, "0.0", "0",
                          label_list$Naked_DNA, total_length, sep = "\t")
        label <- "Naked_DNA"
    }else{
        max_contribution_in_the_intersect <- which.max(indiv_fp_percent)
        max_percentage_overlap <- indiv_fp_percent[
            max_contribution_in_the_intersect]
        length_of_max_contribution <- fp_info_query$fp_len_query[
            max_contribution_in_the_intersect]
        abs_start_loc_of_max <- fp_info_query$fp_start_loci_on_query[
            max_contribution_in_the_intersect]
        abs_fp_start_on_complete_fp <- fp_info_query$abs_fp_start_in_query[
            abs_start_loc_of_max]
        if (total_percent_fp_in_intersect <= 30 &&
            percentage_methylation > 25){
            labelled = paste(overlap_string, intersected_fp_str,
                             max_percentage_overlap,
                             length_of_max_contribution,
                             label_list$Naked_DNA, total_length, sep = "\t")
            label <- "Naked_DNA"
        }else if (total_percent_fp_in_intersect <=30 &&
                  percentage_methylation <=25){
            labelled <- paste(overlap_string, intersected_fp_str,
                              max_percentage_overlap,
                              length_of_max_contribution,
                              label_list$Nuc, total_length, sep = "\t")
            label <- "Nuc"
        }else if (total_percent_fp_in_intersect > 30 &&
                  length_of_max_contribution <=fp_cap){
            
            if ( (abs_fp_start_on_complete_fp +
                  length_of_max_contribution - 1) == total_length
                 || abs_fp_start_on_complete_fp == 1) {
                labelled <- paste(overlap_string,
                                  intersected_fp_str, max_percentage_overlap,
                                  length_of_max_contribution,
                                  label_list$discard, total_length, sep = "\t")
                label <- "discard"
            }else {
                labelled <- paste(overlap_string, intersected_fp_str,
                                  max_percentage_overlap,
                                  length_of_max_contribution,
                                  label_list$TF, total_length, sep = "\t")
                label <- "TF"
            }
        }else if (length_of_max_contribution > fp_cap){
            labelled <- paste(overlap_string, intersected_fp_str,
                              max_percentage_overlap,
                              length_of_max_contribution,
                              label_list$Nuc, total_length, sep = "\t")
            label <- "Nuc"
        }
    }
    label_id = label_list[[label]]
    return (list(labelled = labelled,
                 label = label,
                 label_id = label_id,
                 read_id = read_id,
                 intersected_fp_str = intersected_fp_str,
                 intersected_m_str = intersected_m_str,
                 abs_roi_start = abs_roi_start,
                 abs_roi_end = abs_roi_end,
                 abs_read_start = abs_read_start,
                 abs_read_end = abs_read_end,
                 complete_fp_str = complete_fp_str,
                 complete_m_str = complete_m_str,
                 chrQ = chrQ))
}

convert_fp_char_to_num <- function (char_vector = "MMM...FF...FFEEEEMMM"){
    total_len <- str_length(char_vector)
    num_vec <- c()
    for (i in seq (total_len)){
        char <- str_sub(char_vector, i, i)
        if (char == "."){
            num_vec <- c (num_vec, 0)
        }else if (char == "M"){
            num_vec <- c (num_vec, -1)
        }else if (char == "F" || char == "E"){
            num_vec <- c (num_vec, 1)
        }
    }
    ret = paste0(num_vec, collapse = "\t")
    return (ret)
}

convert_m_char_to_num <- function (char_vector = "MMM...x..X..H...hh...MMMM"){
    total_len <- str_length(char_vector)
    num_vec <- c()
    for (i in seq (total_len)){
        char <- str_sub(char_vector, i, i)
        if (char == "."){
            num_vec <- c(num_vec, 0)
        }else if (char == "M"){
            num_vec <- c (num_vec, -1)
        } else if (char == tolower(char)){
            num_vec <- c (num_vec, 1)
        }else {
            num_vec <- c (num_vec, 2)
        }
    }
    ret = paste0(num_vec, collapse = "\t")
    return (ret)
}
saveToPDF <- function(...) {
    d <- dev.copy(pdf,...);
    dev.off(d);
}

saveToPNG <- function(...) {
    d <- dev.copy(png,...);
    dev.off(d);
}

saveToEPS <- function(...) {
    d <- dev.copy(postscript,...);
    dev.off(d);
}
label_list <- list ("Naked_DNA" = 0, "TF" = 1,
                    "Nuc" = 2, "discard" = 3)
getHubURL <- function (organism = "dmelanogaster",
                       model = "S2", condition = "WT", type = "dSMF",
                       genome_assembly = "dm6", tr = "fp_and_mvec"){
    row.names (hubs) <- paste(hubs$Organism, hubs$model,
                              hubs$Condition, hubs$Type,
                              hubs$Genome, hubs$Track, sep = "@")
    
    identifier <- paste(organism, model, condition, type,
                        genome_assembly, tr, sep = "@")
    track_url <- hubs[identifier, "Hub_URL"]
    track_url <- paste0("https://api.genome.ucsc.edu/getData/track?hubUrl=",
                        track_url)
    return (track_url)
}
remove_duplicate_molecules <- function (fp_file = "peak229.all_fp.bed",
                                        label = "peak229", target_dir = ""){
    # the data is coming from big bed file, so the coordinates
    # are ought to be sorted
    dt <- read.table(fp_file, sep = "\t", header = FALSE,
                     stringsAsFactors = FALSE)
    dt ["srr"] <-  unlist (lapply (dt$V4, function (x) {
        unlist (strsplit (x, split = "\\."))[[1]]}))
    dt["mol"] <- unlist (lapply (dt$V4,
                                 function (x) {unlist (
                                     strsplit (x, split = "`"))[[2]]}))
    idx_duplicate <- which (duplicated(dt[, c ("V1", "V2", "V3",
                                               "srr", "mol")]))
    uniq_mol_df = dt[-idx_duplicate, ]
    dup_mol_df <- dt[idx_duplicate, ]
    uniq_mol_df$srr <- NULL
    uniq_mol_df$mol <- NULL
    dup_mol_df$srr <- NULL
    dup_mol_df$mol <- NULL
    
    subject_file <- paste(target_dir, paste0(label, ".dup_removed.bed"),
                          sep = "/")
    dup_file <- paste(target_dir, paste0(label, ".dup.bed"), sep = "/")
    write.table(uniq_mol_df, file = subject_file, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(dup_mol_df, file = dup_file, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    return(subject_file)
    
}
#' Visualize Single Molecule Footprint Patterns
#'
#' @description Generates footprint plots from NOME-seq data within specified genomic regions.
#' Designed for Drosophila melanogaster analysis with expandable parameters for other organisms.
#'
#' @param organism Organism code (default: "dmelanogaster")
#' @param model Biological model/system (default: "S2" cells)
#' @param condition Experimental condition (default: "WT")
#' @param genome_assembly Genome version (default: "dm6")
#' @param type Data type ("dSMF" = dual enzyme Single Molecule Footprinting)
#' @param chromosome Chromosome ID (e.g., "chr2L")
#' @param start Genomic start position (numeric/character)
#' @param end Genomic end position (numeric/character)
#' @param tr Track name for BigBed file resource (default: "fp_and_mvec")
#' @param label Plot title annotation
#' @param span_left Upstream window size from region (default: 150)
#' @param span_right Downstream window size from region (default: 150)
#' @param remove_dup Remove duplicate reads? (default: FALSE)
#' @param fp_cap Maximum footprint value for y-axis scaling (default: 50)
#' @param target_dir destination directory where outputs will be stored
#'
#' @return Saves a heatmap in pdf, png, and eps file formats. Also, occupancies are saved in a tsv file
#'
#' @details
#' This function retrieves and visualizes DNA accessibility patterns from single-molecule data:
#' - Integrates with UCSC-style track hubs via `tr` parameter
#' - Automatically handles coordinate conversion for dm6 genome
#' - Implements dynamic y-axis scaling using `fp_cap`
#'
#' @examples
#' # Basic usage with default parameters
#' plotFootprints()
#'
#' # Custom genomic region analysis
#' plotFootprints(
#'   chromosome = "chr2L", start = 480290, end = 480320, label = "peak229")
#'
#' @export
plotFootprints <- function(
        organism = "dmelanogaster", model = "S2",
        condition = "WT", genome_assembly = "dm6",
        type = "dSMF", chromosome = "chr2L",
        start = 480290, end = 480320,
        tr = "fp_and_mvec", label = "peak229",
        span_left = 150, span_right = 150, remove_dup = FALSE,
        fp_cap = 50, target_dir = "") {
    # Function implementation
}

plotFootprints <- function (organism = "dmelanogaster", model = "S2",
                            condition = "WT", genome_assembly = "dm6",
                            type = "dSMF", chromosome = "chr2L",
                            start = 480290, end = 480320,
                            tr = "fp_and_mvec", label = "peak229",
                            span_left = 150, span_right = 150,
                            remove_dup = FALSE, fp_cap = 50, target_dir = ""){
    
    # get JSON first
    #print (track)
    if (target_dir == ""){
        target_dir = tempdir()
    }else if (!dir.exists(target_dir)){
        dir.create(target_dir, recursive = TRUE)
    }
    track_url <- getHubURL(organism = organism, model = model,
                           condition = condition, type = type,
                           genome_assembly = genome_assembly, tr = tr)
    #print(track_url)
    json_file <- downloadJson(track_url = track_url, chromosome = chromosome,
                              assembly = genome_assembly, st = start,
                              en = end, label = label, tr = tr, 
                              target_dir = target_dir)
    if (is.na(json_file)){
        print ("Download from UCSC REST API failed!")
        return (NA)
    }else {
    
        tmp_file <- JSON_to_Bed(json_file = json_file, label = label, tr = tr,
                                target_dir = target_dir)
        subject_file <- NA
        if (remove_dup == TRUE){
            subject_file <- remove_duplicate_molecules(fp_file = tmp_file,
                                                       label = label, target_dir = target_dir)
        }else{
            subject_file <- tmp_file
        }
        query_file <- prepareQueryFile(chromosome = chromosome,
                                       start = start, end = end, 
                                       label = label, target_dir = target_dir)
        
        overlap_df <- intersectQueryandSubject(query_file = query_file,
                                               subject_file = subject_file)
        # successful execution of this function will ensure that `overlap` file is
        # generated
        # get the second last column of the overlap_df and parse it using "|" as
        # delimiter. Example:
        # SRR3133329.41244179_41244179/1_adjacent`99~147|.131F52.77F35.5|.19h.33z.30h.14h.9z.20Z.4z.47H.3X.4H.Z.8Z.H.5Z.8X.9Z.18H.9Z.18x.16H.4
        # 2nd field is footprint, 3rd field is methylation vector
        fp_and_mvec_info <- overlap_df[, (length(overlap_df) - 1)]
        read_details <- unlist (lapply(fp_and_mvec_info,
                                       function (x) {
                                           yt = unlist(strsplit(x, split = "\\|")
                                           )[1] }))
        fp_vec <- unlist (lapply(fp_and_mvec_info,
                                 function (x) {
                                     yt = unlist(strsplit(x, split = "\\|")
                                     )[2] }))
        m_vec <- unlist (lapply(fp_and_mvec_info,
                                function (x) {
                                    yt = unlist(strsplit(x, split = "\\|"))[3] }))
        expanded_fp <- unlist(lapply(fp_vec, expand_cigar_string))
        expanded_m_vec <- unlist (lapply(m_vec, expand_cigar_string))
        combined <- data.frame (expanded_m_vec = expanded_m_vec,
                                expanded_fp = expanded_fp)
        with_edges <- unlist (apply(combined, 1,
                                    function(x){
                                        add_footprints_on_edges(x[1], x[2])} ))
        processed_df <- overlap_df[, c(1,2,3,4,5,6)]
        processed_df$read_details <- read_details
        processed_df$score <- 1.0
        processed_df$strand <- "."
        processed_df$expanded_footprint <- with_edges
        processed_df$expanded_m_vec <- expanded_m_vec
        length_and_read_info_df <- writeBindingStatesFiles (
            processed_df = processed_df,
            label = label,fp_cap = fp_cap,
            span_left = span_left, span_right = span_right,
            target_dir = target_dir)
        ordered_by_length <- length_and_read_info_df[with (length_and_read_info_df,
                                                           order(binding_state_vec,
                                                                 flag_vec,
                                                                 abs_start_vec)),]
        
        writeIntermediateFiles(ordered_by_length = ordered_by_length,
                               label = label, target_dir = target_dir)
        
        plot_title = paste0(organism, " ", model, " ", condition)
        x_label = paste0(label, " (", genome_assembly, " ", 
                         chromosome, ":", start, "-", end, ")")
        pdf_file_path = savePlotSMForDSMF(label = label, span_left = span_left, 
                                          span_right = span_right, plot_title = plot_title,
                                          x_label = x_label, target_dir = target_dir)
        #removing the intermediate files:
        deleteIntermediates(label = label, target_dir = target_dir)
        return (pdf_file_path) 
    }
}


#' Visualize Single Molecule Footprint Patterns using Local BigBed file
#'
#' @description Generates footprint plots from SMF/dSMF data within specified genomic regions.
#' Designed for Drosophila melanogaster analysis with expandable parameters for other organisms.
#'
#' @param bigBed BigBed file path (e.g., "demo.bb")
#' @param organism Organism code (default: "dmelanogaster")
#' @param model Biological model/system (default: "S2" cells)
#' @param condition Experimental condition (default: "WT")
#' @param genome_assembly Genome version (default: "dm6")
#' @param chromosome Chromosome ID (e.g., "chr2L")
#' @param start Genomic start position (numeric/character)
#' @param end Genomic end position (numeric/character)
#' @param label Plot title annotation
#' @param span_left Upstream window size from region (default: 150)
#' @param span_right Downstream window size from region (default: 150)
#' @param remove_dup Remove duplicate reads? (default: FALSE)
#' @param fp_cap Maximum footprint value for y-axis scaling (default: 50)
#' @param target_dir destination directory where outputs will be stored
#'
#' @return Saves a heatmap in pdf, png, and eps file formats. Also, occupancies are saved in a tsv file
#'
#' @details
#' This function uses local SMF/dSMF bigBed file and plot heatmap
#'
#'
#' @examples
#' # Basic usage with default parameters
#' plotFootprintsUsingLocalBigBed()
#'
#' # Custom genomic region analysis
#' plotFootprintsUsingLocalBigBed(
#'   bigBed = "inst/extdata/demo.bb", organism = "dmelanogaster",
#'   model = "S2", condition = "WT", genome_assembly = "dm3",
#'   chromosome = "chr2L", start = 480290, end = 480320, label = "peak229")
#'
#' @export
plotFootprintsUsingLocalBigBed <- function(
        bigBed = "inst/extdata/demo.bb", organism = "dmelanogaster", 
        model = "S2", condition = "WT", genome_assembly = "dm3",
        chromosome = "chr2L",start = 480290, end = 480320,
        label = "peak229", span_left = 150, span_right = 150,
        remove_dup = FALSE,fp_cap = 50, target_dir = "") {
    # Function implementation
}


plotFootprintsUsingLocalBigBed <- function (
        bigBed = "", organism = "dmelanogaster",
        model = "S2", condition = "WT", genome_assembly = "dm3",
        chromosome = "chr2L",start = 480290, end = 480320,
        label = "peak229", span_left = 150, span_right = 150,
        remove_dup = FALSE,fp_cap = 50, target_dir = "") {
    if (bigBed == "") {
        cat ("'bigBed' arugment is empty. bigBed file file path requied\n")
        return ("bigBed file path")
    }else if (!file.exists(bigBed)){
        cat ("bigBed file path requied\n")
        return ("bigBed file path")
    }
    if (target_dir == ""){
        target_dir = tempdir()
    }else if (!dir.exists(target_dir)){
        dir.create(target_dir, recursive = TRUE)
    }
    query_file <- prepareQueryFile(chromosome = chromosome,
                                   start = start, end = end, 
                                   label = label, target_dir = target_dir)
    
    overlap_df <- overlapUsingBigBed(bigbed_file = bigBed,
                                     query_file = query_file)
    # successful execution of this function will ensure that `overlap` file is
    # generated
    # get the second last column of the overlap_df and parse it using "|" as
    # delimiter. Example:
    # SRR3133329.41244179_41244179/1_adjacent`99~147|.131F52.77F35.5|.19h.33z.30h.14h.9z.20Z.4z.47H.3X.4H.Z.8Z.H.5Z.8X.9Z.18H.9Z.18x.16H.4
    # 2nd field is footprint, 3rd field is methylation vector
    fp_and_mvec_info <- paste(overlap_df$name, overlap_df$field8, sep = "|")
    read_details <- unlist (lapply(fp_and_mvec_info,
                                   function (x) {
                                       yt = unlist(strsplit(x, split = "\\|")
                                       )[1] }))
    fp_vec <- unlist (lapply(fp_and_mvec_info,
                             function (x) {
                                 yt = unlist(strsplit(x, split = "\\|")
                                 )[2] }))
    m_vec <- unlist (lapply(fp_and_mvec_info,
                            function (x) {
                                yt = unlist(strsplit(x, split = "\\|"))[3] }))
    expanded_fp <- unlist(lapply(fp_vec, expand_cigar_string))
    expanded_m_vec <- unlist (lapply(m_vec, expand_cigar_string))
    combined <- data.frame (expanded_m_vec = expanded_m_vec,
                            expanded_fp = expanded_fp)
    with_edges <- unlist (apply(combined, 1,
                                function(x){
                                    add_footprints_on_edges(x[1], x[2])} ))
    processed_df <- overlap_df[, c(1,2,3,4,5,6)]
    #print (head (processed_df))
    processed_df$read_details <- read_details
    processed_df$score <- 1.0
    processed_df$strand <- "."
    processed_df$expanded_footprint <- with_edges
    processed_df$expanded_m_vec <- expanded_m_vec
    length_and_read_info_df <- writeBindingStatesFiles(
        processed_df = processed_df,
        label = label,fp_cap = fp_cap,
        span_left = span_left, span_right = span_right, 
        target_dir = target_dir)
    ordered_by_length <- length_and_read_info_df[with (length_and_read_info_df,
                                                       order(binding_state_vec,
                                                             flag_vec,
                                                             abs_start_vec)),]
    
    writeIntermediateFiles(ordered_by_length = ordered_by_length,
                           label = label, target_dir = target_dir)
    
    plot_title = paste0(organism, " ", model, " ", condition)
    x_label = paste0(label, " (", genome_assembly, " ", 
                     chromosome, ":", start, "-", end, ")")
    savePlotSMForDSMF(label = label, span_left = span_left, 
                      span_right = span_right, plot_title = plot_title, 
                      target_dir = target_dir)
    #removing the intermediate files:
    deleteIntermediates(label = label, target_dir = target_dir)
}



#' Discover Available Single Molecule Data Tracks
#'
#' @description Retrieves track metadata from a centralized Google Sheet containing
#' organism-specific single-molecule datasets and experimental conditions.
#'
#'
#' @return A tibble containing track metadata with columns:
#' \itemize{
#'   \item organism: Species identifier (e.g. "dmelanogaster")
#'   \item genome_assembly: Reference genome version
#'   \item condition: Experimental condition
#'   \item model: Biological model/system
#'   \item track_name: Identifier for `tr` parameter in `plotFootprints`
#'   \item data_type: NOME-seq protocol variant
#'   \item bigbed_url: Cloud-hosted data location
#' }
#'
#' @examples
#' # List all available tracks from default repository
#' listTracks()
#'
#'
#' @seealso
#' `plotFootprints()`: Use track names from this list to visualize specific datasets
#'
#' @export
listTracks <- function() {
    print(hubs)
}


listTracks <- function(){
    
    print(hubs)
}

############################## Nanopore ########################################


#' Visualize Single Molecule Footprint Patterns
#'
#' @description Generates footprint plots from NOME-seq data within specified genomic regions.
#' Designed for Drosophila melanogaster analysis with expandable parameters for other organisms.
#'
#' @param organism Organism code (default: "scerevisiae")
#' @param model Biological model/system (default: "BY4741 strain")
#' @param condition Experimental condition (default: "WT")
#' @param genome_assembly Genome version (default: "sacCer3")
#' @param type Data type ("Nanopore" = dual enzyme Single Molecule Footprinting)
#' @param chromosome Chromosome ID (e.g., "chrIII")
#' @param start Genomic start position (numeric/character)
#' @param end Genomic end position (numeric/character)
#' @param tr Track name for BigBed file resource (default: "nanopore_meth_calls")
#' @param label Plot title annotation (default: "smac_seq")
#' @param span_left Upstream window size from region (default: 1000)
#' @param span_right Downstream window size from region (default: 1000)
#' @param stride sliding window used for singal aggregation 
#' @param target_dir destination directory where outputs will be stored
#'
#' @return Saves a heatmap in pdf, png, and eps file formats. Also, occupancies are saved in a tsv file
#'
#' @details
#' This function retrieves and visualizes DNA accessibility patterns from single-molecule data:
#' - Integrates with UCSC-style track hubs via `tr` parameter
#' - Automatically handles coordinate conversion for dm6 genome
#' - Implements dynamic y-axis scaling using `fp_cap`
#'
#' @examples
#' # Basic usage with default parameters
#' plotMethylationCallsNanopore()
#'
#' # Custom genomic region analysis
#'
#' plotMethylationCallsNanopore (
#'         organism = "scerevisiae", model = "BY4741 strain",
#'         condition = "WT", genome_assembly = "sacCer3",
#'         type = "SMF", chromosome = "chrIII",
#'         start = 114300, end = 114600,
#'         tr = "nanopore_meth_calls", label = "smac_seq",
#'         span_left = 1000, span_right = 1000, stride = 5, target_dir = "")


#' @export
plotMethylationCallsNanopore <- function(
        organism = "scerevisiae", model = "BY4741 strain",
        condition = "WT", genome_assembly = "sacCer3",
        type = "SMF", chromosome = "chrIII",
        start = 114300, end = 114600,
        tr = "nanopore_meth_calls", label = "smac_seq",
        span_left = 1000, span_right = 1000, stride = 5, target_dir = "") {
    # Function implementation
}


plotMethylationCallsNanopore <- function(
        organism = "scerevisiae", model = "BY4741 strain",
        condition = "WT", genome_assembly = "sacCer3",
        type = "SMF", chromosome = "chrIII",
        start = 114300, end = 114600,
        tr = "nanopore_meth_calls", label = "smac_seq",
        span_left = 1000, span_right = 1000, stride = 5, target_dir = "") {
    
    if (target_dir == ""){
        target_dir = tempdir()
    }else if (!dir.exists(target_dir)){
        dir.create(target_dir, recursive = TRUE)
    }
    
    track_url <- getHubURL(organism = organism, model = model,
                           condition = condition, type = type,
                           genome_assembly = genome_assembly, tr = tr)
    #print(track_url)
    json_file <- downloadJson(track_url = track_url, chromosome = chromosome,
                              assembly = genome_assembly, st = start,
                              en = end, label = label, tr = tr, 
                              target_dir = target_dir)
    if (is.na(json_file)){
        print("JSON download from UCSC REST API failed")
        return (NA)
    } else {
        subject_file <- JSON_to_Bed_for_Nanopore(json_file = json_file,
                                                 label = label, tr = tr,
                                                 target_dir = target_dir)
        query_file <- prepareQueryFile(chromosome = chromosome,
                                       start = start, end = end, 
                                       label = label, target_dir = target_dir)
        
        overlap_df <- intersectQueryandSubject(query_file = query_file,
                                               subject_file = subject_file)
        
        overlap_df_for_plot = overlap_df[, c(1, 2, 3,4,5,6,7,8)]
        read_id <- unlist (lapply(overlap_df$name,
                                  function (x) {
                                      yt = unlist(strsplit(x, split = "\\|")
                                      )[1] }))
        methylation_call <- unlist (lapply(overlap_df$name,
                                           function (x) {
                                               yt = unlist(strsplit(x, split = "\\|")
                                               )[2] }))
        overlap_df_for_plot$name = read_id
        overlap_df_for_plot$score = overlap_df$score
        overlap_df_for_plot$methylation_call = methylation_call
        
        write.table(overlap_df_for_plot, 
                    file = paste(target_dir, paste0(label, ".intersect.bed"), 
                                 sep = "/"),
                    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
        data_to_plot_and_strand_map_df = generatePlotMatrix(overlap_df_for_plot,
                                                            span_left = span_left,
                                                            span_right = span_right)
        data_to_plot = data_to_plot_and_strand_map_df$out_mat
        strand_map_df = data_to_plot_and_strand_map_df$strand_map_df
        pos_strand_reads = names (strand_map_df[strand_map_df$strand == "+", ])
        neg_strand_reads = names (strand_map_df[strand_map_df$strand == "-", ])
        
        center = as.integer(dim(data_to_plot)[2]/2)
        half_width = as.integer((as.integer(end) - as.integer(start))/(2*stride) )
        l_index = center - half_width
        r_index = center + half_width
        sub_df = data_to_plot[, seq (l_index, r_index)]
        sum_by_row = apply(sub_df, 1, sum)
        tmp = strand_map_df
        
        ordered_by_methylation = sum_by_row[order(sum_by_row)]
        tmp$cnt = ordered_by_methylation[row.names(strand_map_df)]
        ordered_tmp = tmp[order(factor (tmp$strand, levels = c("-", "+")), tmp$cnt),]
        data_to_plot = data_to_plot[row.names(ordered_tmp),]
        write.table (table (factor(strand_map_df$strand, levels = c("-", "+"))),
                     file = paste (target_dir, 
                                   paste0(label, ".strand_wise_count.tsv"), sep = "/"),
                     sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table (data_to_plot, file = paste(target_dir, 
                                                paste0(label, ".nanopore_methylation.tsv"), sep = "/"),
                     sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
        plot_title = paste0(organism, " ", model, " ", condition)
        x_label = paste0(label, " (", genome_assembly, " ", 
                         chromosome, ":", start, "-", end, ")")
        savePlotNanopore(label = label, plot_title = plot_title, x_label = x_label,
                     target_dir = target_dir)
    }
}





#' Visualize aggregate methylation calls from ONT sequencing data
#'
#' @description Generates heatmap of methylation calls from likes of SMAC-seq datasets.
#' currently has data from Yeast (sacCer3)
#'
#' @param bigBed BigBed file path (e.g., "demo.bb")
#' @param chromosome Chromosome ID (e.g., "chr2L")
#' @param start Genomic start position (numeric/character)
#' @param end Genomic end position (numeric/character)
#' @param label Plot title annotation
#' @param span_left Upstream window size from region (default: 150)
#' @param span_right Downstream window size from region (default: 150)
#' @param remove_dup Remove duplicate reads? (default: FALSE)
#' @param stride a window size used for averaging the methylation call values
#' @param target_dir destination directory where outputs will be stored
#'
#' @return Saves a methylation calls (binary) heatmap in pdf, png, and eps file formats.
#'
#' @details
#' This function uses local SMF/dSMF bigBed file and plot heatmap
#'
#'
#'
#' @export
plotMethylationCallsNanoporeUsingLocalBigBed <- function(
        bigBed = "", 
        chromosome = "chrIII",start = 114300, end = 114600,
        label = "smac_seq", span_left = 1000, span_right = 1000,
        remove_dup = FALSE, stride = 5, target_dir = "") {
    # Function implementation
}


plotMethylationCallsNanoporeUsingLocalBigBed <- function (
        bigBed = "",
        chromosome = "chrIII",start = 114300, end = 114600,
        label = "smac_seq", span_left = 1000, span_right = 1000,
        remove_dup = FALSE, stride = 5, target_dir = "") {
    
    if (bigBed == "") {
        cat ("'bigBed' arugment is empty. bigBed file file path requied\n")
        return ("bigBed file path")
    }else if (!file.exists(bigBed)){
        cat ("bigBed file path requied\n")
        return ("bigBed file path")
    }
    if (target_dir == ""){
        target_dir = tempdir()
    }else if (!dir.exists(target_dir)){
        dir.create(target_dir, recursive = TRUE)
    }
    
    
    query_file <- prepareQueryFile(chromosome = chromosome,
                                   start = start, end = end,
                                   label = label, target_dir = target_dir)
    
    overlap_df <- overlapUsingBigBed(bigbed_file = bigBed,
                                     query_file = query_file)
    
    write.table(overlap_df, file = paste(target_dir, 
                                         paste0(label, ".intersect.bed"), sep = "/"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    data_to_plot_and_strand_map_df = generatePlotMatrix(overlap_df,
                                                        span_left = span_left,
                                                        span_right = span_right)
    data_to_plot = data_to_plot_and_strand_map_df$out_mat
    strand_map_df = data_to_plot_and_strand_map_df$strand_map_df
    pos_strand_reads = names (strand_map_df[strand_map_df$strand == "+", ])
    neg_strand_reads = names (strand_map_df[strand_map_df$strand == "-", ])
    
    center = as.integer(dim(data_to_plot)[2]/2)
    half_width = as.integer((as.integer(end) - as.integer(start))/(2*stride) )
    l_index = center - half_width
    r_index = center + half_width
    sub_df = data_to_plot[, seq (l_index, r_index)]
    sum_by_row = apply(sub_df, 1, sum)
    tmp = strand_map_df
    
    ordered_by_methylation = sum_by_row[order(sum_by_row)]
    tmp$cnt = ordered_by_methylation[row.names(strand_map_df)]
    ordered_tmp = tmp[order(factor (tmp$strand, levels = c("-", "+")), tmp$cnt),]
    data_to_plot = data_to_plot[row.names(ordered_tmp),]
    write.table (table (factor(strand_map_df$strand, levels = c("-", "+"))),
                 file = paste(target_dir, paste0(label, ".strand_wise_count.tsv"),
                              sep = "/"),
                 sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table (data_to_plot, file = paste(target_dir, paste0(label, 
                                                               ".nanopore_methylation.tsv"), sep = "/"),
                 sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
    savePlotNanopore(label = label, target_dir = target_dir)
}
