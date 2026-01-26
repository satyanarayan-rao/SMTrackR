writeBindingStatesFiles <- function (
        processed_df = NULL,
        label = "peak229",
        fp_cap = 50, span_left = 150, span_right = 150,
        target_dir = ""){
    
    ofp <- file (paste(target_dir, 
                       paste0(label, "_assigned_states.tsv"), sep = "/"), "w")
    ofp_verbose <- file (paste(target_dir, 
                               paste0(label, "_verbose.tsv"), sep = "/"), "w")
    max_length_vec <- c()
    read_info_vec <- c()
    abs_start_vec <- c()
    fp_vec <- c ()
    m_vec <- c()
    binding_state_vec <- c()
    flag_vec <- c()
    
    for (idx in seq (nrow(processed_df)) ){
        line <- paste (processed_df[idx, ], collapse = "\t")
        #print (line)
        assigned_states <- assign_binding_states(line, fp_cap = fp_cap)
        #print (assigned_states)
        cat (assigned_states$labelled, file = ofp, append = TRUE)
        cat("\n", file = ofp, append = TRUE)
        #print (assigned_states$read_id)
        ## Generate the fp extended file
        sub_molecule_fp <- assigned_states$intersected_fp_str
        sub_molecule_m <- assigned_states$intersected_m_str
        complete_molecule_fp <- assigned_states$complete_fp_str
        complete_molecule_m <- assigned_states$complete_m_str
        read_id <- assigned_states$read_id
        abs_roi_start <- assigned_states$abs_roi_start
        abs_roi_end <- assigned_states$abs_roi_end
        abs_read_start <- assigned_states$abs_read_start
        abs_read_end <- assigned_states$abs_read_end
        label_id <- assigned_states$label_id
        chrQ <- assigned_states$chrQ
        marker <- paste0(chrQ,":", abs_roi_start, "-", abs_roi_end, "`",
                         read_id, "#", label_id)
        extended_fp <- extend_fp_from_center(
            sub_moleculue = sub_molecule_fp,
            complete_molecule = complete_molecule_fp,
            abs_roi_start = abs_roi_start,
            abs_roi_end = abs_roi_end,
            abs_complete_molecule_start = abs_read_start,
            abs_complete_molecule_end = abs_read_end,
            lextend = span_left, rextend = span_right)
        extended_m <- extend_fp_from_center(
            sub_moleculue = sub_molecule_m,
            complete_molecule = complete_molecule_m,
            abs_roi_start = abs_roi_start,
            abs_roi_end = abs_roi_end,
            abs_complete_molecule_start = abs_read_start,
            abs_complete_molecule_end = abs_read_end,
            lextend = span_left, rextend = span_right)
        cat(paste0(marker, "\t", extended_fp, "\t", extended_m),
            file = ofp_verbose, append = TRUE)
        cat ("\n", file = ofp_verbose, append =  TRUE)
        
        tmp <- unlist(strsplit(marker, split = "`"))
        flag_and_state <- unlist(strsplit(tmp[length(tmp)], split = "#"))
        flag <- flag_and_state[1]
        state <- flag_and_state[2]
        fp_without_m <- str_replace_all(extended_fp, "M", ".")
        fp_info <- fp_info_for_complete_molecule(fp_without_m)
        if (length(fp_info$just_length_vec) > 0){
            max_length <- max (fp_info$just_length_vec)
            which_max <- which.max (fp_info$just_length_vec)
            read_info_vec <- c(read_info_vec, marker)
            max_length_vec <- c (max_length_vec, max_length)
            abs_start_vec <- c (abs_start_vec,
                                fp_info$fp_start_locus[which_max])
            fp_vec <- c(fp_vec, extended_fp)
            m_vec <- c(m_vec, extended_m)
            flag_vec <- c(flag_vec, flag)
            binding_state_vec <- c(binding_state_vec, state)
        }else {
            read_info_vec <- c(read_info_vec, marker)
            max_length_vec <- c (max_length_vec, 0)
            abs_start_vec <- c(abs_start_vec, (0 - span_left) )
            fp_vec <- c(fp_vec, extended_fp)
            m_vec <- c(m_vec, extended_m)
            flag_vec <- c(flag_vec, flag)
            binding_state_vec <- c(binding_state_vec, state)
        }
        
    }
    close(ofp)
    close(ofp_verbose)
    length_and_read_info_df <- data.frame(read_id = read_info_vec,
                                          binding_state_vec = binding_state_vec,
                                          flag_vec = flag_vec,
                                          max_fp_length = max_length_vec,
                                          abs_start_vec = abs_start_vec,
                                          fp_vec = fp_vec,
                                          m_vec = m_vec)
    return (length_and_read_info_df)
    
    
}
