generatePlotMatrix <- function (methylation_call_df, span_left = 1000,
                                span_right = 1000, stride = 5){
    # define a zero matrix to fill
    n_molecules = dim(methylation_call_df)[1]
    roi_start = methylation_call_df[1, 2]
    roi_end = methylation_call_df[1, 3]
    total_columns = (span_left + span_right + (roi_end - roi_start))/stride + 1
    print (total_columns)
    out_mat = matrix (0L, nrow = 0, ncol = total_columns)
    read_name_vec = c ()
    strand_map_df = data.frame("strand" = c())
    strand_vec = c()
    #methylation_count_df = data.frame("read_id" = )
    for (i in seq (1, n_molecules)){
        # complete vector for a read
        spanned_vector = c()
        read_start = methylation_call_df[i, 5]
        read_end = methylation_call_df[i, 6]

        status_on_read = as.numeric(unlist(
            strsplit(methylation_call_df[i, 11], split = ",")))
        total_length = length(status_on_read)
        abs_start_for_plot = roi_start - span_left
        abs_end_for_plot = roi_end + span_right
        read_name_vec = c(read_name_vec, paste0(i,"_", methylation_call_df[i, 9]))
        strand_vec = c(strand_vec, as.character(methylation_call_df[i, 8]))
        approximate_start_offset = (abs_start_for_plot - read_start)/stride
        approximate_start_offset = as.integer(approximate_start_offset)
        approximate_end_offset = approximate_start_offset + (abs_end_for_plot - abs_start_for_plot)/stride
        approximate_end_offset = as.integer(approximate_end_offset)

        for (j in seq (approximate_start_offset, approximate_end_offset)){
            if (j < 1 || j > total_length){
                spanned_vector = c(spanned_vector, -1)
            } else {
                spanned_vector = c(spanned_vector , status_on_read[j])
            }
        }
        out_mat = rbind (out_mat, spanned_vector)

    }
    row.names (out_mat) = read_name_vec
    strand_map_df = data.frame("strand" = strand_vec)
    row.names (strand_map_df) = read_name_vec
    return (list (out_mat = out_mat, strand_map_df = strand_map_df))
}

