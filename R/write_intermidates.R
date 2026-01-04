writeIntermediateFiles <- function (ordered_by_length = NULL,
                                    label = "peak229", target_dir = "") {
    write.table(ordered_by_length,
                file = paste(target_dir, 
                             paste0(label, "_labelled_and_ordered.tsv"), 
                             sep = "/"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    num_coverted_fp <- unlist (lapply(ordered_by_length$fp_vec,
                                      function (x) {convert_fp_char_to_num(x)}))
    num_coverted_m <- unlist (lapply(ordered_by_length$m_vec,
                                     function (x) {convert_m_char_to_num(x)}))
    num_fp_df <- data.frame(read_id = ordered_by_length$read_id,
                            num_coverted_fp = num_coverted_fp)
    num_mvec_df <- data.frame(read_id = ordered_by_length$read_id,
                              num_coverted_m = num_coverted_m)
    write.table (num_fp_df, 
                 file = paste(target_dir, paste0(label, ".num.fp.tsv"), sep ="/"),
                 sep = "\t", row.names = FALSE,
                 col.names = FALSE, quote = FALSE)
    write.table (num_mvec_df, 
                 file = paste(target_dir, paste0(label, ".num.mvec.tsv"), sep ="/"),
                 sep = "\t", row.names = FALSE,
                 col.names = FALSE, quote = FALSE)
}
