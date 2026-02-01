retrieve_json <- function (track_url = "", assembly = "dm6",
                           chromosome = "chr2L", st = "480290",
                           en = "480320", label = "peak229",
                           tr = "fp_and_mvec", target_dir = ""){
    
    download_command = paste0(track_url, ";", "genome=", assembly,";",
                              "track=", tr, ";",
                              "chrom=", chromosome,";",
                              "start=", st,";", "end=", en)

    tryCatch (
        {
            download.file(download_command, 
                      paste (target_dir, paste0(label, ".json"), sep = "/"))
        }, 
        error = function(e) {
            print (paste0("Download failed with error ", e))
            return (NA)
        },
        warning = function(w) {
            print (paste0("Please see the warning ", w))
            return (NA)
        },
        finally = {
            
            json_path = paste (target_dir, paste0(label, ".json"),sep = "/")
            return (json_path)
        }
    )
}

downloadJson <- function (track_url = "", assembly = "dm6",
                          chromosome = "chr2L", st = "",
                          en = "", label = "peak229",
                          tr = "fp_and_mvec", target_dir = ""){
    cat ("Using UCSC API to get JSON dump from BigBed Track ...\n")
    #print (c(track_url, assembly, chromosome, st, en, label))
    json_name <- retrieve_json (track_url = track_url, assembly = assembly,
                                chromosome = chromosome, st = st, en = en,
                                label = label, tr = tr, target_dir = target_dir)
    return (json_name)
    
}

JSON_to_Bed <- function (json_file = "peak229.json", label = "peak229", 
                         tr = "fp_and_mvec", target_dir = "") {
    x <- fromJSON(json_file)
    flat_names <- names (x[[tr]])
    df_of_json = as.data.frame(x[[tr]], col.names = flat_names)
    
    output_file = file (paste(target_dir, paste0(label, ".all_fp.bed"), 
                              sep = "/"), "w")
    output_file_name = paste(target_dir, paste0(label, ".all_fp.bed"),
                             sep = "/")
    for (fp in seq_len(nrow(df_of_json))) {
        chrom <- df_of_json$chrom[fp]
        st <- df_of_json$chromStart[fp]
        en <- df_of_json$chromEnd[fp]
        fp_info <- df_of_json$name[fp]
        score <- df_of_json$score[fp]
        strand <- "+"
        m_vec <- df_of_json$field8[fp]
        cat(paste(chrom, st, en, paste(fp_info, m_vec, sep = "|"),
                  score, strand, sep = "\t"), file = output_file, append = TRUE)
        cat ("\n", file  = output_file, append = TRUE)
    }
    close(output_file)
    return (output_file_name)
}

JSON_to_Bed_for_Nanopore <- function (json_file = "smac_seq.json",
                                      label = "smac_seq",
                                      tr = "nanopore_meth_calls",
                                      target_dir = "") {
    x <- fromJSON(json_file)
    flat_names <- names (x[[tr]])
    df_of_json = as.data.frame(x[[tr]], col.names = flat_names)
    
    output_file = file (paste(target_dir, paste0(label, ".methylation_calls.bed"),
                              sep = "/"), "w")
    output_file_name = paste(target_dir, paste0(label, ".methylation_calls.bed"),
                             sep = "/")
    for (meth_call in seq_len(nrow(df_of_json))) {
        chrom <- df_of_json$chrom[meth_call]
        st <- df_of_json$chromStart[meth_call]
        en <- df_of_json$chromEnd[meth_call]
        read_id <- df_of_json$name[meth_call]
        score <- df_of_json$score[meth_call]
        strand <- df_of_json$strand[meth_call]
        m_vec <- df_of_json$field8[meth_call]
        
        cat(paste(chrom, st, en, paste(read_id, m_vec, sep = "|"),
                  score, strand, sep = "\t"),
            file = output_file, append = TRUE)
        cat ("\n", file  = output_file, append = TRUE)
    }
    close(output_file)
    return (output_file_name)
}



intersectQueryandSubject <- function(query_file = "peak229.bed",
                                     subject_file = "peak229.all_fp.bed",
                                     overlap = "peak229_overlap.bed") {
    q_obj <- import (query_file)
    s_obj <- import (subject_file)
    overlaps <- findOverlaps(q_obj, s_obj, minoverlap = q_obj@ranges@width,
                             select = "all")
    q_obj_to_write <- resize (q_obj[overlaps@from, ],
                              width(q_obj[overlaps@from, ]) + 1, fix = "end")
    q_df = as.data.frame(q_obj_to_write)
    sub_q_df <- q_df [, seq(3)]
    
    s_obj_to_write <- resize (s_obj[overlaps@to, ],
                              width (s_obj[overlaps@to, ]) + 1, fix = "end")
    s_df <- as.data.frame(s_obj_to_write)
    result_df <- cbind (sub_q_df, s_df)
    return (result_df)
    # now reading files again to prepare
    
}
prepareQueryFile <- function (chromosome = "chr2L",
                              start = 480290,
                              end = 480320,
                              label = "peak229", target_dir = ""){
    query_file <- paste(target_dir, paste0(label, ".bed"), sep = "/")
    cat (paste(chromosome, start, end, sep = "\t"),
         file = query_file, append = FALSE)
    cat ("\n", file = query_file, append = TRUE)
    return (query_file)
}

overlapUsingBigBed <- function (
        bigbed_file = "inst/extdata/demo.bb",
        query_file = "peak229.bed"){
    q_obj <- import (query_file)
    s_obj <- import(bigbed_file, which = q_obj)
    overlaps <- findOverlaps(q_obj, s_obj, minoverlap = q_obj@ranges@width,
                             select = "all")
    q_obj_to_write <- resize (q_obj[overlaps@from, ],
                              width(q_obj[overlaps@from, ]) + 1, fix = "end")
    q_df <- as.data.frame(q_obj_to_write)
    sub_q_df <- q_df [, seq(3)]
    
    s_obj_to_write <- resize (s_obj[overlaps@to, ],
                              width (s_obj[overlaps@to, ]) + 1, fix = "end")
    s_df <- as.data.frame(s_obj_to_write)
    result_df <- cbind (sub_q_df, s_df)
    return (result_df)
    
}
