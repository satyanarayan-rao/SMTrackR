#hubs = read.table("hubs.tsv", sep = "\t", header = T,
#                  stringsAsFactors = F, check.names = F)
#hubs = gsheet2tbl("https://docs.google.com/spreadsheets/d/1jyM4LhhrCPKtrnj81MDizc33Za1yvar80-HWC47kyfQ/edit?gid=1076362671#gid=1076362671")
#row.names (hubs) = paste(hubs$Organism, hubs$model,
#                         hubs$Condition, hubs$Type, hubs$Genome, sep = "@")

retrieve_json <- function (track_url = "", assembly = "dm6",
                           chromosome = "chr2L", st = "480290",
                           en = "480320", label = "peak229",
                           tr = "fp_and_mvec"){

    download_command = paste0(track_url, ";", "genome=", assembly,";",
                              "track=", tr, ";",
                              "chrom=", chromosome,";",
                              "start=", st,";", "end=", en)
    print (download_command)
    download.file(download_command, paste0(label, ".json"))
    return (paste0(label, ".json"))
}

downloadJson <- function (track_url = "", assembly = "dm6",
                          chromosome = "chr2L", st = "",
                          en = "", label = "peak229",
                          tr = "fp_and_mvec"){
    cat ("Using UCSC API to get JSON dump from BigBed Track ...\n")
    #print (c(track_url, assembly, chromosome, st, en, label))
    json_name <- retrieve_json (track_url = track_url, assembly = assembly,
                               chromosome = chromosome, st = st, en = en,
                               label = label, tr = tr)
    return(json_name)

}

JSON_to_Bed <- function (json_file = "peak229.json", label = "peak229", tr = "fp_and_mvec") {
    x <- fromJSON(json_file)
    flat_names <- names (x[[tr]])
    df_of_json = as.data.frame(x[[tr]], col.names = flat_names)


    output_file = file (paste0(label, ".all_fp.bed"), "w")
    output_file_name = paste0(label, ".all_fp.bed")
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
                                      tr = "nanopore_meth_calls") {
    x <- fromJSON(json_file)
    flat_names <- names (x[[tr]])
    df_of_json = as.data.frame(x[[tr]], col.names = flat_names)

    output_file = file (paste0(label, ".methylation_calls.bed"), "w")
    output_file_name = paste0(label, ".methylation_calls.bed")
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
                              start = "480290",
                              stop = "480320",
                              label = "peak229"){
    query_file <- paste0(label, ".bed")
    cat (paste(chromosome, start, stop, sep = "\t"),
         file = paste0(label, ".bed"))
    cat ("\n", file = paste0(label, ".bed"), append = TRUE)
    return (query_file)
}

overlapUsingBigBed <- function (
        bigbed_file = "local_bigbed/demo.bb",
        query_file = "peak229.bed"){
    q_obj <- import (query_file)
    s_obj <- import(bigbed_file, which = q_obj)
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
    print (head (result_df))
    return (result_df)

}
