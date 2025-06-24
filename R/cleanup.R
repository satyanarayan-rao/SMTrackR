deleteIntermediates <- function (label = "peak229"){
    file.remove(paste0(label, "_assigned_states.tsv"))
    file.remove(paste0(label, "_labelled_and_ordered.tsv"))
    file.remove(paste0(label, "_verbose.tsv"))
    if (file.exists(paste0(label, "all_fp.bed"))){
        file.remove(paste0(label, ".all_fp.bed"))
    }
    file.remove(paste0(label, ".bed"))
    if (file.exists(paste0(label, ".json"))){
        file.remove(paste0(label, ".json"))
    }

    file.remove(paste0(label, ".num.fp.tsv"))
    file.remove(paste0(label, ".num.mvec.tsv"))
}
