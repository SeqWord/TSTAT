args = commandArgs(trailingOnly=TRUE)

if (length(args) == 5) {
  source_path = args[1]
  out_path    = args[2]
  ref_file    = args[3]
  info_file   = args[4]
  moltypes    = args[5]
} else {
  stop("Incorrect number of arguments was supplied (expecting 5 arguments).", call. = FALSE)
}

library(Rsubread)

bam_files = readLines(info_file)
for (fname in bam_files) {
  # Determine if file is paired-end or single-end based on file extension
  if (grepl("\\.pe\\.bam$", fname)) {
    is_paired = TRUE
  } else {
    is_paired = FALSE
  }

  count <- featureCounts(
    files = file.path(source_path, fname),
    annot.ext = ref_file,
    GTF.attrType = 'ID',
    isGTFAnnotationFile = TRUE,
    GTF.featureType = moltypes,
    isPairedEnd = is_paired
  )

  write.table(count$counts, file = file.path(out_path, paste0(fname, ".count.tmp")))
}
