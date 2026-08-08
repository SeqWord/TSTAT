args = commandArgs(trailingOnly=TRUE)

if (length(args)==7) {
  source_path = args[1]
  workdir_path = args[2]
  out_path = args[3]
  ref_file = args[4]
  info_file = file.path(workdir_path, args[5])
  reads = args[6]
  title = args[7]
} else {
  stop("Incorrect number of arguments was supplied (input file).", call.=FALSE)
}

library(Rsubread)
setwd(workdir_path)
buildindex(basename=title, reference=ref_file)
basename = file.path(workdir_path, title)

fastq = readLines(info_file)

for (name in fastq) {
  if (reads != "single") {
    fnames = unlist(strsplit(name, "[,]"))
    fname1 = file.path(source_path, fnames[1])
    fname2 = ifelse(fnames[2] != " ", file.path(source_path, fnames[2]), NULL)
    bam_suffix = ".pe.bam"
    file_base = tools::file_path_sans_ext(basename(fnames[1]))
  } else {
    fname1 = file.path(source_path, name)
    fname2 = NULL
    bam_suffix = ".se.bam"
    file_base = tools::file_path_sans_ext(basename(name))
  }

  out_bam = file.path(out_path, paste0(file_base, bam_suffix))
  align(index=basename, readfile1=fname1, readfile2=fname2, output_file=out_bam)
}
