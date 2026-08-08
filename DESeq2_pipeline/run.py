import sys, os, argparse
path = os.getcwd()
sys.path.append(os.path.join(path, "lib"))
import main

__version__ = "1.0"
date_of_creation = "June 7, 2025"

help_text = ("\n" + f"""
{'='*80}
■ TSTAT v{__version__} ■
{'='*80}
📅 Created: {date_of_creation}

🔬 Purpose:
Pipeline to process FASTQ transcripts.

💻 Tested Environments:
- CentOS Linux 7.3.1611
- Ubuntu 20.04 LTS
- Windows 10.0.22631.5039

⚙️ Dependencies:
┌──────────────────────────────┬───────────────────────┐
│ Tool                         │ Minimum Version       │
├──────────────────────────────┼───────────────────────┤
│ Python                       │ 3.9                   │
│ R                            │ R-4.4.3               │
│ Bioconductor                 │ 3.20                  │
│ Bioconductor::DESeq2         │ 1.46.0                │
│ Bioconductor::Rsubread       │ 2.20.0                │
│ Bioconductor::GenomicFeatures│ 1.58.0                │
└──────────────────────────────┴───────────────────────┘

R module installation:

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    # Now install any Bioconductor package, e.g., DESeq2
    BiocManager::install("DESeq2")
    
    # Now install any Bioconductor package, e.g., Rsubread
    BiocManager::install("Rsubread")
    
    
    # Now install any Bioconductor package, e.g., GenomicFeatures
    BiocManager::install("GenomicFeatures")


🚀 Usage:
    python3 run.py - you may be prompted to select R version:
        Select R version:
        1       x64
        Select R version number: 1
        
    Next, an interactive command-prompt menu will be shown:
    
        Settings for this run:
        
          I    Project folder           :
          F    Subfolder (optional)     :
          R    Reference GBK file       :
          M    Moltypes to include      : CDS
          T    Single or paired end?    : single
          C    Control sample marker    :
          E    Experiment sample marker :
          B    Project basename         :
          P    Plot results             : No
          S    Start point              : 0.scratch
        
        Press L-Enter to load the last run options.
        Y to accept these settings, type the letter for one to change or Q to quit
    
Alternatively, set the arguments in the command line:
 
    python3 run.py <arguments>
    
⚡ Required Arguments:
    -i, --project           Project folder in 'input' folder (default: '')
    -r, --reference         Reference GBK file in input/project/embedded/ (default: '')
    -b, --basename          Locus tag base (default: '')
                            If /locus_tag="K8B78_00005"
                            set --basename 'K8B78'
    -c, --control           Marker of control read files, like '_ctrl_' (default: '')
    -e, --experimen         Marker of experimet read files, like '_exp_' (default: '')
    
🔧 Optional Arguments:
    -u, --input             Input folder (default: 'input')
    -o, --output            Output folder (default: 'output')
                            output/project an tmp/project folders will be created
    -f, --embedded          Subfolder within input/project folder (default: '')
    -m, --moltype           Comma-separated gene types to include (default: 'CDS')
    -t, --type              Read type single|paired-end|smart (default: 'single')
                            If paired-end, files must have '_1' and '_2' marks,
                            if type is 'smart', unpaired files are treated as 'single' 
    -p, --plot              Comma-separated graphical plots to create: 
                                Volcano plot,BaseMean plot,Gene expression (default: '')
    -s, --step              Step to start the program with (default: '0.scratch')
                            - 0.scratch - start from creation folders tmp/project and output/project
                            - 1.mapping - mapping reads against reference, 
                                store resulting *.bam and *.vcf files to output/project/alignments
                            - 2.counting - count reads against genomic features (CDS);
                                store resulting *.count and coldata files to output/project/counts
                            - 3.normalization - DESeq2 statistics;
                            - 4.plotting - creation output CSV graphs and store in output/project;
    
🆘 Help Options:
    -help, --help           Show this help message
    -version, --version     Show version information\n""" +
          
f"{'='*80}")


def show_banner():
    print("""
 ███████╗  ██████╗  ███████╗  █████╗  ███████╗
 ╚══██╔╝  ██╔════╝  ╚══██╔╝  ██╔══██╗ ╚══██╔╝
    ██║   ╚█████╗     ██║    ███████║    ██║
    ██║    ╚═══██╗    ██║    ██╔══██║    ██║
    ██║   ██████╔╝    ██║    ██║  ██║    ██║

 Pipeline to process FASTQ transcripts.
""")

def show_help():
    print(help_text)
    
def parse_arguments():

    # Define parser and disable automatic help
    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument("-u", "--input", default="input", help="Input folder")
    parser.add_argument("-o", "--output", default="output", help="Output folder")
    parser.add_argument("-i", "--project", default="", help="Project folder name")
    parser.add_argument("-f", "--embedded", default="", help="Embedded folder name")
    parser.add_argument("-r", "--reference", default="", help="Reference GBK file")
    parser.add_argument("-m", "--moltypes", nargs="+", default=["CDS"], help="Moltypes to include")
    parser.add_argument("-t", "--type", choices=["single", "paired-end", "smart"], default="single", help="Sequencing type")
    parser.add_argument("-b", "--basename", default="", help="Basename")
    parser.add_argument("-c", "--control", default="", help="Control sample marker")
    parser.add_argument("-e", "--experiment", default="", help="Experiment sample marker")
    parser.add_argument("-p", "--plots", default="", help="Graphical plots: Volcano plot,BaseMean plot,Gene expression")
    parser.add_argument("-s", "--step", default="0.scratch", help="Pipeline step to start")
    parser.add_argument("-help","--help","-h", action="store_true", help="Show this help message and exit")
    parser.add_argument("--version", "-v", action="version",
                    version=f"\nversion {__version__} created on {date_of_creation}\n")

    return parser.parse_args()

###############################################################################
if __name__ == "__main__":
    args = parse_arguments()
    
    # Show help and exit
    if args.help:
        show_help()
        exit()

    options = {
        "-u": args.input,
        "-o": args.output,
        "-i": args.project,
        "-f": args.embedded,
        "-r": args.reference,
        "-m": args.moltypes,
        "-t": args.type,
        "-b": args.basename,
        "-c": args.control,
        "-e": args.experiment,
        "-p": args.plots,
        "-s": args.step,
    }

    show_banner()

    oMain = main.Interface(options)
