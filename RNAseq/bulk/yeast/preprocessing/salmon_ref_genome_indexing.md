wget ftp://ftp.ensembl.org/pub/release-97/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
salmon index -t Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa -i yeast_salmon_index -k 31 
