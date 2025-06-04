sbatch -p short --mem=16G --cpus-per-task=8 --time=02:00:00 --wrap="salmon index -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i salmon_index_grch38_r97 -k 31 --threads 8"
