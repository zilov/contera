def download_blast(outdir):
    blast_nt_url = "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz"
    blastnt_data_dir = os.path.join(outdir, "blastnt")
    blast_db = os.path.join(blastnt_data_dir, "nt")
    if not os.path.exists(blastnt_data_dir):
        os.makedirs(blastnt_data_dir)
    elif os.path.exists(blast_db):
        print("\nBLAST database already downloaded, continue...\n")
        return blast_db
    cmd = (
        f'cd {blastnt_data_dir} && '
        f'echo Download blast database..'
        f'wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O nt.gz {blast_nt_url} && '
        f'echo Decompressing... && '
        f'gunzip nt.gz && '
        f'echo Build database... &&'
        f'makeblastdb -dbtype nucl -in nt'
    )
    os.system(cmd)
    return blast_db