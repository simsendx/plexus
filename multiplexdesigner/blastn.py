import subprocess

def blastn(query, database, outfile, threads = 4):
    """
    Run blastn locally.

    https://biopython.org/docs/dev/Tutorial/chapter_blast.html#introduction
    """

    cmd = "blastn -task blastn -query " + query + " -db " + database + " -num_threads " + threads + " -out " + outfile

    subprocess.run(cmd, capture_output=True)