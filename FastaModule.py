## This is a module to read a fasta file into a dictionary

def FaAsDict(fasta):
        seqs = {}
        for line in fasta:
                line = line.rstrip("\n")
                if line.startswith(">"):
                        header = line.replace(">", "")
                        seqs[header] = ""
                else:
                        seqs[header] = seqs[header] + line
        return seqs

