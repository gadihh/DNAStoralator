import numpy as np
from data.ipuac_toy.deletion_distorter import DeletionDistorter
from data.ipuac_toy.insertion_distorter import InsertionDistorter
from data.ipuac_toy.mismatch_distorter import MismatchDistorter

IUPAC = {
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'R': 'AG',
        'Y': 'CT',
        'S': 'GC',
        'W': 'AT',
        'K': 'GT',
        'M': 'AC',
        'B': 'CGT',
        'D': 'AGT',
        'H': 'ACT',
        'V': 'ACG',
        'N': 'ACGT'
    }


def char_array_to_string(c_array):
    return ''.join(c_array)

def generate_fasta_from_iupac(word):
    variants = []
    def extract_code(word, i):
        if i == len(word):
            variants.append(char_array_to_string(word))
            return

        for letter in IUPAC[word[i]]:
            word[i] = letter
            extract_code(word.copy(), i + 1)

    extract_code(word, 0)
    return variants


sequences = generate_fasta_from_iupac(list("YYYYYACGACATCGTACAGGTCACTGTACGTACGTTAAGAGATGAGA"))
sequence_template = """@MN00713:7:000H2JMTL:1:11101:8035:1104 1:N:0:2
{}
+
AFFFAFFFAFFFFFFFFF#FFFFFFFFFFF#FFFFFFFFFFFFFFFFFFFF/FFFF#AFFF/IIIIIIIIIIIIIIIIIIIIIIIIIF/IIIIIIII/IIIIFIIIIIIIIIIIIFII/IIIIIIII/IIIIIIIIIIFI=IIIIIFII/I/F////##FFFFFFF/FFF6FFF//F/A/F/#F6#FFFFF/F//6A////F#AFF///#A"""

# print(sequence.format(var[0]))

with open('generated.fastq', 'w') as file:
    samples_a = np.random.choice(sequences, 100000 , replace=True)

    barcode_end = 5
    # disorter = DeletionDistorter("A", 1000, start_from=barcode_end)
    # disorter.distort(samples_a)

    disorter = DeletionDistorter("G", 2000, start_from=barcode_end)
    # disorter.distort(samples_a)

    disorter = MismatchDistorter("T", "AC", 1000, start_from=barcode_end)
    # disorter.distort(samples_a)

    disorter = InsertionDistorter("G", 1000, start_from=barcode_end)
    disorter.distort(samples_a)


    for sequence in samples_a:
        file.write(sequence_template.format(sequence))
        file.write("\n")
