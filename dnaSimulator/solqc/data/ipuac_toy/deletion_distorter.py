
import numpy as np
from data.ipuac_toy.distorter import Distorter


class DeletionDistorter(Distorter):
    def __init__(self, letter="A", rate=100, start_from=-1):
        self.letter = letter
        self.rate = rate
        self.start_from = start_from

    def distort(self, sequences):
        letters_counter = self.count_letters(sequences)
        l_count = letters_counter[self.letter]

        num_to_del = l_count // self.rate
        print(num_to_del)

        deleted = 0

        while deleted != num_to_del:
            for index, seq in enumerate(sequences):
                np_seq = np.array(list(seq))
                l_indices = np.where(np_seq == self.letter)[0]
                l_indices = l_indices[l_indices > self.start_from]

                # We sample one of the indices where we have our letter.
                sampled_index = np.random.choice(l_indices)

                # split the sequence to 'delete' the speicifc index
                seq = self.split_string_by_index(seq, sampled_index)

                # replace the old index with the new one.
                sequences[index] = seq

                deleted += 1

                # Once we reach the desired number of deletion break the loop.
                if deleted == num_to_del:
                    break


    @staticmethod
    def split_string_by_index(s, index):
        # Notice this also removes the value at the index.
        return s[:index] + s[index + 1:]


if __name__ == '__main__':
    # sequnces = generate_fasta_from_iupac(list("AAANAAA"))
    # ds = DeletionDistorter('A', 3)
    # ds.distort(sequnces)
    # print(sequnces)
    pass