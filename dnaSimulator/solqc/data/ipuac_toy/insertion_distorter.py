"""
Mismatch Distorter
------------------
This class adds insertions errors in generated data.
"""

from data.ipuac_toy.distorter import Distorter
import numpy as np


class InsertionDistorter(Distorter):
    def __init__(self, inserted_letter, rate, start_from=0):
        self.inserted_letter = inserted_letter
        self.rate = rate
        self.start_form = start_from

    def distort(self, sequences):
        letters_counter = self.count_letters(sequences)
        from_letter_count = letters_counter[self.inserted_letter]

        insertion_to_produce = self.get_events_number(from_letter_count, self.rate)

        while insertion_to_produce != 0:
            for index, seq in enumerate(sequences):
                # Get random position for insertion
                insertion_index = np.random.randint(self.start_form, len(seq)-1)

                # Create a new sequence with the inserted letter in the sampled index
                new_seq = seq[:insertion_index] + self.inserted_letter + seq[insertion_index:]

                # replace the old sequence with the new one.
                sequences[index] = new_seq

                insertion_to_produce -= 1

                # Once we reach the desired number of deletion break the loop.
                if insertion_to_produce == 0:
                    break



