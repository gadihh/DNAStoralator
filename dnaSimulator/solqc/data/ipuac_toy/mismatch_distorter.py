"""
Mismatch Distorter
------------------
This class inserts mismatch errors in generated data.
"""
import numpy as np

from data.ipuac_toy.distorter import Distorter


class MismatchDistorter(Distorter):
    def __init__(self, from_letter, to_letter, rate, start_from=0):
        '''

        :param from_letter: The letter you want to change
        :param to_letter: A string of possible letter to replace from_letter
                Example : "G", "AC", "GCT" etc..
        :param rate: The rate of events (example - rate=100 meaning every 100
               occureneces of from_letter preform a missmatch)
        :param start_from:
        '''
        self.from_letter = from_letter
        self.to_letter = np.asarray(list(to_letter))
        self.rate = rate
        self.start_from = start_from

    def distort(self, sequences):
        letters_counter = self.count_letters(sequences)
        from_letter_count = letters_counter[self.from_letter]

        missmatch_to_produce = self.get_events_number(from_letter_count, self.rate)

        while missmatch_to_produce != 0:
            for index, seq in enumerate(sequences):
                np_seq = np.array(list(seq))
                l_indices = np.where(np_seq == self.from_letter)[0]
                l_indices = l_indices[l_indices > self.start_from]

                # We sample one of the indices where we have our letter.
                sampled_index = np.random.choice(l_indices)

                # Preform the missmatch.
                np.put(np_seq, sampled_index, np.random.choice(self.to_letter))

                # replace the old sequence with the altered one.
                sequences[index] = ''.join(np_seq)

                missmatch_to_produce -= 1

                # Once we reach the desired number of deletion break the loop.
                if missmatch_to_produce == 0:
                    break
