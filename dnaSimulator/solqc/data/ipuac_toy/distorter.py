"""
Distorter Class
---------------
This is an abstract class for the distorter family.
A distorter class is intended to "distort" generated data.
"""
from collections import Counter


class Distorter():
    def distort(self):
        print("This is an abstract method which your distorted should implement.")

    @staticmethod
    def get_events_number(count, rate):
        '''
        Given a count and a rate return the number of events that should happen
        to produce the given rate.
        :param count:
        :param rate:
        :return: The number of events that should happen to produce the given rate.
        '''
        return count // rate

    @staticmethod
    def count_letters(sequences):
        '''
        Given a list of sequences returns a dictionary counter with
        key = letter
        value = the number of occurrences of the letter in all the sequecnes.
        :param sequences: The sequences we want to build the counter for.
        :return: A dictionary counter with key = letter, value = # of occurrences of the letter in the sequences.
        '''
        c = Counter()
        for seq in sequences:
            c.update(seq)

        return c