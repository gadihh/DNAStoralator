from src.utils.biology import get_reverse_compliment


class FixPreProcessor(object):

    PREFIX = ""
    SUFFIX = ""

    def __init__(self, data):
        if data is None:
            print ("| No pre processing")
            return

        if 'prefix' in data and len(data['prefix']) > 0:
            self.PREFIX = data['prefix']
        else:
            self.PREFIX = False

        if 'suffix' in data and len(data['suffix']) > 0:
            self.SUFFIX = data['suffix']
        else:
            self.SUFFIX = False

        if 'length' in data and data['length'] > 0:
            self.LENGTH = data['length']
        else:
            self.LENGTH = False

    def process(self, sequence):
        original = self.check_original(sequence)
        if original:
            return original

        reversed_compliment = self.check_reverse_compliment(sequence)
        if reversed_compliment:
            return reversed_compliment

        return None

    def check_sequence(self, sequence):
        prefix_index, suffix_index = 0, len(sequence)
        if self.PREFIX:
            prefix_index = sequence.index(self.PREFIX) if self.PREFIX in sequence else None
            if prefix_index is None:
                return None

        if self.SUFFIX:
            suffix_index = sequence.index(self.SUFFIX) + len(self.SUFFIX) if self.SUFFIX in sequence else None
            if suffix_index is None:
                return None
        #print(prefix_index)
        #print(suffix_index)

        if prefix_index is not None and suffix_index is not None:
            #print(prefix_index+len(self.PREFIX))
            #print(suffix_index-len(self.SUFFIX))
            if self.SUFFIX and self.PREFIX:
                trimmed_sequence = sequence[prefix_index+len(self.PREFIX):suffix_index-len(self.SUFFIX)]
            else:
                trimmed_sequence = sequence[prefix_index:suffix_index]
            #print(trimmed_sequence)
            if not self.LENGTH or self.LENGTH - 5 <= len(trimmed_sequence) <= self.LENGTH + 5:
                return trimmed_sequence

        else:
            return None

    def check_original(self, sequence):
        return self.check_sequence(sequence)


    def check_reverse_compliment(self, sequence):
        reversed_compliment = get_reverse_compliment(sequence)
        return self.check_sequence(reversed_compliment)
