"""
This class represents a variant in our design.
You can query it to get any information you want about a specific variant.
"""
import pandas as pd

class Variant(object):
    def __init__(self, row):
        self.data = row.to_dict()

        # Converting keys to lower case case for maintenance
        self.data = {k.lower(): v for k, v in self.data.items()}

        # And array of statistics objects for the variant.
        self.statistics = []

    def get_name(self):
        return self.get_attribute('name')

    def get_sequence(self):
        return self.get_attribute('sequence')

    def get_index(self):
        return self.get_attribute('index')

    def get_barcode(self):
        return self.get_attribute('barcode')

    def has_attribute(self, attribute):
        if attribute in self.data:
            return True
        return False

    def get_attribute(self, attribute_name):
        attribute_name = attribute_name.lower()
        if self.has_attribute(attribute_name):
            return self.data[attribute_name]
        else:
            self.print_attribute_not_in_data(attribute_name)

    def __call__(self, *args, **kwargs):
        return self.get_attribute('sequence').upper()

    @staticmethod
    def print_attribute_not_in_data(self, attribute):
        print ("The variant do not have {} attribute".format(attribute))


def test_variant(row):
    variant = get_sample_variant()
    print ("Name is : {}.".format(variant.get_name()))
    print ("Id is : {}.".format(variant.get_index()))
    print ("Barcode is : {}".format(variant.get_barcode()))
    print ("Group is : {}".format(variant.get_attribute('group')))


def get_sample_variant():
    data = pd.read_csv("/Users/yoavorlev/Desktop/library/data/Twist_zohar_design.csv", nrows=1)
    row = data.iloc[0]
    variant = Variant(row)
    return variant

if __name__ == "__main__":
    print("\n\t==-- Testing variant.py-==")
    test_variant()
    print("\n\t==-- Done Testing variant.py-==")