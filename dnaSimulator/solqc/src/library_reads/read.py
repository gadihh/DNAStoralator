class Read():
    def __init__(self, row):
        self.data = row.to_dict()

        # Converting keys to lower case for maintenance
        self.data = {k.lower(): v for k, v in self.data.items()}

    def get_cigar_path(self):
        return self.get_attribute('cigar_path')

    def get_variant_id(self):
        return self.get_attribute('variant_id')

    def get_row_count(self):
        return self.get_attribute('count')

    def has_attribute(self, attribute):
        if attribute in self.data:
            return True
        return False

    def is_matched(self):
        return self.get_variant_id() != -1

    def get_attribute(self, attribute_name):
        attribute_name = attribute_name.lower()
        if self.has_attribute(attribute_name):
            return self.data[attribute_name]
        else:
            self.print_attribute_not_in_data(attribute_name)

    def __call__(self, *args, **kwargs):
        return self.get_attribute('sequence').upper()

    @staticmethod
    def print_attribute_not_in_data(attribute):
        print ("The variant do not have {} attribute".format(attribute))