# This is the SARA2 data structure class
# this is my first update push

class DataStructure(object):
    def __init__(self):
        self.value = None
        self.name = None
        self.description = None

    def set_name(self, name):
        self.name = name

    def get_name(self):
        return self.name

    def set_description(self, description):
        self.description = description

    def get_description(self):
        return self.description

    def set_data_structure(self, data_structure):
        self.data_structure = data_structure

    def get_data_structure(self):
        return self.data_structure

    def set_value(self, value):
        self.value = value

    def get_value(self):
        return self.value


class DataCollection:

    def __init__(self, name):
        self.name = name
        self.dict = {}

    def init_dict(self):
        self.dict.clear()

    # vdd entry to dictionary
    def add_entry(self, key_name, value):
        temp_struct = DataStructure
        temp_struct.set_name(self,key_name)
        temp_struct.set_value(self, value)
        temp_dict = {key_name: temp_struct}
        self.dict.update(temp_dict)

    def get_entry(self, key_name):
        value = self.dict[key_name].get_name(self)
        return value




