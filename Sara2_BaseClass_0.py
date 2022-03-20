# This is the SARA2 data structure class
# this is my first update push

class DataStructure:
    def __init__(self, name, description, data_structure, value):
        self.name = name
        self.description = description
        self.data_structure = data_structure
        self.value = value

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


class DataCollection(DataStructure):

    def __init__(self, name):
        self.name = name
        self.dict = {}

    def init_dict(self):
        self.dict.clear()

    # vdd entry to dictionary
    def add_entry(self, key_name, value):
        taco = DataStructure
        taco.set_value(self, value)
        taco.set_data_structure(self, "work")
        temp_dict = {key_name, taco}
        self.dict.update(key_name=taco)
        tst=1
