# This is the SARA2 data structure class

class DataStructure:
    def __init__(self, name, description, data_structure, value, upLimit):
        self.name = name
        self.description = description
        self.data_structure = data_structure
        self.value = value
        self.upLimit = upLimit

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
    def __init__(self):





