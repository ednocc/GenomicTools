from . import json

class GenomeEncoder(json.JSONEncoder):
    def default(self):
        pass

class GenomeDecoder(json.JSONDecoder):
    def object_hook(self):
        pass
