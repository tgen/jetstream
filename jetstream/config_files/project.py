""" Project model and orm definitions """

class Project(object):
    """ Basic model for a project """
    def __init__(self, samples, properties, format=None, source=None):
        self.samples = samples
        self.properties = properties
        self.format = format
        self.source = source

