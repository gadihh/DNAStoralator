"""
This class represents a pdf content.
Types : Image, Text (Maybe more in the future).
It contains a type and data.

If content type is text than data is a string.
If content type is image than data is the image name. 
"""

# TODO can be an abstract class and have the image and text be sub classes of it.
class Content (object):
    from enum import Enum

    class Type(Enum):
        TEXT = 1
        IMAGE = 2
        NEW_PAGE = 3

    def __init__(self, type, data):
        self.type = type
        self.data = data

    def __call__(self):
        return self.data

