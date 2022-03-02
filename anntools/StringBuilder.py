from io import StringIO

class StringBuilder:

    def __init__(self):
        self.string = StringIO()

    def append(self, str):
        self.string.write(str)
        return self

    def __str__(self):
        return self.string.getvalue()

    def to_string(self):
        return self.string.getvalue()