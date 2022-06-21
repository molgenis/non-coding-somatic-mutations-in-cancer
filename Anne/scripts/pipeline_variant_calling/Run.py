#!/usr/bin/env python3


class Run:
    """

    """

    def __init__(self, sample_num, type_sample, name, category):
        """
        Constructor
        :param sample_num:      Sample number of the participant
        :param type_sample:     Sample type (GL, FL, tFL)
        :param name:            Specific number of the tissue
        :param category:        Is the tissue tumor tissue or healthy tissue
        """
        self.sample_num = sample_num
        self.type_sample = type_sample
        self.name = name
        self.category = category
