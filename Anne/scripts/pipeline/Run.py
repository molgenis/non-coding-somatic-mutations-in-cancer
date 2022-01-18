#!/usr/bin/env python3


class Run:
    """

    """

    def __init__(self, sample_num, type_sample, name, category):
        """
        Constructor
        :param sample_num:      sample number of the participant
        :param type_sample:     sample type (GL, FL, tFL)
        :param name:            specific number of the tissue
        :param category:        is the tissue tumor tissue or healthy tissue
        """
        self.sample_num = sample_num
        self.type_sample = type_sample
        self.name = name
        self.category = category
