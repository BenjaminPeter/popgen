#!/usr/bin/env python
"""
    a population designates  population sample
"""

import numpy as np



class Population:
    """attributes:
        location: physical location (coordinates) of the population
        id: identifier
        name: name of pop
        sample_size: number of sampled haplotypes
    """

    pop_counter = 0 
    def __init__(self, **kwargs):
        self.id = Population.pop_counter
        Population.pop_counter += 1
        if 'file' in kwargs:
            self = Population.load(kwargs['file'])
            return
        elif 'line' in kwargs:
            self.load_line(kwargs['line'])
            return

        location = None
        self.name = "Anon"
        self.sample_size = 0
        print "run init"

    def load(self,f):
        """
            if file is a string, try opening and reading first line
            otherwise, assume it's a file handle
        """
        if type(f) is str:
            f = open(f,"r")
        line = f.readline()
        self.load_line(line)


    def load_line(self,line):
        line = line.split()
        self.name, self.sample_size = line[0], int(line[1])
        self.location = [float(f) for f in line[2:]]



