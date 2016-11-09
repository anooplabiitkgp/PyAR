'''
gamma.py - to make gamma[sys] and gamma[rand] in PyAR
'''
'''
Copyright (C) 2016 by Surajit Nandi, Anoop Ayyappan, and Mark P. Waller
Indian Institute of Technology Kharagpur, India and Westfaelische Wilhelms
Universitaet Muenster, Germany

This file is part of the PyAR project.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
'''

import random

def get_systematic_gamma(lower_limit, upper_limit, step):
    gamma = []
    while lower_limit <= upper_limit:
        gamma.append(lower_limit)
        lower_limit += step
    return gamma


def get_random_low_gamma(gamma_max):
    allgamma = []
    sN = round(random.uniform(0.0, 1.0), 3)
    gamma = gamma_max * sN
    while gamma <= gamma_max:
        allgamma.append(gamma)
        gamma += 0.1*gamma_max
    return allgamma


def get_gamma(lower_limit, upper_limit):
    """ This function will generate a list of gamma values
        according to the new formula.
    """
    gamma = []
    #gamma_random_up_limit = random.randint(lower_limit, upper_limit)
    gamma_random_up_limit = upper_limit
    no_of_gamma_block = random.randint(1, 10)
    block_width = int(gamma_random_up_limit/no_of_gamma_block)
    for block_index in range(no_of_gamma_block):
        sN = round(random.uniform(0.0, 1.0), 3)
        igamma = (block_index * block_width) + (block_width * sN) + lower_limit
        gamma.append(igamma)
    return gamma


def main():
    pass


if __name__ == "__main__":
    main()
