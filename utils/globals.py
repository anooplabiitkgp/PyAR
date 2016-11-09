'''
globals.py - definition of some constants in PyAR
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

ELEMENTS = { 'X ': 0.00, \
        'h' : 0.32, 'he': 0.46, 'li': 1.20, 'be': 0.94, 'b' : 0.77, \
        'c' : 0.75, 'n' : 0.71, 'o' : 0.63, 'f' : 0.64, 'ne': 0.67, \
        'na': 1.40, 'mg': 1.25, 'al': 1.13, 'si': 1.04, 'p' : 1.10, \
        's' : 1.02, 'cl': 0.99, 'ar': 0.96, 'k' : 1.76, 'ca': 1.54, \
        'sc': 1.33, 'ti': 1.22, 'v' : 1.21, 'cr': 1.10, 'mn': 1.07, \
        'fe': 1.04, 'co': 1.00, 'ni': 0.99, 'cu': 1.01, 'zn': 1.09, \
        'ga': 1.12, 'ge': 1.09, 'as': 1.15, 'se': 1.10, 'br': 1.14, \
        'kr': 1.17, 'rb': 1.89, 'sr': 1.67, 'y' : 1.47, 'zr': 1.39, \
        'nb': 1.32, 'mo': 1.24, 'tc': 1.15, 'ru': 1.13, 'rh': 1.13, \
        'pd': 1.08, 'ag': 1.15, 'cd': 1.23, 'in': 1.28, 'sn': 1.26, \
        'sb': 1.26, 'te': 1.23, 'i' : 1.32, 'xe': 1.31, 'cs': 2.09, \
        'ba': 1.76, 'la': 1.62, 'ce': 1.47, 'pr': 1.58, 'nd': 1.57, \
        'pm': 1.56, 'sm': 1.55, 'eu': 1.51, 'gd': 1.52, 'tb': 1.51, \
        'dy': 1.50, 'ho': 1.49, 'er': 1.49, 'tm': 1.48, 'yb': 1.53, \
        'lu': 1.46, 'hf': 1.37, 'ta': 1.31, 'w' : 1.23, 're': 1.18, \
        'os': 1.16, 'ir': 1.11, 'pt': 1.12, 'au': 1.13, 'hg': 1.32, \
        'tl': 1.30, 'pb': 1.30, 'bi': 1.36, 'po': 1.31, 'at': 1.38, \
        'rn': 1.42, 'fr': 2.01, 'ra': 1.81, 'ac': 1.67, 'th': 1.58, \
        'pa': 1.52, 'u' : 1.53, 'np': 1.54, 'pu': 1.55 }
AUTOANG  = 0.52917726
AUTOKJ   = 2625.5
AUTOKCAL = 627.509541
