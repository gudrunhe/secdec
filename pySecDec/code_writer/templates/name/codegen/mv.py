'''
Move a file with the following modifications:
 o remove leading and trailing whitespaces, tabs,
   and backspaces (line continuation characters)
   from each line
 o remove existing newline characters
 o replace "#@SecDecInternalNewline@#" by
   newline characters

This is necessary because FORMs formatting
is incompatible with c++.
The source is the first, the destination
is the second command line argument.

'''

from sys import argv
from os import remove

src_filename, dest_filename = argv[-2:]

with open(src_filename, 'r') as src:
    with open(dest_filename, 'w') as dest:
        for line in src:
            dest.write(line.strip('\\\t\n ').replace("#@SecDecInternalNewline@#",'\n'))

remove(src_filename)
