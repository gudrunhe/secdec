'''
Move a file with the following modifications:
 o remove leading and trailing whitespaces, tabs,
   and backspaces (line continuation characters)
   from each line
 o remove existing newline characters
 o replace "#@SecDecInternalNewline@#" by
   newline characters
 o replace ";#@no_split_expression@# +=" by
   "+"
 o replace "#@SecDecInternalDblquote@#" by
   the double-quote character '"'
 o replace "#@SecDecInternalEscapedDblquote@#"
   by an escaped double-quote '\"'

This is necessary because FORMs formatting
is incompatible with c++.
The source is the first, the destination
is the second command line argument.

'''

from sys import argv
from os import remove, path

dest_dirname = argv[-1]
src_filenames = argv[1:-1]

for src_filename in src_filenames:
    dest_filepath = path.join(dest_dirname, src_filename)
    txt = []
    with open(src_filename, 'r') as src:
        for line in src:
            txt.append(line.strip('\\\t\n '))
    txt = ''.join(txt).replace("#@SecDecInternalNewline@#",'\n')
    txt = txt.replace('#@SecDecInternalSpace@#', ' ')
    txt = txt.replace(';#@no_split_expression@# +=', '+')
    txt = txt.replace('#@SecDecInternalDblquote@#', '"')
    txt = txt.replace('#@SecDecInternalEscapedDblquote@#', r'\"')
    with open(dest_filepath, 'w') as dest:
        dest.write(txt)

    remove(src_filename)
