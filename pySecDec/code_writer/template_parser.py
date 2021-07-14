"""
Functions to generate c++ sources from template files.

"""

import os
import re

def validate_pylink_qmc_transforms(pylink_qmc_transforms):
    '''
    Check if `pylink_qmc_transforms` are valid options and remove duplicates

    :param pylink_qmc_transforms:
        list or None;
        Required qmc integral transforms, options are:

        * ``korobov<i>x<j>`` for 1 <= i,j <= 6
        * ``korobov<i>`` for 1 <= i <= 6 (same as ``korobov<i>x<i>``)
        * ``sidi<i>`` for 1 <= i <= 6

    :return:
        Sorted set of pylink_qmc_transforms
    '''
    pylink_qmc_transforms_available_options = set(
        ['korobov'+str(i)+'x'+str(j) for i in range(1,7) for j in range(1,7)] +
        ['sidi'+str(i) for i in range(1,7)]
    )
    if pylink_qmc_transforms != None:
        # korobov%i -> korobov%ix%i
        korobov_symmetric_transforms = ['korobov%i' % i for i in range(1,7)]
        pylink_qmc_transforms = [ x if x not in korobov_symmetric_transforms else 'korobov'+x[7:]+'x'+x[7:] for x in pylink_qmc_transforms]
        # remove duplicates
        pylink_qmc_transforms = set(pylink_qmc_transforms)
    else:
        pylink_qmc_transforms = set(['korobov3x3']) # Default
    assert pylink_qmc_transforms.issubset(pylink_qmc_transforms_available_options), \
        '"%s" found in `pylink_qmc_transforms` but not in `pylink_qmc_transforms_available_options`' % \
        pylink_qmc_transforms.difference(pylink_qmc_transforms_available_options)
    return sorted(pylink_qmc_transforms)

def generate_pylink_qmc_macro_dict(macro_function_name):
    '''

    Generate translation from transform short names 'korobov#x#' and 'sidi#' to C++ macros

    :param macro_function_name:
        string;
        Name of the macro function to consider

    :return:
        dict;
        A mapping between the transform short names and C++ macros

    '''
    pylink_qmc_translation = {}
    for i in range(1, 7):
        for j in range(1, 7):
            pylink_qmc_translation['korobov' + str(i) + 'x' + str(j)] = macro_function_name + '_KOROBOV_QMC(%i,%i)' % (i, j)
    for i in range(1, 7):
        pylink_qmc_translation['sidi' + str(i)] = macro_function_name + '_SIDI_QMC(%i)' % i
    return pylink_qmc_translation

def parse_template_file(src, dest, replacements={}):
    '''
    Copy a file from `src` to `dest` replacing
    ``%(...)`` instructions in the standard python
    way.

    .. warning::
        If the file specified in `dest` exists, it
        is overwritten without prompt.

    .. seealso::
        :func:`.parse_template_tree`

    :param src:
        str;
        The path to the template file.

    :param dest:
        str;
        The path to the destination file.

    :param replacements:
        dict;
        The replacements to be performed.
        The standard python replacement rules
        apply:

            >>> '%(var)s = %(value)i' % dict(
            ...     var = 'my_variable',
            ...     value = 5)
            'my_variable = 5'

    '''
    # read template file
    with open(src, 'r') as src_file:
        string = src_file.read()

    # write parsed file
    with open(dest, 'w') as dest_file:
        def write_large(file,bigtext):
            textsize = len(bigtext)
            chunksize = 1000*10**6
            chunks = textsize//chunksize+1
            for chunk in range(chunks):
                dest_file.write(bigtext[chunk*chunksize:(chunk+1)*chunksize])
                dest_file.flush()

        def recursive_write(file, text, unflushedlen = [0]):
            if type(text) == list:
                for part in text:
                    recursive_write(file,part,unflushedlen)
                return
            if type(text) is not str:
                text = str(text)
            # file.write(text)
            write_large(file,text)
            unflushedlen[0] += len(text)
            if unflushedlen[0] > 10**9:
                file.flush()
                unflushedlen[0] = 0

        dest_file_parts = re.split("%\(([^)]*)\)s",string)
        dest_file.write(dest_file_parts[0] % replacements)
        for n in range(1,len(dest_file_parts),2):
            bigtext = replacements[dest_file_parts[n]]
            recursive_write(dest_file, bigtext)
            dest_file.write(dest_file_parts[n+1] % replacements)

def parse_template_tree(src, dest, replacements_in_files={}, filesystem_replacements={}):
    '''
    Copy a directory tree from `src` to `dest` using
    :func:`.parse_template_file` for each file and
    replacing the filenames according to
    `filesystem_replacements`.

    .. seealso::
        :func:`.parse_template_file`

    :param src:
        str;
        The path to the template directory.

    :param dest:
        str;
        The path to the destination directory.

    :param replacements_in_files:
        dict;
        The replacements to be performed in the
        files. The standard python replacement
        rules apply:

            >>> '%(var)s = %(value)i' % dict(
            ...     var = 'my_variable',
            ...     value = 5)
            'my_variable = 5'

    :param filesystem_replacements:
        dict;
        Renaming rules for the destination files. and
        directories. If a file or directory name in
        the source tree `src` matches a key in this
        dictionary, it is renamed to the corresponding
        value. If the value is ``None``, the
        corresponding file is ignored.

    '''
    # walk through tree
    for dirpath, dirnames, filenames in os.walk(src):
        # create target directory
        this_source_directory = os.path.abspath(dirpath)
        target_relpath_without_renamings = os.path.relpath(this_source_directory, src)

        # rename/ignore if desired
        # do not mody top-level name
        if this_source_directory != os.path.abspath(src):
            omit_directory = False
            target_relpath_building_blocks = []
            remainder, last_block = os.path.split(target_relpath_without_renamings)
            while last_block != '':
                last_block = filesystem_replacements.get(last_block, last_block) # rename
                if last_block is None:
                    omit_directory = True
                    break
                target_relpath_building_blocks.append(last_block)
                remainder, last_block = os.path.split(remainder)
            if omit_directory:
                continue
            target_relpath_building_blocks.reverse()
            target_relpath = os.path.join(*target_relpath_building_blocks)

            this_target_directory = os.path.join(dest, target_relpath)

        else:
            this_target_directory = dest

        os.mkdir(this_target_directory)

        # parse files
        for source_filename in filenames:
            # rename/ignore if desired
            target_filename = filesystem_replacements.get(source_filename, source_filename)
            if target_filename is None:
                continue

            # parse the file using `parse_template_file`
            source_file_path = os.path.join(this_source_directory, source_filename)
            target_file_path = os.path.join(this_target_directory, target_filename)

            try:
                parse_template_file(source_file_path, target_file_path, replacements_in_files)
            except Exception as error:
                error.args = tuple([arg for arg in error.args] + ['while parsing "' + source_file_path + '"'])
                raise
