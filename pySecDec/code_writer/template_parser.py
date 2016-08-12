"""
Functions to generate c++ sources from template files.

"""

import os

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

    # apply replacements
    string = string % replacements

    # write parsed file
    with open(dest, 'w') as dest_file:
        dest_file.write(string)

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
