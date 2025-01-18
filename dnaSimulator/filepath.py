import os
import time

"""
Some helpers to deal with files, filepath and filenames manipulations.

Need to mention that all dataset files are structured into a pattern:
type_tech1_tech2_*.txt or type_*.txt
Where type is either a string 'error' or 'clustered', represents the general algorithm that created it.
tech1, tech2 words are short codes for the technologies that were used in creating the dataset. See tech_dict
"""

tech_dict = {'miseq_twist': 'Twist Bioscience + Ilumina miSeq',
             'miseq_custom': 'CustomArray + Ilumina miSeq',
             'nextseq_twist': 'Twist Bioscience + Ilumina NextSeq',
             'minion_idt': 'Integrated DNA Technology (IDT) + MinION',
             'stutter_def': 'Stutter',
             'other': 'Other'}


def rename_file(old_name, new_name):
    if os.path.exists(old_name):
        os.rename(old_name, new_name)


def get_base_name(path):
    return os.path.basename(path).split('.')[0]


def get_show_name(path):
    """
    Given a path, return a display string for the table.
    """
    splits = os.path.basename(path).split('_')
    start_idx = 1
    if len(splits) > 2:
        tech = '_'.join(splits[1:3])
        if tech in tech_dict:
            start_idx = 3
    show_name = '_'.join(splits[start_idx:]).split('.')[0]
    return show_name


def get_tech_from_path(path):
    """
    Given a path to a dataset, return a technology label using which it was created.
    """
    split = os.path.basename(path).split('_')
    if len(split) > 2:
        key = '_'.join(split[1:3])
        if key not in tech_dict:
            key = 'other'
    else:
        key = 'other'
    return tech_dict[key]


def get_dir_from_path(path):
    dir = os.path.dirname(path)
    return dir


def replace_suffix(old_base_name, new_suffix):
    splits = old_base_name.split('_')
    new_prefix = splits[0]
    if len(splits) > 2:
        tech = '_'.join(splits[1:3])
        if tech in tech_dict:
            new_prefix = '_'.join(splits[:3])
    new_base_name = new_prefix + '_' + new_suffix + '.txt'
    return new_base_name


def replace_suffix_path(old_filepath, new_suffix):
    dir = get_dir_from_path(old_filepath)
    old_base_name = get_base_name(old_filepath)
    new_base_name = replace_suffix(old_base_name, new_suffix)
    new_filepath = join_path(dir, new_base_name)
    # if new name is not unique - return empty string.
    return new_filepath if not os.path.exists(new_filepath) else ''


def delete(filepath):
    if os.path.exists(filepath):
        os.remove(filepath)


def sizeof_fmt(num):
    suffix = "B"
    for unit in ["", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"]:
        if abs(num) < 1024.0:
            return f"{num:3.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}Yi{suffix}"


def get_shuffled_path_from_orig(orig_path):
    # inserting 'shuffled' after 'error' prefix
    splitted = orig_path.split('_')
    splitted.insert(2, 'shuffled')
    return '_'.join(splitted)


def remove_extension(path):
    return path.split('.')[0]


def join_path(dir, filename):
    return os.path.join(dir + os.sep, filename).replace("\\", "/")


def files(path):
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            yield file


def time_stamp():
    return round(time.time())

