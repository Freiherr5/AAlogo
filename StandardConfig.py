import os
import pathlib


def find_folderpath():
    """
    Identification of Operating System (OS) and current directory of files

    Returns
    _______
    path : current directory path the .py file is located on the OS
    notation : either Linux "/" symbol or Windows "\" symbol for path string
    """

    path = str(pathlib.Path().absolute())
    # generate directory annotation
    if bool(path.find("C:") == 0):  # Windows path
        notation = "\\"
    else:
        notation = "/"  # Linux path

    return path, notation


def make_directory(name_dir, path_dir=None):
    """
    For Directory creation utility

    Parameters
    __________
    name_dir : directory name
    path_dir : directory path --> if None, creation of name_dir in .py folder
    """

    if path_dir is None:
        path_dir = find_folderpath()[0]

    sep = find_folderpath()[1]  # either seperator for Linux or Windows

    # check and create new target directory if input given
    path_name_dir = f"{path_dir}{sep}{name_dir}"
    if os.path.exists(path_name_dir):
        print("Path already exists, skip folder creation...")
    else:
        os.mkdir(path_name_dir)
        print("Path " + str(path_name_dir) + " is created...")  # read files from target folder "input"
