import numpy as np
from PIL import Image
import StandardConfig


def convert_image_color(folder_path, image_name, tuple_rgb):
    """
    Transforms white pixels into colored pixels based on a rgb-tuple

    Parameters
    __________
    folder_path : global directory of folder
    image_name : .png file name that for recoloring (white surfaces only, no alpha)
    tuple_rgb : Tuple of RGB values --> (Red, Green, Blue)
                values ranging from 0 (no color) to 255 (full color)

    Returns
    _______
    im_recolor : picture with recolored white surfaces
    """

    sep = StandardConfig.find_folderpath()[1]
    im = Image.open(f"{folder_path}{sep}{image_name}.png")
    im_rgba = im.convert("RGBA")
    color_channels = np.array(im_rgba)
    red, green, blue, alpha = color_channels.T  # alpha is not used
    white_area = (red == 255) & (blue == 255) & (green == 255)

    color_channels[..., :-1][white_area.T] = tuple_rgb  # new RGB-values
    im_recolor = Image.fromarray(color_channels)

    return im_recolor


def check_input_string(input_str):
    """
    Checks if input string has allowed characters

    Parameters
    __________
    input_str : Amino Acid sequence

    Returns
    _______
    input_str_modify : Amino Acid sequence replaces U with S, all letters upper case
    """
    if not isinstance(input_str, str):
        print(f"Input Sequence is not str type")
        exit()

    input_str_modify = input_str.upper().replace("U", "S")  # replace U (Selenocysteine) with S (Cysteine)
    aa_matching_list = ["W", "F", "Y", "V", "L", "I", "M", "A", "P", "C",
                        "G", "S", "T", "N", "Q", "D", "E", "R", "K", "H"]

    for letters in input_str_modify:
        if not (letters in aa_matching_list):
            print(f"{letters} is not allowed, canonical amino acids are: {aa_matching_list}")
            break

    return input_str_modify


def rgb_to_hex(r, g, b):
    """
    Parameters
    __________
    r, g, b : red, green, blue values ranging from 0 (black) to 255 (white)

    Returns
    _______
    hex-converted color code e.g. #FFFFFF
    """
    return '#{:02x}{:02x}{:02x}'.format(r, g, b)


def hex_to_rgb(hex_code):
    """
    Paramters
    _________
    hex : color format e.g. #FFFFFF

    Returns
    _______
    rgb-converted color code as tuple (r, g, b) : red, green, blue values ranging from 0 (black) to 255 (white)
    """
    rgb = []
    for i in (0, 2, 4):
        decimal = int(hex_code[i:i + 2], 16)
        rgb.append(decimal)

    return tuple(rgb)
