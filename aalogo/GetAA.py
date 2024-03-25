# standard libs
import ast
from configparser import ConfigParser
import numpy as np
from PIL import Image
# intern
import LogoUtil
import StandardConfig


def color_fader(position_number, c_top, c_bottom, color_advance=None):
    # fade (linear interpolate) from color c1 (at mix = 0) to c2 (mix = 1)
    """
    generates a color gradient between two given colors for each AA

    Paramters
    _________
    position_number : order of AA from 1 to 20
    c_top : top color of the image (first AA)
    c_bottom : bottom color of the image (last AA)

    Returns
    _______
    gradient color, needs to be iterated!
    """
    # color advance settings
    if color_advance is not None:
        value_color = color_advance[position_number]
    else:
        value_color = ((position_number+1)/20)
    c1 = np.array(c_top)
    c2 = np.array(c_bottom)
    return ((1-value_color) * c1 + (value_color * c2)).tolist()


def aa_image_colorizer(aa_config_section_name, font_type="bold_AA_fonts", config_set=True, color_grad=None,
                       order_aa=None, color_advance=None):
    """
    Amino Acid (aa) .png images reorder and recolor for AALogo generation

    Paramters
    _________
    aa_config_section_name : LogoStyle.ini entry of listed Amino Acids and their corresponding RGB-color
    font_type : font-folder package (bold_AA_fonts, classic_AA_fonts, modern_AA_fonts)

    Returns
    _______
    list_recolor_aa : list of [Amino Acid tag, recolored Amino Acid image]
    aa_compare : true order of Amino Acids post transformation (will be generated top to bottom)
    color_check_box_list : based on LogoStyle.ini, color categories of Amino Acids, can be customized freely
    """


    # get config file
    path_file, sep = StandardConfig.find_folderpath()

    config_file = "LogoStyle.ini"
    config = ConfigParser()
    config.read(config_file)
    config.sections()

    aa_compare = []  # were all AA mentioned?, if not append white AA
    list_recolor_aa = []

    # colors for gradient
    # all 20 canonical AA in order --> hydrophobicity to hydrophilicity (Kyte and Doolittle, 1982)
    color_top = [223, 130, 48]  # orange (hydrophobicity)
    color_bottom = [51, 154, 205]  # light blue (hydrophobicity)
    aa_matching_list = ["W", "F", "Y", "V", "L", "I", "M", "A", "P", "C",
                        "G", "S", "T", "N", "Q", "D", "E", "R", "K", "H"]
    color_check_box_list = []

    # folder paths for the asset generation:
    path_fonts = f"{path_file.split("aalogo")[0]}fonts"
    aa_font_type_path = f"{path_fonts}{sep}{font_type}"
    aa_index_box_path = f"{path_fonts}{sep}AA_letters_common"

    if config_set:  # if config is set to true

        # get all Amino Acid (AA) categories (define color categories of the AA-images)
        aa_category = list(config[aa_config_section_name])

        for category in aa_category:
            entry_category = ast.literal_eval(config[aa_config_section_name][category])
            aa_list = entry_category[0].split(",")
            for items in aa_list:
                if not (items in aa_matching_list):
                    print(f"{items} is not a valid character for amino acid nomenclature")
                    aa_list.remove(items)

            aa_compare.extend(aa_list)
            r, g, b = entry_category[1]  # assign new RGB values form .ini file entry
            for aa in aa_list:
                im_recolor = LogoUtil.convert_image_color(aa_font_type_path, aa, (r, g, b))
                list_recolor_aa.append([aa, im_recolor])

            # color boxes for index_box
            im_box_recolor = LogoUtil.convert_image_color(aa_index_box_path, "color_box_index", (r, g, b))
            color_check_box_list.append([category, im_box_recolor])

        list_non_specified_aa = [AA for AA in aa_matching_list if AA not in aa_compare]
        aa_compare.extend(list_non_specified_aa)
        for aa in list_non_specified_aa:
            im = Image.open(f"{aa_font_type_path}{sep}{aa}.png")
            list_recolor_aa.append([aa, im])

    else:
        if color_grad is not None:
            color_top = color_grad[0]
            color_bottom = color_grad[1]
        if order_aa is not None:
            aa_list = list(dict.fromkeys(order_aa))
            aa_compare = [str(amino).upper() for amino in aa_list if str(amino).upper() in aa_matching_list]

        list_non_specified_aa = [AA for AA in aa_matching_list if AA not in aa_compare]
        aa_compare.extend(list_non_specified_aa)

        # normalize color advance
        if color_advance is not None:
            color_advance = [(1-(float(i) - min(color_advance)) / (max(color_advance) - min(color_advance))) for i in
                             color_advance]
        for aa in aa_compare:
            im_recolor = LogoUtil.convert_image_color(aa_font_type_path, aa, color_fader(aa_compare.index(aa),
                                                                                         color_top, color_bottom,
                                                                                         color_advance))
            list_recolor_aa.append([aa, im_recolor])
        color_check_box_list = None             # color_check_box set false since it makes no sense as gradient

    return list_recolor_aa, aa_compare, color_check_box_list
