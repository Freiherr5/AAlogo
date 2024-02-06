from PIL import Image
from configparser import ConfigParser
import StandardConfig
import ast
import LogoUtil


def aa_image_colorizer(aa_config_section_name, font_type="classic_AA_fonts"):
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
    path, sep = StandardConfig.find_folderpath()  # get the current directory and the directory-nomenclature (/ or \)

    config_file = "LogoStyle.ini"
    config = ConfigParser()
    config.read(config_file)
    config.sections()

    # get all Amino Acid (AA) categories (define color categories of the AA-images)
    aa_category = list(config[aa_config_section_name])

    aa_compare = []  # were all AA mentioned?, if not append white AA
    list_recolor_aa = []
    aa_matching_list = ["W", "F", "Y", "V", "L", "I", "M", "A", "P", "C", "G", "S", "T", "N", "Q", "D", "E", "R", "K", "H"]  # all 20 canonical AA in order --> hydrophobicity to hydrophilicity
    color_check_box_list = []

    # folder paths for the asset generation:
    aa_font_type_path = f"{path}{sep}{font_type}"
    aa_index_box_path = f"{path}{sep}AA_letters_common"

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

    return list_recolor_aa, aa_compare, color_check_box_list
