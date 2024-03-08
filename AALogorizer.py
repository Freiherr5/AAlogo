import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
from PIL import Image
import pandas as pd
import numpy as np
import GetAA
import StandardConfig
from configparser import ConfigParser
import ast
import LogoUtil


class AALogoGenerator:

    # Intilialize with AALogoGenerator()
    def __init__(self, set_legend: bool, list_columns: list):
        self.set_legend = set_legend               # bool value for legend asset
        self.list_columns = list_columns           # shape: [AA_sequence, TMD_start, TMD_stop]

    def list_slicer(self, df, length_tmd, length_jmd):
        """
        create 2 DataFrames:
        1. DataFrame with JMD_N / TMD (list_n)
        2. DataFrame with TMD / JMD_C (list_c)

        Parameters
        __________
        df : pandas.DataFrame
        length_tmd : window size of tmd (number of amino acid residues shown)
        length_jmd : window size of jmd (number of amino acid residues shown)

        Returns
        _______
        list_n : (JMD_N sequence part, TMD sequence part)
        list_c : (TMD sequence part, JMD_C sequence part)
        """

        list_n = []
        list_c = []
        df = df.reset_index()  # make sure indexes pair with number of rows

        for index, row in df.iterrows():
            # N-term
            start_pos_tmd = int(row[self.list_columns[1]])-1
            seq_slice_n = [row[self.list_columns[0]][start_pos_tmd-int(length_jmd):start_pos_tmd],
                           row[self.list_columns[0]][start_pos_tmd:start_pos_tmd + int(length_tmd)]]
            seq_slice_n_check = [LogoUtil.check_input_string(seq) for seq in seq_slice_n]  # check sequence
            list_n.append(seq_slice_n_check)

            # C-term
            stop_pos_tmd = int(row[self.list_columns[2]])
            seq_slice_c = [row[self.list_columns[0]][stop_pos_tmd-int(length_tmd):stop_pos_tmd],
                           row[self.list_columns[0]][stop_pos_tmd:stop_pos_tmd + int(length_jmd)]]
            seq_slice_c_check = [LogoUtil.check_input_string(seq) for seq in seq_slice_c]
            list_c.append(seq_slice_c_check)

        return list_n, list_c

    @staticmethod
    def data_frame_aa_propensities(list_n_or_c, length_tmd, length_jmd, aa_list):
        """
        Position resolved propensities of amino acids of all sequences

        Parameters
        __________
        list_n_or_c : return list from list_slicer, has amino acid sequence of both domains as list
        length_tmd : window size of tmd (number of amino acid residues shown)
        length_jmd : window size of jmd (number of amino acid residues shown)
        aa_list : order of the amino acids (assigned by the LogoStyle.ini file), from top to bottom

        Returns
        _______
        aa_propensity_df : creates dataframe with the propensity of the amino acids (from 0 to 1 normalized)
        """

        list_concat = []
        for items in list_n_or_c:
            list_concat.append(f"{items[0]}{items[1]}")
        df_intermediate = pd.Series(list_concat).str.split(r"", expand=True).drop([0, length_jmd + length_tmd + 1],
                                                                                  axis=1)

        aa_list_transform = [aa_list]
        for slices in df_intermediate:
            aa_20 = []
            list_aa_slice = df_intermediate[slices].to_numpy().tolist()
            for AA in aa_list:
                aa_20.append(list_aa_slice.count(AA)/len(df_intermediate))
            aa_list_transform.append(aa_20)
        aa_propensity_df = pd.DataFrame(aa_list_transform, columns=aa_list).transpose().drop(0, axis=1)
        return aa_propensity_df

    def make_logo(self, df, name: str, length_tmd: int, length_jmd: int, aa_config_section_name: str = "OG_AA_config",
                  font_type: str = "bold_AA_fonts", config_set: bool = True, color_grad: list = None,
                  order_aa: list = None):    # name of the plot, length of the domains for JMD and TMD
        """
        Generates the final plot/logo

        Parameters
        __________
        df : pd.DataFrame
        name : title of generated sequence window
        length_tmd : window size of tmd (number of amino acid residues shown)
        length_jmd : window size of jmd (number of amino acid residues shown)
        aa_config_section_name : LogoStyle.ini configuration file section name with the amino acid categories
        font_type : classic_AA_fonts, bold_AA_fonts, modern_AA_fonts (writing style package)
        config
        """

        path_file, sep = StandardConfig.find_folderpath()
        assets_path = f"{path_file}{sep}AA_letters_common{sep}"

        # check external parameters for config / gradient
        # _____________________________________________________________________________________________________________
        list_n, list_c = AALogoGenerator.list_slicer(self, df, length_tmd, length_jmd)
        # check for get_aa_list if input is correct!
        if not isinstance(config_set, bool):
            raise ValueError(f"config parameter needs to be bool")
        elif config_set is False:
            if color_grad is not None:
                if not isinstance(color_grad, list):
                    raise ValueError(f"color_grad needs to have the following shape: [[r1,g1,b1],[r1,g1,b1]],"
                                     f"where r,g,b is int or float")
                elif np.array(color_grad).shape != (2, 3):
                    raise ValueError(f"color_grad needs to have the following shape: [[r1,g1,b1],[r1,g1,b1]],"
                                     f"where r,g,b is int or float")
                test_color_grad = color_grad[0]
                test_color_grad.extend(color_grad[1])
                for values in test_color_grad:
                    if not isinstance(values, (int, float)):
                        raise ValueError(f"given list value {values} is not an allowed datatype,"
                                         f"please enter an int or float value")

            if order_aa is not None:
                if not isinstance(order_aa, list):
                    raise ValueError(f"order parameter needs to be list type")
        # _____________________________________________________________________________________________________________

        get_aa_list = GetAA.aa_image_colorizer(aa_config_section_name, font_type, config_set, color_grad, order_aa)
        df_n = AALogoGenerator.data_frame_aa_propensities(list_n, length_tmd, length_jmd, get_aa_list[1])
        df_c = AALogoGenerator.data_frame_aa_propensities(list_c, length_tmd, length_jmd, get_aa_list[1])

        df_plot_list = [df_n, df_c]

        for dfs in df_plot_list:

            if dfs is df_c:
                length_tmd, length_jmd = length_jmd, length_tmd

            matplotlib.rcParams["axes.linewidth"] = 3
            font = {"weight": "bold",
                    "size": 18}
            matplotlib.rc("font", **font)

            # matplotlib implementation
            fig, ax = plt.subplots(dpi=254)
            fig.set_size_inches((length_tmd + length_jmd)*2, 10, forward=True)  # 1 inch = 100 pixels

            ax.set_xlim(-length_jmd-0.5, length_tmd-0.5)
            ax.set_ylim(0, 1)
            ax.set_xticks(np.arange(-length_jmd, length_tmd, 1))

            # define the figure
            ax.set_xlabel("sequence position", fontsize=18, weight="bold")
            ax.set_ylabel("AA frequency", fontsize=18, weight="bold")
            ax.spines[['right', 'top']].set_visible(False)

            # draw areas of the domains

            # color code tmd / jmd
            color_tmd = "#d9bd82"
            color_jmd = "#99c0de"
            color_gradient = None

            config_file = "LogoStyle.ini"
            config = ConfigParser()
            config.read(config_file)
            config.sections()

            # get color codes if mentioned
            bg_style = list(config["bg_style"])

            if bg_style[0] == "tmd":
                r, g, b = ast.literal_eval(config["bg_style"][bg_style[0]])
                color_tmd = LogoUtil.rgb_to_hex(r, g, b)                            # hex
            if bg_style[1] == "jmd":
                r, g, b = ast.literal_eval(config["bg_style"][bg_style[1]])
                color_jmd = LogoUtil.rgb_to_hex(r, g, b)                            # hex
            if bg_style[2] == "gradient":
                color_gradient = ast.literal_eval(config["bg_style"][bg_style[2]])  # list with rgb values

            # N-term or C-term
            if dfs is df_n:
                ax.add_patch(Rectangle((-length_jmd-0.5, 0), length_jmd, 1, color=color_jmd))  # blue
                ax.add_patch(Rectangle((-0.5, 0), length_tmd, 1, color=color_tmd))  # orange
                name_fig = "N_term"

                # right gradient
                if color_gradient is not list:
                    img = Image.open(f"{assets_path}white_r_grad.png")
                else:
                    img = LogoUtil.convert_image_color(assets_path, "white_r_grad", tuple(color_gradient))
                imagebox = OffsetImage(img, zoom=0.554)
                # box_alignment position is upper middle of picture
                r_grad = AnnotationBbox(imagebox, xy=(-0.5, 0), box_alignment=(1, 0), frameon=False)
                imagebox.image.axes = ax
                ax.add_artist(r_grad)

                # text region sequences
                ax.text((-length_jmd-1.6) / 2, 1.02, "JMD-N", fontsize=22, weight="bold")
                ax.text((length_tmd-1.5) / 2, 1.02, "TMD", fontsize=22, weight="bold")

            else:
                ax.add_patch(Rectangle((-length_jmd-0.5, 0), length_jmd, 1, color=color_tmd))  # orange
                ax.add_patch(Rectangle((-0.5, 0), length_tmd, 1, color=color_jmd))  # blue
                name_fig = "C_term"

                # left gradient
                if color_gradient is not list:
                    img = Image.open(f"{assets_path}white_l_grad.png")
                else:
                    img = LogoUtil.convert_image_color(assets_path, "white_l_grad", tuple(color_gradient))
                imagebox = OffsetImage(img, zoom=0.554)
                # box_alignment position is upper middle of picture
                l_grad = AnnotationBbox(imagebox, xy=(-0.5, 0), box_alignment=(0, 0), frameon=False)
                imagebox.image.axes = ax
                ax.add_artist(l_grad)

                # text regions sequence
                ax.text((-length_jmd-1.5) / 2, 1.02, "TMD", fontsize=22, weight="bold")
                ax.text((length_tmd-1.6) / 2, 1.02, "JMD-C", fontsize=22, weight="bold")

            """
            code part to insert letters of clustered sequences
            must be a nested loop (in a loop)
            get dataframe saved in list for each position in the sequence
            dataframe contains amino acid frequency
            """

            # if gradient is used, switch off index_color_boxes automatically
            if get_aa_list[2] is None:
                self.set_legend = False

            # add index box
            if self.set_legend:
                list_index_color_boxes = get_aa_list[2]
                i = 1
                for entries in list_index_color_boxes:
                    index_image = entries[1]
                    index_text_string = entries[0]
                    x, y = (length_tmd-0.50), i - 0.03
                    # important for converting it into an usable format for AnnotationBbox
                    imagebox = OffsetImage(index_image, zoom=0.3)
                    # AnnotationBbox for translation
                    # box_alignment position is upper middle of picture
                    index_box = AnnotationBbox(imagebox, xy=(x, y), box_alignment=(-0.6, 0.35), frameon=False)
                    imagebox.image.axes = ax
                    ax.add_artist(index_box)
                    ax.text(x+0.51, y, index_text_string, fontsize=14, weight="bold")
                    i = i-0.05

            # add the AA letters
            for columns in dfs:
                i = 0
                concat_distance = 0
                aa_pos_column_list = dfs[columns].tolist()
                while i < len(dfs):
                    concat_distance += aa_pos_column_list[i]
                    # Python program to change the ratio of height and width of an image
                    # Taking image as input
                    if aa_pos_column_list[i] > 0:
                        img = get_aa_list[0][i][1]
                        # Changing the height and width of the image
                        factor = aa_pos_column_list[i]  # get info from dataframe!
                        width = 110
                        height = int(554*factor)+1  # conserved height
                        # Resizing the image
                        img = img.resize((width, height))
                        # important for converting it into an usable format for AnnotationBbox
                        imagebox = OffsetImage(img, zoom=1)
                        # AnnotationBbox for translation
                        # box_alignment position is upper middle of picture
                        ab = AnnotationBbox(imagebox, xy=(columns - (length_jmd+1), 1-concat_distance),
                                            box_alignment=(0.5, 0), frameon=False)
                        imagebox.image.axes = ax
                        ax.add_artist(ab)
                    i += 1
            plt.savefig(f"{path_file}{sep}Output{sep}{name}_{name_fig}.png", bbox_inches='tight')
