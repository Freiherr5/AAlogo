# standard libs
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
from PIL import Image
import pandas as pd
import numpy as np
from configparser import ConfigParser
import ast
# intern
import GetAA
import StandardConfig
from StandardConfig import timingmethod
import LogoUtil


class _AALogoGenerator:

    # Intilialize with AALogoGenerator()
    def __init__(self, set_legend: bool, list_columns: list, start_pos: bool = True):
        self.set_legend = set_legend               # bool value for legend asset
        self.list_columns = list_columns           # shape: [aa_sequence: str, pos_seq: str]
        self.start_pos = start_pos

    # Internal Processes for AAlogo generation
    # __________________________________________________________________________________________________________________
    def _list_slicer(self, df, length_right, length_left):
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
        list_seq_pos = []

        df = df.reset_index()  # make sure indexes pair with number of rows

        for index, row in df.iterrows():
            # N-term
            if self.start_pos:
                start_pos_tmd = int(row[self.list_columns[1]])-1
                if start_pos_tmd-int(length_left) >= 0:
                    seq_slice_n = [row[self.list_columns[0]][start_pos_tmd - int(length_left):start_pos_tmd],
                                   row[self.list_columns[0]][start_pos_tmd:start_pos_tmd + int(length_right)]]
                else:
                    seq_slice_n_start_pre = row[self.list_columns[0]][0:start_pos_tmd]
                    overhang_start = "-"*(int(length_left)-start_pos_tmd)  # "-" take in account shorter seq
                    seq_slice_n_start = f"{overhang_start}{seq_slice_n_start_pre}"
                    seq_slice_n = [seq_slice_n_start,
                                   row[self.list_columns[0]][start_pos_tmd:start_pos_tmd + int(length_right)]]
                seq_slice_n_check = [LogoUtil.check_input_string(seq) for seq in seq_slice_n]  # check sequence
                list_seq_pos.append(seq_slice_n_check)

            # C-term
            else:
                stop_pos_tmd = int(row[self.list_columns[1]])
                seq_slice_c = [row[self.list_columns[0]][stop_pos_tmd-int(length_left):stop_pos_tmd],
                               row[self.list_columns[0]][stop_pos_tmd:stop_pos_tmd + int(length_right)]]
                seq_slice_c_check = [LogoUtil.check_input_string(seq) for seq in seq_slice_c]
                list_seq_pos.append(seq_slice_c_check)

        return list_seq_pos

    @staticmethod
    def _data_frame_aa_propensities(list_seq_pos, length_right, length_left, aa_list):
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
        for items in list_seq_pos:
            list_concat.append(f"{items[0]}{items[1]}")
        df_intermediate = (pd.Series(list_concat).str.split(r"", expand=True)
                           .drop([0, length_left + length_right + 1], axis=1))
        df_intermediate.replace("-", None, inplace=True)

        aa_list_transform = [aa_list]
        for slices in df_intermediate:
            aa_20 = []
            list_aa_slice = df_intermediate[slices].to_numpy().tolist()
            len_slice = len(list_aa_slice)
            for AA in aa_list:
                aa_20.append(list_aa_slice.count(AA)/len_slice)  # before df_intermediate
            aa_list_transform.append(aa_20)
        aa_propensity_df = pd.DataFrame(aa_list_transform, columns=aa_list).transpose().drop(0, axis=1)
        return aa_propensity_df
    # __________________________________________________________________________________________________________________

    def _make_logo(self, df, name: str, length_right: int, length_left: int,
                   aa_config_section_name: str = "OG_AA_config", font_type: str = "bold_AA_fonts",
                   config_set: bool = True, color_grad: list = None, order_aa_grad: list = None,
                   color_advance: list = None, list_title_sides: list = None):
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
        config_set : turn on and off the configuration color schema --> enables color_grad
        color_grad : list of two colors, color gradient from top (first list color) to bottom (second list color)
                     colors are written as RGB lists --> [R, G, B]
        order_aa_grad : ordering amino acids for gradient function
        color_advance : int or float values which are normalized and applied to the color gradient
        list_title_sides : list of the titles for both sides of the plot (separated by the start/stop position)
        """

        path_file, sep = StandardConfig.find_folderpath()
        assets_path = f"{path_file.split("aalogo")[0]}fonts{sep}AA_letters_common{sep}"

        list_seq_pos = _AALogoGenerator._list_slicer(self, df, length_right, length_left)

        # get the necessary dataframes and image lists for AAlogo generation
        # ______________________________________________________________________________________________________________
        get_aa_list = GetAA.aa_image_colorizer(aa_config_section_name, font_type, config_set, color_grad,
                                               order_aa_grad, color_advance)
        df_propensity = _AALogoGenerator._data_frame_aa_propensities(list_seq_pos, length_right, length_left,
                                                                     get_aa_list[1])

        # matplotlib.pyplot based visualization
        # ______________________________________________________________________________________________________________
        matplotlib.rcParams["axes.linewidth"] = 3
        font = {"weight": "bold",
                "size": 18}
        matplotlib.rc("font", **font)

        # matplotlib implementation
        fig, ax = plt.subplots(dpi=254)
        fig.set_size_inches((length_right + length_left)*2, 10, forward=True)  # 1 inch = 100 pixels

        ax.set_xlim(-length_left-0.5, length_right-0.5)
        ax.set_ylim(0, 1)
        ax.set_xticks(np.arange(-length_left, length_right, 1))

        # define the figure
        ax.set_xlabel("sequence position", fontsize=18, weight="bold")
        ax.set_ylabel("AA frequency", fontsize=18, weight="bold")
        ax.spines[['right', 'top']].set_visible(False)

        # draw areas of the domains

        # color code tmd / jmd --> original membrane focused utility of visualization
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

        # generates the background with the gradient png on top
        # is it a start position or a stop position (important differentiation)
        if self.start_pos:
            ax.add_patch(Rectangle((-length_left-0.5, 0), length_left, 1, color=color_jmd))  # blue
            ax.add_patch(Rectangle((-0.5, 0), length_right, 1, color=color_tmd))  # orange

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

        else:
            ax.add_patch(Rectangle((-length_left-0.5, 0), length_left, 1, color=color_tmd))  # orange
            ax.add_patch(Rectangle((-0.5, 0), length_right, 1, color=color_jmd))  # blue

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

        """
        code part to insert letters of clustered sequences
        must be a nested loop (in a loop)
        get dataframe saved in list for each position in the sequence
        dataframe contains amino acid frequency
        """

        if list_title_sides is not None:
            if isinstance(list_title_sides, list) and (len(list_title_sides) == 2):
                # text region sequences
                ax.text((-length_left - 1.55) / 2, 1.02, str(list_title_sides[0]), fontsize=22, weight="bold")
                ax.text((length_right - 1.55) / 2, 1.02, str(list_title_sides[1]), fontsize=22, weight="bold")

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
                x, y = (length_right-0.50), i - 0.03
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
        for columns in df_propensity:
            i = 0
            concat_distance = 0
            aa_pos_column_list = df_propensity[columns].tolist()
            while i < len(df_propensity):
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
                    ab = AnnotationBbox(imagebox, xy=(columns - (length_left+1), 1-concat_distance),
                                        box_alignment=(0.5, 0), frameon=False)
                    imagebox.image.axes = ax
                    ax.add_artist(ab)
                i += 1

        # last naming differentiation depending on start/stop position (make better)
        if self.start_pos:
            start_tag = "set_start_true"
        else:
            start_tag = "set_start_false"
        plt.savefig(f"{path_file}{sep}Output{sep}{name}_{start_tag}.png", bbox_inches='tight')


class AAlogoMaker:

    def __init__(self, df: pd.DataFrame, name: str, column_seq: str, *args_position: str):
        self.df = df
        self.name = name
        self.column_seq = column_seq
        self.args_position = args_position

    def _check_self(self):
        # check df
        if not isinstance(self.df, pd.DataFrame):
            raise TypeError("needs to be pd.DataFrame")
        # check if column names are in df
        col_df_list = self.df.columns.tolist()
        if self.column_seq not in col_df_list:
            raise ValueError(f"{self.column_seq} not in pd.DataFrame.columns")
        for arg_pos in self.args_position:
            if arg_pos not in col_df_list:
                raise ValueError(f"{arg_pos} not in pd.DataFrame.columns")

    @staticmethod
    def _check_function_inputs(dict_input):
        dict_exchange = {"start_pos": True, "aa_right": 5, "aa_left": 5, "font_type": "bold_AA_fonts",
                         "custom_color": None, "config_name": None, "headers": None}

        # simple checking of variables
        # ______________________________________________________________________________________________________________
        dict_typing = {"start_pos": bool, "aa_right": int, "aa_left": int, "font_type": str,
                       "custom_color": list, "config_name": str, "headers": list}
        dict_input_keys = dict_input.keys()
        for keys in dict_input_keys:
            if not isinstance(dict_input[keys], dict_typing[keys]):
                dict_input[keys] = dict_exchange[keys]

        # more complex variables check
        # ______________________________________________________________________________________________________________
        for aa_lengths in ["aa_right", "aa_left"]:
            if dict_input[aa_lengths] > 40:
                print(f"{dict_input[aa_lengths]} exceeds the length limit for each side (max = 40 amino acids)")
                dict_input[aa_lengths] = 40
            if dict_input[aa_lengths] < 1:
                print(f"{dict_input[aa_lengths]} is below the minimal length for each side (min = 1 amino acids)")
                dict_input[aa_lengths] = 1

        if "custom_colors" in dict_input_keys:
            if dict_input["headers"] is not None:
                if np.array(dict_input["custom_colors"]).shape != (2, 3):
                    raise ValueError(f"custom_colors needs to have the following shape: [[r1,g1,b1],[r1,g1,b1]],"
                                     f"where r,g,b is int or float")
                else:
                    for sublists in dict_input["custom_colors"]:
                        for items in sublists:
                            if not isinstance(items, (int, float)):
                                ValueError(f"custom_colors needs to have the following shape: [[r1,g1,b1],[r1,g1,b1]],"
                                           f"where r,g,b is int or float")

        if "headers" in dict_input_keys:
            if dict_input["headers"] is not None:
                if np.array(dict_input["headers"]).shape != (2,):
                    raise ValueError(f"headers must bei [left_side_title, right_side_title] in shape!")

        if "config_name" in dict_input_keys:
            config_file = "LogoStyle.ini"
            config = ConfigParser()
            config.read(config_file)
            config_sections = config.sections()
            if dict_input["config_name"] not in config_sections:
                dict_input["config_name"] = None
        return dict_input

    @staticmethod
    def help():
        print("""
        Descriptionstart_pos_tmd - int(length_left)
        ___________
        AAlogoMaker is meant to ease the usage of AALogoGenerator
                        ________________________________________________________________________________________________
                        df: pd.DataFrame, name: str, column_seq: str, *args_position: str
                        df --> needs to contain the amino acid sequences = column_seq
                           --> needs at least one *args_position = name of column with alignment positions for 
                               sequences in column_seq
                        ________________________________________________________________________________________________
                        
        Callable Functions
        __________________
        single_mode() : application for sequence-propensity visualization based on a single alignment position
                        ________________________________________________________________________________________________
                        start_pos: bool = True, aa_right: int = 5, aa_left: int = 5, font_type: str = "bold_AA_fonts",
                        theme: str = "Kyte-Doolittle", custom_colors: list = None, config_name: str = None,
                        headers: list = None
                        ________________________________________________________________________________________________
        tmd_mode() : application for sequence-propensity visualization based on start and stop position of a tmd 
                     (usage for transmembrane proteins)
                     ___________________________________________________________________________________________________
                     start_pos: bool = True, aa_jmd: int = 5, aa_tmd: int = 5, font_type: str = "bold_AA_fonts",
                     theme: str = "Kyte-Doolittle", custom_colors: list = None, config_name: str = None,
                     headers: list = None
                     ___________________________________________________________________________________________________
        
        Note
        ____
        Inputting a config setting is dominant over a gradient setting!
        
        Config
        ______
        Structure: [amino acids, RGB-color-code]
        Example:
        [OG_AA_config]
        Aromatic =  ["F,Y,W", [216,116,45]]
        Aspartic = ["V,A,L,I,M", [223,168,70]]
        Neutral = ["S,T,N,Q", [255,255,255]]
        Acidic = ["D,E", [224,60,60]]
        Basic = ["H,K,R", [60,60,224]]
        Special = ["P,G", [50,50,50]]
        Cysteine = ["C", [224,90,224]]
        
        Gradient
        ________
        AAlogo_Util.single_mode().available_themes
        --> see all available scale-themes
        """)

    @timingmethod
    def single_mode(self, start_pos: bool = True, aa_right: int = 5, aa_left: int = 5, font_type: str = "bold_AA_fonts",
                    theme: str = "Kyte-Doolittle", custom_colors: list = None, config_name: str = None,
                    headers: list = None):

        # check inputs
        # ______________________________________________________________________________________________________________
        AAlogoMaker._check_self(self)
        dict_inputs = {"start_pos": start_pos, "aa_right": aa_right, "aa_left": aa_left, "font_type": font_type,
                       "custom_color": custom_colors, "config_name": config_name, "headers": headers}
        dict_inputs = AAlogoMaker._check_function_inputs(dict_input=dict_inputs)


        if dict_inputs["config_name"] is not None:
            config_set = True
        else:
            config_set = False

        # get AA order / theme for AAlogo
        # ______________________________________________________________________________________________________________
        data_hydrophobicity_scales = pd.read_excel("scales_hydrophobicity.xlsx")
        if theme not in data_hydrophobicity_scales.columns.tolist():
            theme = "Kyte-Doolittle"
        theme_df = (data_hydrophobicity_scales[["aa_code", theme]].sort_values(by=theme, ascending=False)
                    .reset_index().drop("index", axis=1))

        order_aa_grad = theme_df["aa_code"]
        color_advance = theme_df[theme]
        # config is always dominant if name is given!
        set_legend = False
        # config is always dominant if name is given!
        if config_set:
            order_aa_grad, color_advance = None, None
            set_legend = True
        # have an attribute with all sorting options
        # ______________________________________________________________________________________________________________
        AAlogoMaker.single_mode.available_themes = [item for item in
                                                    data_hydrophobicity_scales.columns.tolist()[2:]]

        for arg_pos in self.args_position:
            init_aalogo = _AALogoGenerator(set_legend=set_legend, list_columns=[self.column_seq, arg_pos],
                                           start_pos=start_pos)

            # Plot generation command
            # make_logo(df, name, length_tmd, length_jmd, aa_config_section_name, font_type="classic_AA_fonts")
            init_aalogo.make_logo(df=self.df, name=str(self.name), length_right=dict_inputs["aa_right"],
                                  length_left=dict_inputs["aa_left"], font_type=dict_inputs["font_type"],
                                  config_set=config_set, aa_config_section_name=dict_inputs["config_name"],
                                  order_aa_grad=order_aa_grad, color_advance=color_advance,
                                  list_title_sides=dict_inputs["headers"], color_grad=custom_colors)

    @timingmethod
    def tmd_mode(self, start_pos: bool = True, aa_jmd: int = 5, aa_tmd: int = 5, font_type: str = "bold_AA_fonts",
                 theme: str = "Kyte-Doolittle", custom_colors: list = None, config_name: str = None,
                 headers: list = None):
        # check inputs
        # ______________________________________________________________________________________________________________
        AAlogoMaker._check_self(self)
        dict_inputs = {"start_pos": start_pos, "aa_right": aa_tmd, "aa_left": aa_jmd, "font_type": font_type,
                       "custom_color": custom_colors, "config_name": config_name, "headers": headers}
        dict_inputs = AAlogoMaker._check_function_inputs(dict_input=dict_inputs)

        if dict_inputs["config_name"] is not None:
            config_set = True
        else:
            config_set = False

        # get AA order / theme for AAlogo
        # ______________________________________________________________________________________________________________
        data_hydrophobicity_scales = pd.read_excel("scales_hydrophobicity.xlsx")
        if theme not in data_hydrophobicity_scales.columns.tolist():
            theme = "Kyte-Doolittle"
        theme_df = (data_hydrophobicity_scales[["aa_code", theme]].sort_values(by=theme, ascending=False)
                    .reset_index().drop("index", axis=1))

        order_aa_grad = theme_df["aa_code"]
        color_advance = theme_df[theme]
        set_legend = False
        # config is always dominant if name is given!
        if config_set:
            order_aa_grad, color_advance = None, None
            set_legend = True
        # have an attribute with all sorting options
        # ______________________________________________________________________________________________________________
        AAlogoMaker.tmd_mode.available_themes = [item for item in
                                                 data_hydrophobicity_scales.columns.tolist()[2:]]

        # algorithm with TMD mode, check only two *args allowed for start and stop position
        # ______________________________________________________________________________________________________________
        if len(self.args_position) > 2:
            self.args_position = self.args_position[:2]
            print("Warning! Given columns exceeds limit of allowed TMD positions; max = 2")
        elif len(self.args_position) < 2:
            raise Warning("Only 1 position given, TMD has a start and a stop position!")

        for arg_pos in self.args_position:
            init_aalogo = _AALogoGenerator(set_legend=set_legend, list_columns=[self.column_seq, arg_pos],
                                           start_pos=start_pos)

            # Plot generation command
            # make_logo(df, name, length_tmd, length_jmd, aa_config_section_name, font_type="classic_AA_fonts")
            init_aalogo.make_logo(df=self.df, name=str(self.name), length_right=dict_inputs["aa_right"],
                                  length_left=dict_inputs["aa_left"], font_type=dict_inputs["font_type"],
                                  config_set=config_set, aa_config_section_name=dict_inputs["config_name"],
                                  order_aa_grad=order_aa_grad, color_advance=color_advance,
                                  list_title_sides=dict_inputs["headers"], color_grad=custom_colors)
            start_pos = False  # change orientation for stop position
            if dict_inputs["headers"] is not None:
                dict_inputs["headers"].reverse()
            dict_inputs["aa_right"], dict_inputs["aa_left"] = dict_inputs["aa_left"], dict_inputs["aa_right"]
