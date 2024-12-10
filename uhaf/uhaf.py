import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import re


class uHAF:
    def __init__(self, uhaf_ex, target_sheetnames=None):
        """
        Initialize the uHAF class.

        Args:
            uhaf_ex (dict): Dictionary of dataframes extracted from the Excel sheets.
            target_sheetnames (list, optional): List of target sheet names to process. Default is None.
        """
        self.uhaf_ex = uhaf_ex
        self.sheet_names = list(self.uhaf_ex.keys())
        self.df_uhafs = {}
        self.dict_uhafs = {}
        self.df_uhaf_marker = {}
        self.generate_dict_uhafs(target_sheetnames)
        self.oran_list = self.sheet_names

    def track_cell_from_uHAF(self, sheet_name, cell_type_target):
        """
        Track the lineage of a cell type back to its root.

        Args:
            sheet_name (str): Name of the sheet.
            cell_type_target (str): Target cell type to trace.

        Returns:
            list: List of cell types from the root to the target.
        """
        trace = []
        while cell_type_target != 'Cell':
            if cell_type_target not in self.df_uhafs[sheet_name].index:
                print(f"Warning: {cell_type_target} not found in {sheet_name}")
                return trace
            tar_father = self.df_uhafs[sheet_name].loc[cell_type_target].father
            trace.append(cell_type_target)
            cell_type_target = tar_father
        trace.append('Cell')
        return trace[::-1]

    def cut_level_annotation(self, sheet_name, query_cell_type, level):
        """
        Get the annotation of a cell type at a specific level. level 3 is the celltype itself.

        Args:
            sheet_name (str): Name of the sheet.
            query_cell_type (str): Query cell type.
            level (int): Annotation level.

        Returns:
            str: Cell type at the specified level.
        """
        if level == 3:
            return query_cell_type
        trace = self.track_cell_from_uHAF(sheet_name, query_cell_type)
        if len(trace) <= level:
            return trace[-1]
        return trace[level]

    def set_annotation_level(self, query_cell_types, sheet_name, annotation_level):
        """
        Map cell types to a specified annotation level.

        Args:
            query_cell_types (list): List of query cell types.
            sheet_name (str): Name of the sheet.
            annotation_level (int): Annotation level.

        Returns:
            dict: Mapping of query cell types to their annotation levels.
        """
        annotation_level_map = {}
        for query_cell_type in query_cell_types:
            cell_type_retrieved = self.cut_level_annotation(sheet_name, query_cell_type, annotation_level)
            annotation_level_map[query_cell_type] = cell_type_retrieved
        return annotation_level_map

    def generate_dict_uhafs(self, sheetnames=None):
        """
        Generate dictionaries for uHAF structures.

        Args:
            sheetnames (list, optional): List of sheet names to process. Default is None.
        """
        if sheetnames:
            print('Generating uHAF for specified organs')
            for sheet_name in sheetnames:
                self.decode_child_father_pairs(sheet_name)
        else:
            print('No sheetnames specified, generating uHAF of every organ')
            for sheet_name in self.sheet_names:
                if contains_chinese(sheet_name):
                    continue
                if 'Sheet' in sheet_name:
                    continue
                self.decode_child_father_pairs(sheet_name)

    def create_child_father_pairs(self, start_level_cell_types):
        """
        Create child-father pairs for top-level cell types.

        Args:
            start_level_cell_types (list): List of top-level cell types.

        Returns:
            DataFrame: DataFrame containing child-father pairs.
        """
        return pd.DataFrame([(j, 'Cell') for j in start_level_cell_types], columns=['child', 'father'])

    def decode_child_father_pairs(self, sheet_name):
        """
        Decode child-father pairs from an Excel sheet.

        Args:
            sheet_name (str): Name of the sheet.
        """
        ex = self.uhaf_ex[sheet_name]
        celltypes_maker_excel = self.extract_celltype_and_markers(ex)
        df_0 = self.excel_to_paired_df_with_marker(celltypes_maker_excel)
        df_0.index = df_0['child']
        df_0['organ'] = sheet_name
        df_0 = df_0[['child', 'father', 'marker', 'organ']]
        self.df_uhafs[sheet_name] = df_0[['child', 'father', 'organ']]
        self.df_uhaf_marker[sheet_name] = df_0
        self.dict_uhafs[sheet_name] = self.convert_to_nested_dict(df_0[['child', 'father']])

    def find_marker_column(self, columns):
        """
        Find the marker column in the dataset.

        Args:
            columns (Index): Columns of the DataFrame.

        Returns:
            int: Index of the marker column.

        Raises:
            KeyError: If no marker column is found.
        """
        possible_columns = ['[Marker]', '|Marker|', '[Makers]', '|Makers|']
        for col in possible_columns:
            if col in columns:
                return columns.get_loc(col)
        raise KeyError("None of the marker columns found in the data.")

    def extract_celltype_and_markers(self, ex):
        """
        Extract cell types and markers from the sheet.

        Args:
            ex (DataFrame): DataFrame containing the sheet data.

        Returns:
            DataFrame: Extracted data with cell types and markers.
        """
        ex.columns = ex.columns.str.strip()
        try:
            marker_column_index = self.find_marker_column(ex.columns)
            celltype_data = ex.iloc[:, list(ex.columns).index('[Tissue]') + 1: marker_column_index]
            marker_data = ex.iloc[:, marker_column_index]
            celltype_data = celltype_data.dropna(axis=0, how='all').dropna(axis=1, how='all')
            celltype_data['marker'] = marker_data
        except KeyError as e:
            print(f"Error extracting cell type and marker columns: {e}")
            return pd.DataFrame()
        return celltype_data

    def excel_to_paired_df_with_marker(self, data):
        """
        Convert cell type data to paired child-father DataFrame with markers.

        Args:
            data (DataFrame): DataFrame with cell types and markers.

        Returns:
            DataFrame: Paired data with child-father relationships and markers.
        """
        if 'marker' not in data.columns:
            raise ValueError("The input DataFrame must contain a 'marker' column.")
        ans = []
        rows, cols = data.shape
        seen_children = set()
        for col in range(0, cols - 1):
            for row in range(rows):
                if not pd.isna(data.iloc[row, col]):
                    current_cell = data.iloc[row, col]
                    if current_cell in seen_children:
                        continue
                    father = None
                    for left_col in range(col - 1, -1, -1):
                        for up_row in range(row, -1, -1):
                            if not pd.isna(data.iloc[up_row, left_col]):
                                father = data.iloc[up_row, left_col]
                                break
                        if father:
                            break
                    marker = data.iloc[row, -1]
                    ans.append((current_cell, father or 'Cell', marker))
                    seen_children.add(current_cell)
        top_level_cells = data.iloc[:, 0].dropna().tolist()
        for row, cell in enumerate(top_level_cells):
            if cell not in seen_children:
                marker = data.iloc[row, -1]
                ans.append((cell, 'Cell', marker))
                seen_children.add(cell)
        df = pd.DataFrame(ans, columns=['child', 'father', 'marker'])
        df = df.drop_duplicates()
        df.index = df['child']
        return df

    def build_nested_dict(self, df, current_node):
        """
        Build a nested dictionary from the child-father pairs.

        Args:
            df (DataFrame): DataFrame with child-father relationships.
            current_node (str): Current node to process.

        Returns:
            dict: Nested dictionary representing the hierarchy.
        """
        children = df[df['father'] == current_node]['child'].tolist()
        if not children:
            return {}
        return {child: self.build_nested_dict(df, child) for child in children}

    def convert_to_nested_dict(self, df):
        """
        Convert the child-father DataFrame into a nested dictionary.

        Args:
            df (DataFrame): DataFrame with child-father relationships.

        Returns:
            dict: Nested dictionary with 'Cell' as the root.
        """
        return {'Cell': self.build_nested_dict(df, 'Cell')}

    def generate_uhaf_GPTs_prompts(self, sheet_name, custom_cell_types):
        """
        Generate uHAF mapping prompts for GPT tools.

        Args:
            sheet_name (str): Name of the sheet.
            custom_cell_types (list): List of custom cell types.

        Returns:
            str: Prompt for uHAF mapping.
        """
        prompts = f"The cell types are: \n{custom_cell_types}.\n"
        prompts += f'Copy the above cell types and paste them on the website (https://uhaf.unifiedcellatlas.org/#/uHAFMapping) to get the corresponding mapping dictionary.'
        return prompts

def contains_chinese(s):
    """
    Check if a string contains any Chinese characters.

    Args:
        s (str): Input string.

    Returns:
        bool: True if Chinese characters are found, False otherwise.
    """
    for char in s:
        if '\u4e00' <= char <= '\u9fff':
            return True
    return False


def get_latest_version(folder_path):
    version_pattern = r'uHAF(\d+\.\d+\.\d+)\.xlsx$'
    files = os.listdir(folder_path)
    versions = []

    for file in files:
        match = re.search(version_pattern, file)
        if match:
            versions.append(match.group(1))

    if versions:
        latest_version = max(versions, key=lambda v: [int(x) for x in v.split('.')])
        return latest_version
    else:
        return None

def build_uhaf(latest: bool = True, 
               uhaf_xlsx_version: str = None, 
               target_sheetnames: list = None, 
               uhaf_path: str = None) -> uHAF:
    """
    Build a uHAF instance from an Excel file.

    Args:
        latest (bool, optional): Flag to fetch the latest uHAF version. Default is True.
        uhaf_xlsx_version (str, optional): Version of the uHAF Excel file. Default is None.
        target_sheetnames (list, optional): List of target sheet names. Default is None.
        uhaf_path (str, optional): Path to the uHAF Excel file. Default is None.

    Returns:
        uHAF: Instance of the uHAF class.
    """
    if not uhaf_path:
        if latest:
            folder_path = os.path.join(os.path.dirname(__file__), 'reference')
            uhaf_xlsx_version = get_latest_version(folder_path)
            print('Using the latest uHAF version:', uhaf_xlsx_version)
        uhaf_path = os.path.join(os.path.dirname(__file__), 'reference', f'uHAF{uhaf_xlsx_version}.xlsx')
    uhaf_ex = pd.read_excel(uhaf_path, sheet_name=None)
    return uHAF(uhaf_ex, target_sheetnames)