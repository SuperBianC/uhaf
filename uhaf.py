import os
import pandas as pd
import numpy as np
from tqdm import tqdm


class uHAF:
    def __init__(self, uhaf_ex,target_sheetnames = None):
        self.uhaf_ex = uhaf_ex
        self.sheet_names = list(self.uhaf_ex.keys())
        self.df_uhafs = {}
        self.dict_uhafs = {}
        self.generate_dict_uhafs(target_sheetnames)

    def track_cell_from_uHAF(self, sheet_name, cell_type_target):
        trace = []
        while cell_type_target != 'Cell':
            tar_father = self.df_uhafs[sheet_name].loc[cell_type_target].father
            trace.append(cell_type_target)
            cell_type_target = tar_father
        trace.append('Cell')
        return trace[::-1]
    
    def cut_level_annotation(self, sheet_name, query_cell_type, level):
        if level == 3:
            return query_cell_type
        trace = self.track_cell_from_uHAF(sheet_name, query_cell_type)
        if len(trace) <= level:
            return trace[-1]
        return trace[level]

    def set_annotation_level(self, query_cell_types, sheet_name, annotation_level):
        annotation_level_map = {}
        for query_cell_type in query_cell_types:
            cell_type_retrieved = self.cut_level_annotation(sheet_name, query_cell_type, annotation_level)
            annotation_level_map[query_cell_type] = cell_type_retrieved
        return annotation_level_map

    def generate_dict_uhafs(self, sheetnames = None):
        if sheetnames:
            print('Generating uHAF for specified organs')
            for sheet_name in tqdm(sheetnames):
                self.decode_child_father_pairs(sheet_name)
        else:
            print('No sheetnames specified, generating uHAF for all organs')
            for sheet_name in tqdm(self.sheet_names, desc = 'Generating uHAF'):
                if contains_chinese(sheet_name):continue
                if 'Sheet' in sheet_name:continue
                self.decode_child_father_pairs(sheet_name)
            
    def decode_child_father_pairs(self,sheet_name):
        # print(sheet_name)
        ex = self.uhaf_ex[sheet_name]
        df_0=pd.DataFrame()
        try:
            top_level_cell_types=ex.iloc[:,list(ex.columns).index('[Tissue]'):list(ex.columns).index('[Marker]')]
        except:
            top_level_cell_types=ex.iloc[:,list(ex.columns).index('[Tissue]'):list(ex.columns).index('|Marker|')]
        l0=[]
        if '[Subtissue]' in top_level_cell_types.columns:
            start_level_cell_types = top_level_cell_types['[Subtissue]']
        elif '|Subtissue|' in top_level_cell_types.columns:
            start_level_cell_types = top_level_cell_types['|Subtissue|']
        else:print('no subtissue')

        for j in start_level_cell_types:
            l0.append([j,'Cell'])
        l0=pd.DataFrame(l0)
        l0.columns=['child','father']
        df_0=pd.concat([df_0,l0])
        try:
            celltypes=ex.iloc[:,list(ex.columns).index('[Tissue]')+1:list(ex.columns).index('|Marker|')] ## +1 skip tissues
        except:
            celltypes=ex.iloc[:,list(ex.columns).index('[Tissue]')+1:list(ex.columns).index('[Marker]')]
        celltypes=celltypes.dropna(axis=0,how='all').dropna(axis=1,how='all')
        cf_pairs=self.pairs_to_df(celltypes)
        df_0=pd.concat([df_0,cf_pairs])
        df_0=df_0.drop_duplicates()
        df_0.index=range(len(df_0))
        df_0.organ = sheet_name
        df_0 = df_0.dropna()
        df_0.index = df_0.child
        self.df_uhafs[sheet_name] = df_0
        self.dict_uhafs[sheet_name] = self.convert_to_nested_dict(df_0[['child','father']])
        
    def pairs_to_df(self, data):
        lis2=[]
        for i in range(len(data.nunique())):
            if data.nunique().iloc[i]==0:
                lis2.append(data.columns[i])
        for i in lis2:
            data.drop(i,axis = 1,inplace = True)
        a,b=data.shape
        lis=[]
        ans=[]
        markers=range(len(data))
        for i in range(a-1,-1,-1):
            for j in range(b-1,-1,-1):
                if str(data.values[i][j]) !='nan':
                    lis.append([data.values[i][j],i,j])
        def zi_fu(lis,markers):
            cis=lis[0][2]
            for i in lis:
                if cis>i[2]:
                    ans.append(([lis[0][0],i[0],markers[0]]))
                    markers=markers[1:]
                    break
            return lis[1:],markers
        while len(lis)>=1:
            lis,markers=zi_fu(lis,markers)
        df=pd.DataFrame(ans)
        df=df[::-1]
        df.columns=['child','father','Description']
        df.index=df['child']
        df=df.drop('Description',axis=1)
        return df
    
    def build_nested_dict(self, df, current_node):
        children = df[df['father'] == current_node]['child'].tolist()
        if not children:
            return {}
        return {child: self.build_nested_dict(df, child) for child in children}

    def convert_to_nested_dict(self, df):
        return { 'Cell': self.build_nested_dict(df, 'Cell') }
    

    def generate_uhaf_GPTs_prompts(self, sheet_name, custom_cell_types):
        prompts = f"The cell types are: \n{custom_cell_types}."
        return prompts

def contains_chinese(s):
    for char in s:
        if '\u4e00' <= char <= '\u9fff':
            return True
    return False

def build_uhaf(uhaf_xlsx_version, target_sheetnames = None):
    uhaf_path = os.path.join(os.path.dirname(__file__), 'reference', uhaf_xlsx_version + '.xlsx')
    uhaf_ex = pd.read_excel(uhaf_path,sheet_name=None)
    return uHAF(uhaf_ex,target_sheetnames)