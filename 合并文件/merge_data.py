import pandas as pd
from tqdm import tqdm
if __name__ == '__main__':
    import os
    source_dir = 'merge_input_file/'
    li = []
    for file_name in tqdm(os.listdir(source_dir),desc='read singe file...'):
        file_dir = os.path.join(source_dir,file_name)
        name = file_name.split('.')[0]
        df = pd.read_csv(file_dir,sep='\t',index_col=0,skiprows=1,header=None,names=[name])
        li.append(df)
        del(df)
    target_df = pd.concat(li,axis=1)
    target_df = target_df.sort_index()
    target_df.to_csv('merge_output.csv',sep='\t')

