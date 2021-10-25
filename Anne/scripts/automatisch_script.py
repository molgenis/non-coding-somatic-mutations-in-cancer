
import pandas as pd

path_file = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/EGAD00001000292_metadata/delimited_maps/Sample_File.map"
head_path = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/chr22/"

df = pd.read_csv(path_file, sep='\t', header=None)
df_selection = pd.DataFrame()
df_selection[['Sample','type']] = df[0].str.split("_",expand=True,)
df_selection['Name'] = df[2].str.split(".",expand=True,)[0]
df_selection['type2'] = df_selection['type'].replace(['tFL','FL', 'FL1', 'FL2', 'FL3'],'tumor').replace(['GL'],'normal')
dict_sample = dict()

# Make a dictionary
for index, row in df_selection.iterrows():
    if row['Sample'] in dict_sample.keys():
        if row['type2'] == 'normal':
            dict_sample[row['Sample']]['normal'].append(row['Name'])
        elif row['type2'] == 'tumor':
            dict_sample[row['Sample']]['tumor'].append(row['Name'])
    else:
        tumor = list()
        normal = list()
        if row['type2'] == 'normal':
            normal.append(row['Name'])
        elif row['type2'] == 'tumor':
            tumor.append(row['Name'])
        dict_sample[row['Sample']] = {'normal': normal, 'tumor': tumor}

# write files
normal_tumor = pd.DataFrame(dict_sample).T.reset_index()
for index, row in normal_tumor.iterrows():
    argument = ''
    name_file = f'{row["index"]}_'
    output = f'{row["index"]}_'
    for tum in row['tumor']:
        number_tum = tum.replace("SS600", "")
        argument+=f'-I {head_path}{number_tum}/SN_{tum}.DR.bam '
        name_file+=f'{tum}_'
        output+=f'{tum}_'
    for norm in row['normal']:
        number_norm = norm.replace("SS600", "")
        argument+=f'-I {head_path}{number_norm}/SN_{norm}.DR.bam '
        argument+=f'-normal {norm}.DR.bam '
        name_file+=f'{norm}.txt'
        output+=f'{norm}_'
        for i in row['tumor']:
            number_i = i.replace("SS600", "")
            argument2 = ''
            name_file2 = f'{row["index"]}_'
            output2 = f'{row["index"]}_'
            argument2+=f'-I {head_path}{number_i}/SN_{i}.DR.bam '
            argument2+=f'-I {head_path}{number_norm}/SN_{norm}.DR.bam '
            argument2+=f'-normal {norm}.DR.bam '
            name_file2+=f'{i}_'            
            name_file2+=f'{norm}.txt'
            output2+=f'{i}_'
            output2+=f'{norm}_'
            #output2+=f'somatic.vcf.gz'
            argument2+=f'-O {head_path}{output2[:-1]}/unfiltered_{output2}somatic.vcf.gz'
            #print(argument2)
            #print(name_file2)
            f = open(name_file2, "w")
            f.write(argument2)
            f.write(f'\n{head_path}{output2[:-1]}/unfiltered_{output2}')
            f.write(f'\n{head_path}{output2[:-1]}/filtered_{output2}')
            f.close()
    
    #output+=f'somatic.vcf.gz'
    argument+=f'-O {head_path}{output[:-1]}/unfiltered_{output}somatic.vcf.gz'
    #print(argument)
    #print(name_file)
    f = open(name_file, "w")
    f.write(argument)
    f.write(f'\n{head_path}{output[:-1]}/unfiltered_{output}')
    f.write(f'\n{head_path}{output[:-1]}/filtered_{output}')
    f.close()
    
    
