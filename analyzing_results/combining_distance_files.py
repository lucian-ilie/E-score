import pandas as pd
import os


def main():

    path = "path/to/folder/with/multiple/distance.csv/files"

    files = os.listdir(str(path))
    first_completed = False

    for element in files:

        if element.endswith(".csv"):


            if first_completed == True:

                MSA = element.split("_")[0]
                data_2 = pd.read_csv(path + element)

                data_2 = data_2.rename(columns={'SSP': str(MSA) + ' SSP'})
                data_2 = data_2.rename(columns={'SEQ': str(MSA) + ' SEQ'})
                data_2 = data_2.rename(columns={'POS': str(MSA) + ' POS'})
                data_2 = data_2.rename(columns={'DD': str(MSA) + ' DD'})
                data_2 = data_2.rename(columns={'CC': str(MSA) + ' CC'})

                selected_columns2 = data_2.iloc[:, 1:6]
                selected_columns2.reset_index(drop=True, inplace=True)

                selected_columns = pd.concat([selected_columns, selected_columns2], axis=1)
                

            if first_completed == False:
                first_completed = True
                
        
                MSA = element.split("_")[0]
                Shift_Down = element.split("_")[2]
                Gap_Open = element.split("_")[3]
                Gap_Extension = element.split("_")[4]

                data_1 = pd.read_csv(str(path) +"/" + str(element))

                data_1 = data_1.rename(columns={'SSP': str(MSA) + ' SSP'})
                data_1 = data_1.rename(columns={'SEQ': str(MSA) + ' SEQ'})
                data_1 = data_1.rename(columns={'POS': str(MSA) + ' POS'})
                data_1 = data_1.rename(columns={'DD': str(MSA) + ' DD'})
                data_1 = data_1.rename(columns={'CC': str(MSA) + ' CC'})

                selected_columns = data_1.iloc[:, 1:7]
        

                #create the new columns with constant values
                selected_columns.insert(0, 'Shift', Shift_Down)
                selected_columns.insert(1, 'Gap Open', Gap_Open)
                selected_columns.insert(2, 'Gap Extension', Gap_Extension)


        #un-comment which one you need, or define another order
        #CDD 8 domains
        reorder = ['Shift', 'Gap Open', 'Gap Extension', 
        'cd00012 SSP', 'cd00012 SEQ', 'cd00012 POS', 'cd00012 DD', 'cd00012 CC',
        'cd00083 SSP', 'cd00083 SEQ', 'cd00083 POS', 'cd00083 DD', 'cd00083 CC', 
        'cd00173 SSP', 'cd00173 SEQ', 'cd00173 POS', 'cd00173 DD', 'cd00173 CC',
        'cd00637 SSP', 'cd00637 SEQ', 'cd00637 POS', 'cd00637 DD', 'cd00637 CC',
        'cd01040 SSP', 'cd01040 SEQ', 'cd01040 POS', 'cd01040 DD', 'cd01040 CC', 
        'cd01068 SSP', 'cd01068 SEQ', 'cd01068 POS', 'cd01068 DD', 'cd01068 CC', 
        'cd10140 SSP', 'cd10140 SEQ', 'cd10140 POS', 'cd10140 DD', 'cd10140 CC',
        'cd21116 SSP', 'cd21116 SEQ', 'cd21116 POS', 'cd21116 DD', 'cd21116 CC' ]

        # #balibase natural RV11
        # reorder = ['Shift', 'Gap Open', 'Gap Extension',
        # 'BBS11009 SSP', 'BBS11009 SEQ', 'BBS11009 POS', 'BBS11009 DD', 'BBS11009 CC',
        # 'BBS11010 SSP', 'BBS11010 SEQ', 'BBS11010 POS', 'BBS11010 DD', 'BBS11010 CC',
        # 'BBS11016 SSP', 'BBS11016 SEQ', 'BBS11016 POS', 'BBS11016 DD', 'BBS11016 CC',
        # 'BBS11018 SSP', 'BBS11018 SEQ', 'BBS11018 POS', 'BBS11018 DD', 'BBS11018 CC',
        # 'BBS11024 SSP', 'BBS11024 SEQ', 'BBS11024 POS', 'BBS11024 DD', 'BBS11024 CC',
        # 'BBS11026 SSP', 'BBS11026 SEQ', 'BBS11026 POS', 'BBS11026 DD', 'BBS11026 CC',
        # 'BBS11034 SSP', 'BBS11034 SEQ', 'BBS11034 POS', 'BBS11034 DD', 'BBS11034 CC',
        # 'BBS11037 SSP', 'BBS11037 SEQ', 'BBS11037 POS', 'BBS11037 DD', 'BBS11037 CC',
        # 'BBS11038 SSP', 'BBS11038 SEQ', 'BBS11038 POS', 'BBS11038 DD', 'BBS11038 CC']


        # #balibase natural RV12
        # reorder = ['Shift', 'Gap Open', 'Gap Extension',
        # 'BBS12008 SSP', 'BBS12008 SEQ', 'BBS12008 POS', 'BBS12008 DD', 'BBS12008 CC',
        # 'BBS12011 SSP', 'BBS12011 SEQ', 'BBS12011 POS', 'BBS12011 DD', 'BBS12011 CC',
        # 'BBS12030 SSP', 'BBS12030 SEQ', 'BBS12030 POS', 'BBS12030 DD', 'BBS12030 CC',
        # 'BBS12033 SSP', 'BBS12033 SEQ', 'BBS12033 POS', 'BBS12033 DD', 'BBS12033 CC',
        # 'BBS12037 SSP', 'BBS12037 SEQ', 'BBS12037 POS', 'BBS12037 DD', 'BBS12037 CC']

        # #balibase natural RV30
        # reorder = ['Shift', 'Gap Open', 'Gap Extension',
        # 'BBS30001 SSP', 'BBS30001 SEQ', 'BBS30001 POS', 'BBS30001 DD', 'BBS30001 CC',
        # 'BBS30004 SSP', 'BBS30004 SEQ', 'BBS30004 POS', 'BBS30004 DD', 'BBS30004 CC',
        # 'BBS30005 SSP', 'BBS30005 SEQ', 'BBS30005 POS', 'BBS30005 DD', 'BBS30005 CC',
        # 'BBS30006 SSP', 'BBS30006 SEQ', 'BBS30006 POS', 'BBS30006 DD', 'BBS30006 CC',
        # 'BBS30008 SSP', 'BBS30008 SEQ', 'BBS30008 POS', 'BBS30008 DD', 'BBS30008 CC',
        # 'BBS30009 SSP', 'BBS30009 SEQ', 'BBS30009 POS', 'BBS30009 DD', 'BBS30009 CC',
        # 'BBS30010 SSP', 'BBS30010 SEQ', 'BBS30010 POS', 'BBS30010 DD', 'BBS30010 CC',
        # 'BBS30013 SSP', 'BBS30013 SEQ', 'BBS30013 POS', 'BBS30013 DD', 'BBS30013 CC',
        # 'BBS30016 SSP', 'BBS30016 SEQ', 'BBS30016 POS', 'BBS30016 DD', 'BBS30016 CC',
        # 'BBS30030 SSP', 'BBS30030 SEQ', 'BBS30030 POS', 'BBS30030 DD', 'BBS30030 CC']

        # #gpcr
        # reorder = ['Shift', 'Gap Open', 'Gap Extension', 
        # 'alicarboxylic SSP', 'alicarboxylic SEQ', 'alicarboxylic POS', 'alicarboxylic DD', 'alicarboxylic CC',
        # 'aminergic SSP', 'aminergic SEQ', 'aminergic POS', 'aminergic DD', 'aminergic CC',
        # 'lipid SSP', 'lipid SEQ', 'lipid POS', 'lipid DD', 'lipid CC',
        # 'nucleotide SSP', 'nucleotide SEQ', 'nucleotide POS', 'nucleotide DD', 'nucleotide CC', 
        # 'orphan SSP', 'orphan SEQ', 'orphan POS', 'orphan DD', 'orphan CC',
        # 'peptide SSP', 'peptide SEQ', 'peptide POS', 'peptide DD', 'peptide CC',
        # 'protein SSP', 'protein SEQ', 'protein POS', 'protein DD', 'protein CC',
        # 'sensory SSP', 'sensory SEQ', 'sensory POS', 'sensory DD', 'sensory CC']


        # #RV911
        # reorder =  ['Shift', 'Gap Open', 'Gap Extension',
        # 'BOSX001 SSP', 'BOSX001 SEQ', 'BOSX001 POS', 'BOSX001 DD', 'BOSX001 CC',
        # 'BOSX022 SSP', 'BOSX022 SEQ', 'BOSX022 POS', 'BOSX022 DD', 'BOSX022 CC',
        # 'BOSX032 SSP', 'BOSX032 SEQ', 'BOSX032 POS', 'BOSX032 DD', 'BOSX032 CC',
        # 'BOSX034 SSP', 'BOSX034 SEQ', 'BOSX034 POS', 'BOSX034 DD', 'BOSX034 CC',
        # 'BOSX035 SSP', 'BOSX035 SEQ', 'BOSX035 POS', 'BOSX035 DD', 'BOSX035 CC',
        # 'BOSX046 SSP', 'BOSX046 SEQ', 'BOSX046 POS', 'BOSX046 DD', 'BOSX046 CC',
        # 'BOSX060 SSP', 'BOSX060 SEQ', 'BOSX060 POS', 'BOSX060 DD', 'BOSX060 CC',
        # 'BOSX063 SSP', 'BOSX063 SEQ', 'BOSX063 POS', 'BOSX063 DD', 'BOSX063 CC',
        # 'BOSX076 SSP', 'BOSX076 SEQ', 'BOSX076 POS', 'BOSX076 DD', 'BOSX076 CC',
        # 'BOSX096 SSP', 'BOSX096 SEQ', 'BOSX096 POS', 'BOSX096 DD', 'BOSX096 CC',
        # 'BOSX115 SSP', 'BOSX115 SEQ', 'BOSX115 POS', 'BOSX115 DD', 'BOSX115 CC',
        # 'BOSX121 SSP', 'BOSX121 SEQ', 'BOSX121 POS', 'BOSX121 DD', 'BOSX121 CC',
        # 'BOSX122 SSP', 'BOSX122 SEQ', 'BOSX122 POS', 'BOSX122 DD', 'BOSX122 CC',
        # 'BOSX123 SSP', 'BOSX123 SEQ', 'BOSX123 POS', 'BOSX123 DD', 'BOSX123 CC',
        # 'BOSX142 SSP', 'BOSX142 SEQ', 'BOSX142 POS', 'BOSX142 DD', 'BOSX142 CC',
        # 'BOSX172 SSP', 'BOSX172 SEQ', 'BOSX172 POS', 'BOSX172 DD', 'BOSX172 CC',
        # 'BOSX175 SSP', 'BOSX175 SEQ', 'BOSX175 POS', 'BOSX175 DD', 'BOSX175 CC',
        # 'BOSX177 SSP', 'BOSX177 SEQ', 'BOSX177 POS', 'BOSX177 DD', 'BOSX177 CC',
        # 'BOSX180 SSP', 'BOSX180 SEQ', 'BOSX180 POS', 'BOSX180 DD', 'BOSX180 CC',
        # 'BOSX181 SSP', 'BOSX181 SEQ', 'BOSX181 POS', 'BOSX181 DD', 'BOSX181 CC',
        # 'BOSX192 SSP', 'BOSX192 SEQ', 'BOSX192 POS', 'BOSX192 DD', 'BOSX192 CC',
        # 'BOSX212 SSP', 'BOSX212 SEQ', 'BOSX212 POS', 'BOSX212 DD', 'BOSX212 CC',
        # 'BOSX214 SSP', 'BOSX214 SEQ', 'BOSX214 POS', 'BOSX214 DD', 'BOSX214 CC',
        # 'BOSX222 SSP', 'BOSX222 SEQ', 'BOSX222 POS', 'BOSX222 DD', 'BOSX222 CC',
        # 'BOSX240 SSP', 'BOSX240 SEQ', 'BOSX240 POS', 'BOSX240 DD', 'BOSX240 CC',
        # 'BOSX246 SSP', 'BOSX246 SEQ', 'BOSX246 POS', 'BOSX246 DD', 'BOSX246 CC',
        # 'BOSX258 SSP', 'BOSX258 SEQ', 'BOSX258 POS', 'BOSX258 DD', 'BOSX258 CC',
        # 'BOSX270 SSP', 'BOSX270 SEQ', 'BOSX270 POS', 'BOSX270 DD', 'BOSX270 CC',
        # 'BOSX284 SSP', 'BOSX284 SEQ', 'BOSX284 POS', 'BOSX284 DD', 'BOSX284 CC']

        # #RV912
        # reorder =  ['Shift', 'Gap Open', 'Gap Extension',
        # 'BOSX011 SSP', 'BOSX011 SEQ', 'BOSX011 POS', 'BOSX011 DD', 'BOSX011 CC',
        # 'BOSX032 SSP', 'BOSX032 SEQ', 'BOSX032 POS', 'BOSX032 DD', 'BOSX032 CC',
        # 'BOSX045 SSP', 'BOSX045 SEQ', 'BOSX045 POS', 'BOSX045 DD', 'BOSX045 CC',
        # 'BOSX047 SSP', 'BOSX047 SEQ', 'BOSX047 POS', 'BOSX047 DD', 'BOSX047 CC',
        # 'BOSX049 SSP', 'BOSX049 SEQ', 'BOSX049 POS', 'BOSX049 DD', 'BOSX049 CC',
        # 'BOSX050 SSP', 'BOSX050 SEQ', 'BOSX050 POS', 'BOSX050 DD', 'BOSX050 CC',
        # 'BOSX054 SSP', 'BOSX054 SEQ', 'BOSX054 POS', 'BOSX054 DD', 'BOSX054 CC',
        # 'BOSX060 SSP', 'BOSX060 SEQ', 'BOSX060 POS', 'BOSX060 DD', 'BOSX060 CC',
        # 'BOSX075 SSP', 'BOSX075 SEQ', 'BOSX075 POS', 'BOSX075 DD', 'BOSX075 CC',
        # 'BOSX076 SSP', 'BOSX076 SEQ', 'BOSX076 POS', 'BOSX076 DD', 'BOSX076 CC',
        # 'BOSX096 SSP', 'BOSX096 SEQ', 'BOSX096 POS', 'BOSX096 DD', 'BOSX096 CC',
        # 'BOSX121 SSP', 'BOSX121 SEQ', 'BOSX121 POS', 'BOSX121 DD', 'BOSX121 CC',
        # 'BOSX122 SSP', 'BOSX122 SEQ', 'BOSX122 POS', 'BOSX122 DD', 'BOSX122 CC',
        # 'BOSX142 SSP', 'BOSX142 SEQ', 'BOSX142 POS', 'BOSX142 DD', 'BOSX142 CC',
        # 'BOSX149 SSP', 'BOSX149 SEQ', 'BOSX149 POS', 'BOSX149 DD', 'BOSX149 CC',
        # 'BOSX154 SSP', 'BOSX154 SEQ', 'BOSX154 POS', 'BOSX154 DD', 'BOSX154 CC',
        # 'BOSX175 SSP', 'BOSX175 SEQ', 'BOSX175 POS', 'BOSX175 DD', 'BOSX175 CC',
        # 'BOSX177 SSP', 'BOSX177 SEQ', 'BOSX177 POS', 'BOSX177 DD', 'BOSX177 CC',
        # 'BOSX187 SSP', 'BOSX187 SEQ', 'BOSX187 POS', 'BOSX187 DD', 'BOSX187 CC',
        # 'BOSX192 SSP', 'BOSX192 SEQ', 'BOSX192 POS', 'BOSX192 DD', 'BOSX192 CC',
        # 'BOSX202 SSP', 'BOSX202 SEQ', 'BOSX202 POS', 'BOSX202 DD', 'BOSX202 CC',
        # 'BOSX214 SSP', 'BOSX214 SEQ', 'BOSX214 POS', 'BOSX214 DD', 'BOSX214 CC',
        # 'BOSX240 SSP', 'BOSX240 SEQ', 'BOSX240 POS', 'BOSX240 DD', 'BOSX240 CC',
        # 'BOSX246 SSP', 'BOSX246 SEQ', 'BOSX246 POS', 'BOSX246 DD', 'BOSX246 CC',
        # 'BOSX258 SSP', 'BOSX258 SEQ', 'BOSX258 POS', 'BOSX258 DD', 'BOSX258 CC',
        # 'BOSX259 SSP', 'BOSX259 SEQ', 'BOSX259 POS', 'BOSX259 DD', 'BOSX259 CC',
        # 'BOSX281 SSP', 'BOSX281 SEQ', 'BOSX281 POS', 'BOSX281 DD', 'BOSX281 CC',
        # 'BOSX290 SSP', 'BOSX290 SEQ', 'BOSX290 POS', 'BOSX290 DD', 'BOSX290 CC']

        # #RV913
        # reorder =  ['Shift', 'Gap Open', 'Gap Extension',
        # 'BOSX012 SSP', 'BOSX012 SEQ', 'BOSX012 POS', 'BOSX012 DD', 'BOSX012 CC',
        # 'BOSX017 SSP', 'BOSX017 SEQ', 'BOSX017 POS', 'BOSX017 DD', 'BOSX017 CC',
        # 'BOSX032 SSP', 'BOSX032 SEQ', 'BOSX032 POS', 'BOSX032 DD', 'BOSX032 CC',
        # 'BOSX034 SSP', 'BOSX034 SEQ', 'BOSX034 POS', 'BOSX034 DD', 'BOSX034 CC',
        # 'BOSX036 SSP', 'BOSX036 SEQ', 'BOSX036 POS', 'BOSX036 DD', 'BOSX036 CC',
        # 'BOSX043 SSP', 'BOSX043 SEQ', 'BOSX043 POS', 'BOSX043 DD', 'BOSX043 CC',
        # 'BOSX049 SSP', 'BOSX049 SEQ', 'BOSX049 POS', 'BOSX049 DD', 'BOSX049 CC',
        # 'BOSX063 SSP', 'BOSX063 SEQ', 'BOSX063 POS', 'BOSX063 DD', 'BOSX063 CC',
        # 'BOSX075 SSP', 'BOSX075 SEQ', 'BOSX075 POS', 'BOSX075 DD', 'BOSX075 CC',
        # 'BOSX076 SSP', 'BOSX076 SEQ', 'BOSX076 POS', 'BOSX076 DD', 'BOSX076 CC',
        # 'BOSX079 SSP', 'BOSX079 SEQ', 'BOSX079 POS', 'BOSX079 DD', 'BOSX079 CC',
        # 'BOSX082 SSP', 'BOSX082 SEQ', 'BOSX082 POS', 'BOSX082 DD', 'BOSX082 CC',
        # 'BOSX121 SSP', 'BOSX121 SEQ', 'BOSX121 POS', 'BOSX121 DD', 'BOSX121 CC',
        # 'BOSX122 SSP', 'BOSX122 SEQ', 'BOSX122 POS', 'BOSX122 DD', 'BOSX122 CC',
        # 'BOSX126 SSP', 'BOSX126 SEQ', 'BOSX126 POS', 'BOSX126 DD', 'BOSX126 CC',
        # 'BOSX132 SSP', 'BOSX132 SEQ', 'BOSX132 POS', 'BOSX132 DD', 'BOSX132 CC',
        # 'BOSX133 SSP', 'BOSX133 SEQ', 'BOSX133 POS', 'BOSX133 DD', 'BOSX133 CC',
        # 'BOSX142 SSP', 'BOSX142 SEQ', 'BOSX142 POS', 'BOSX142 DD', 'BOSX142 CC',
        # 'BOSX146 SSP', 'BOSX146 SEQ', 'BOSX146 POS', 'BOSX146 DD', 'BOSX146 CC',
        # 'BOSX158 SSP', 'BOSX158 SEQ', 'BOSX158 POS', 'BOSX158 DD', 'BOSX158 CC',
        # 'BOSX180 SSP', 'BOSX180 SEQ', 'BOSX180 POS', 'BOSX180 DD', 'BOSX180 CC',
        # 'BOSX183 SSP', 'BOSX183 SEQ', 'BOSX183 POS', 'BOSX183 DD', 'BOSX183 CC',
        # 'BOSX214 SSP', 'BOSX214 SEQ', 'BOSX214 POS', 'BOSX214 DD', 'BOSX214 CC',
        # 'BOSX222 SSP', 'BOSX222 SEQ', 'BOSX222 POS', 'BOSX222 DD', 'BOSX222 CC',
        # 'BOSX246 SSP', 'BOSX246 SEQ', 'BOSX246 POS', 'BOSX246 DD', 'BOSX246 CC',
        # 'BOSX258 SSP', 'BOSX258 SEQ', 'BOSX258 POS', 'BOSX258 DD', 'BOSX258 CC',
        # 'BOSX290 SSP', 'BOSX290 SEQ', 'BOSX290 POS', 'BOSX290 DD', 'BOSX290 CC']



    selected_columns = selected_columns[reorder]
    display(selected_columns)

    directory = "path/to/where/the/file/should/be/saved/to"
    file_name = str(Gap_Open) + "_" + str(Gap_Extension) + "_single_test.csv"

    #file_name = str(Gap_Open) + "_" + str(Gap_Extension) + "_single_test_blosum.csv"

    print(file_name)
    file_path = os.path.join(directory, file_name)
    selected_columns.to_csv(file_name, index=True)



if __name__ == "__main__":
    main()
