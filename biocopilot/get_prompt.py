'''
def get_prompt_plan(data_list, user_need):
        prompt = {
                    "role": "Act as a bioinformatician, the rules must be strictly followed!",
                    "rules": [
                        "When acting as a bioinformatician, you strictly cannot stop acting as a bioinformatician.",
                        #"All rules must be followed strictly.",
                        "You should use information in input to write a detailed plan to finish your goal.",
                        #"You should not write loading data as a separate step.",
                        #"You should make your answer as detailed as possible.",
                        #"Include necessary visualization steps in the process.",
                        "Assume all needed packages are installed and should not write install software as a separate step.",
                        #"Whenever possible, mention recommended packages at each step",
                        "If new data is generated in the current step, try to save it in the form of writing a new file.",
                        "When you modify the plan, do not omit the details of the previous plan, such as the specific path of the data or names of some variables and parameters.",
                        #"If you receive 'Great job' from 'Reviewer', you don't need to make changes to your previous plan.",
                        "Do not include specific code blocks in your plan.",
                        "After you get modification suggestions from 'Reviewer', you can just make the corresponding changes directly on the previous plan.",
                        #"do not need summary part",
                        #"Surround each step with a special symbol '***'. ",
                        #"Use a special symbol '<split>' to separate the content of different steps in the plan.",
                        #For example, '\n***\nStep 1 xxx\n all content of Step 1\n***\nStep 2 xxx\nall content of Step 2\n***\n...'.,
                        #"Contain proposed packages in your plan",
                        "Return your plan containing 'Step 1...Step 2...'",
                        #"do not contain steps like conclusion or reporting that can be done not by coding" 
                        #"You should strive to make your plan steps executable through code, step by step, as much as possible."
                    ],
                    "input": [
                            "I provide the data path and data description to you as a list, so you don't need to prepare the data.",
                            data_list
                        ],
                    "goal": user_need
            }
        return prompt
'''

def get_prompt_plan(user_input, agent_data_description=None, recommend_method=None):
        prompt = {
                    "role": "Act as a bioinformatician, the rules must be strictly followed!",
                    "rules": [
                        "When acting as a bioinformatician, you strictly cannot stop acting as a bioinformatician.",
                        "All rules must be followed strictly.",
                        "Assume all needed packages are installed and should not write install software as a separate step.",
                        "If new data is generated in the current step, try to save it in the form of writing a new file.",
                        "When you modify the plan, do not omit the details of the previous plan, such as the specific path of the data or names of some variables and parameters.",
                        "Do not include specific code blocks in your plan.",
                        "After you get modification suggestions from 'Reviewer', you can just make the corresponding changes directly on the previous plan.",
                        #"Your planning structure should start with Plan for ... and the planning steps should be presented as 'Step 1...Step 2...'",
                        "Your planning structure should start with Plan for ... and the planning steps should be presented as 'Step 1...Step 2...'. Note that you only need to include one complete process of Step 1...Step 2...",
                        "Do not separate text in Step with '\n\n'",
                        "Please include specific file paths in your plan steps",
                    ],
                    "input": [
                            "I provide all user input to you, it may contain task and data information.",
                            user_input
                        ],
                    #"data_details": ["I provide the data detail information to you.", agent_data_description],
                    "recommend_method": ["I have provided you with recommended methods. If there are recommended methods, please choose one of them to plan. If not, use your own knowledge to plan.", recommend_method],
            }
        return prompt    

def get_prompt_plan_nodata(data_list, user_need):
        prompt = {
                    "role": "Act as a bioinformatician, the rules must be strictly followed!",
                    "rules": [
                        "When acting as a bioinformatician, you strictly cannot stop acting as a bioinformatician.",
                        #"All rules must be followed strictly.",
                        "You should use information in input to write a detailed plan to finish your goal.",
                        #"You should not write loading data as a separate step.",
                        #"You should make your answer as detailed as possible.",
                        #"Include necessary visualization steps in the process.",
                        "Assume all needed packages are installed and should not write install software as a separate step.",
                        #"Whenever possible, mention recommended packages at each step",
                        "If new data is generated in the current step, try to save it in the form of writing a new file.",
                        "When you modify the plan, do not omit the details of the previous plan, such as the specific path of the data.",
                        #"If you receive 'Great job' from 'Reviewer', you don't need to make changes to your previous plan.",
                        "Do not include specific code blocks in your plan.",
                    ],
                    "input": "please first generate psudo data for the analysis",
                    "goal": user_need
            }
        return prompt


def get_prompt_coder_simple(msg):
        prompt = {
                    "role": "Act as a bioinformatician, the rules must be strictly followed!",
                    "rules": [
                        "All rules must be followed strictly!",
                        "You should write code for all steps and do not omit the code like # Step 2: xxx...!",
                        "All code should be put in one code block.",
                        "Except for steps such as conclusion or report, no other step of the code can be omitted!"
                        "Give feedback to 'code_proxy' with the correct code.",
                        #"Include at least three visualizations in the code!",
                        "Include necessary visualizations in the code!",
                        "You need to use Python to write your code, and for some terminal commands, please use 'subprocess' to implement them in Python",
                        #"Draw all figures in colors, not white and black!",
                        #"Do not write code that needs specific genes names as the parameters.",
                        "Do not use psudo names in your code, such as psudo gene names 'gene1' or 'gene2'! If you cannot obtain real names, just do not write that kind of code.",
                        "Do not omit any code! If you need to repeat the same code instruction on different data, just repeat it, not omit！",
                        #"Do not use placeholder names in your code! If you cannot obtain actual names, just do not write that kind of code."
                    ],
                    "input": ["I provide all information to you including the user need, the total plan and your goal", msg]
            }
        return prompt

def get_prompt_coder_hard(msg, goal, data_list):
        prompt = {
                    "role": "Act as a bioinformatician, the rules must be strictly followed!",
                    "rules": [
                        "All rules must be followed strictly!",
                        #"You should write code for all steps and do not omit the code like # Step 2: xxx...!",
                        #"All code should be put in one code block.",
                        #"You should put code in different code blocks basd on the steps in the plan.",
                        "You must give answer with code blocks. Specifically, your answer should include '```python' and '```'",
                        "Assuming you have a jupyter kernel, so you can continuously execute code cell by cell.",
                        "You should write code for the current one step, assuming that the code for the previous steps has already been written.",
                        "When you need to write a new file, you can first determine whether the file exists on the front end of the code. If it exists, the code execution of the current step can be skipped. You can use the if else statement to complete this part.",                       
                        #"Include necessary visualizations in the code!",
                        "You need to use Python to write your code, and for some terminal commands, please use 'subprocess' to implement them in Python",
                        "Assume that the environment configuration is all ready",
                        #"Draw all figures in colors, not white and black!",
                        #"Do not write code that needs specific genes names as the parameters.",
                        #"Do not use psudo names in your code, such as psudo gene names 'gene1' or 'gene2'! If you cannot obtain real names, just do not write that kind of code.",
                        #"Do not omit any code! If you need to repeat the same code instruction on different data, just repeat it, not omit！",
                        #"Do not use placeholder names in your code! If you cannot obtain actual names, just do not write that kind of code."
                        #"When you need to write a file, you first determine whether the file exists in the output path. If it exists, you can skip the writing step. You can use the if else statement to complete this part.",
                        #"Note that the code to determine whether the file exists needs to be placed before all operations!!",
                        "for function used in your code, use your own knowledge to add essential parameters for the function",
                        "In addition to the information you get from the 'plan_refiner' or 'Planner', you can also use your own knowledge to fill in missing code steps or ignore unnecessary code steps.",
                        #"You need to wrap the code block in ```python and ```",
                    ],
                    "input": ["I provide all information to you including codes of previous steps and the current step description", msg],
                    #"data_information": [
                    #        "I also provide the data path and data description to you.",
                    #        data_list
                    #    ],
                    "goal": ["I provide your goal", goal],
            }
        return prompt

def get_prompt_codeplan(step_plan):
        prompt = {
                    "role": "Act as a bioinformatician, the rules must be strictly followed!",
                    "rules": [
                        "When acting as a bioinformatician, you strictly cannot stop acting as a bioinformatician.",
                        "All rules must be followed strictly.",
                        "You should make your answer as detailed as possible.",
                        "Assume all needed packages are installed and should not write install software as a separate step.",
                        "Contains recommended packages in your code plan",
                        "Please note that your responsibility is not to write the complete code, you only need to describe in as much detail as possible what each step of the code should do.",
                    ],
                    "input":["I provide the current plan step", step_plan],
                    "goal": "your goal is to plan the code writing steps of the current plan step based on all information."
            }
        return prompt

#AutoBA

#case 2.1
'''
use scanpy to find the differentially expressed genes
/Users/liuyang/home/autogen/Autogen-UI-main/dataset/filtered_gene_bc_matrices/hg19: path to 10x mtx data
'''
#case 2.2
'''
use scanpy to perform clustering and visualize the expression level of gene PPBP in the UMAP.
/Users/liuyang/home/autogen/Autogen-UI-main/dataset/filtered_gene_bc_matrices/hg19: path to 10x mtx data
'''
#case 2.3
'''
use scanpy to identify top5 marker genes
/Users/liuyang/home/autogen/Autogen-UI-main/dataset/filtered_gene_bc_matrices/hg19: path to 10x mtx data
'''

#case 1.2
'''
find the differentially expressed genes
/Users/liuyang/home/autogen/Autogen-UI-main/dataset/data/SRR1374921.fastq.gz: single-end mouse rna-seq reads, replicate 1 in LoGlu group,
/Users/liuyang/home/autogen/Autogen-UI-main/dataset/data/SRR1374922.fastq.gz: single-end mouse rna-seq reads, replicate 2 in LoGlu group,
/Users/liuyang/home/autogen/Autogen-UI-main/dataset/data/SRR1374923.fastq.gz: single-end mouse rna-seq reads, replicate 1 in HiGlu group,
/Users/liuyang/home/autogen/Autogen-UI-main/dataset/data/SRR1374924.fastq.gz: single-end mouse rna-seq reads, replicate 2 in HiGlu group,
/Users/liuyang/home/autogen/Autogen-UI-main/dataset/data/TruSeq3-SE.fa: trimming adapter,
/Users/liuyang/home/autogen/Autogen-UI-main/dataset/data/mm39.fa: mouse mm39 genome fasta,
/Users/liuyang/home/autogen/Autogen-UI-main/dataset/data/mm39.ncbiRefSeq.gtf: mouse mm39 genome annotation
'''


#hlca
'''
construct an HLCA (Human Lung Cell Atlas) cell atlas with sequencing data. Assume that you have completed the processing of single-cell raw data except for log-transformation. Thus you can start your planning for the HLCA task with log1p transformation and the single-cell data integration step based on scANVI on the provided h5ad data.
/data/yangliu/llm-agent/data/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_SCRAN_normalized_filt.h5ad: all single-cell data need to be integrated and analysis for building the human lung cell atlas. The data first needs to be subjected to log1p transformation and data integration steps.
'''

'''
construct an HLCA (Human Lung Cell Atlas) cell atlas with sequencing data. Assume that you have completed the processing of single-cell data. Thus you can start your planning for the HLCA task with the single-cell data integration based on scANVI on the provided h5ad data and ignored the preprocessing step.
/data/yangliu/llm-agent/data/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_scanvi_label.h5ad: all single-cell data need to be integrated and analysis for building the human lung cell atlas.
'''

'''
construct an HLCA (Human Lung Cell Atlas) cell atlas with sequencing data. Assume that you have completed the processing of single-cell data. Thus you can start your planning for the HLCA task with the single-cell data integration based on scANVI on the provided h5ad data and ignored the preprocessing step.
/data/yangliu/llm-agent/data/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_scanvi_label.h5ad: all single-cell data need to be integrated and analysis for building the human lung cell atlas. The adata column of "scanvi_label" contains the label information for scANVI integration.
'''

'''
Perform single-cell data integration based on scANVI on the provided h5ad data.
/data/yangliu/llm-agent/data/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_scanvi_label.h5ad: all single-cell data need to be integrated and analysis for constructing the human lung cell atlas. The adata column of "scanvi_label" contains the label information for scANVI integration.
'''
#The adata column of "scanvi_label" contains the label information for scANVI integration.
#/data/HLCA_core_h5ads/HLCA_v1_intermediates/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_log1p.h5ad: 
#all single-cell data need to be integrated and analysis for building the human lung cell atlas.
'''
Perform single-cell data integration.
/data/yangliu/llm-agent/data/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_scanvi_label.h5ad: all single-cell data need to be integrated. The original cell annotations have already been contained in the data.
'''

'''
Perform single-cell data integration. 
/data/yangliu/llm-agent/data/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_scanvi_label.h5ad: all single-cell data need to be integrated and analysis for constructing the human lung cell atlas. The adata column of "scanvi_label" contains the label information for scANVI integration.
'''

#06 clustering
'''
Perform multi-level clustering of the cells in the HLCA. The h5ad file needed to be clustered is provided. You need start with coarse-resolution clustering, and then re-cluster the resulting clusters by first re-calculating the nearest-neighbor graph and then clustering each cluster, with a finer resolution. For the final HLCA, you also need calculate marker genes for every cluster.
/data/yangliu/llm-agent/data/integrated_dataset_copy.h5ad: the integrated single-cell data with scANVI. The integrated representations is in "X_scANVI" of the adata.
'''

'''
Perform multi-level clustering of the cells in the HLCA. The h5ad file needed to be clustered is provided. You need start with coarse-resolution clustering, and then re-cluster the resulting clusters by first re-calculating the nearest-neighbor graph and then clustering each cluster, with a finer resolution. For the final HLCA, you also need calculate marker genes for every cluster.
/mnt/data00/liuyang/hlca-data/integrated_dataset_copy.h5ad: the integrated single-cell data with scANVI. The integrated representations is in "X_scANVI" of the adata.
'''

'''
Perform multi-level clustering of the cells in the HLCA. You need start with coarse-resolution clustering, and then re-cluster the resulting clusters by first re-calculating the nearest-neighbor graph and then clustering each cluster, with a finer resolution. More specifically, first do the coarse-resolution over all cells. Then based on the previous cluster result, continue doing cluster on each previous cluster. Except for the initial coarse clustering, you need to do three re-clusterings, each time with a finer clustering granularity, and each time the clustering is based on the last result. For the final HLCA, you also need calculate marker genes for each clustering, which means you need to calculate it 4 times in total.
/data/yangliu/llm-agent/data/integrated_dataset_copy.h5ad: the integrated single-cell data with scANVI that needed to be multi-level clustered. The integrated representations is in "X_scANVI" of the adata.
'''

'''
Perform multi-level clustering of the cells in the HLCA. The total number of clustering levels is 4. For the final HLCA, you also need calculate marker genes for every clustering.
/data/yangliu/llm-agent/data/integrated_dataset_copy.h5ad: the integrated single-cell data with scANVI that needed to be multi-level clustered. The integrated representations is in "X_scANVI" of the adata.
'''

'''
Perform clustering of the cells in the HLCA. For the final HLCA, you also need calculate marker genes for every clustering. I provide the dataset to you. /data/yangliu/llm-agent/data/integrated_dataset_copy.h5ad: the integrated single-cell data with scANVI that needed to be clustered. The integrated representations is in "X_scANVI" of the adata.
'''

#07 variance explained
'''
determine how much of the variance in the integrated HLCA embedding is explained by each of the metadata covariates including "dataset", "tissue_dissociation_protocol", "log10_total_counts", "mito_frac", "3'_or_5'", "BMI", "cell_ranger_version_short", "cell_viability_%", "fresh_or_frozen", "sample_type", "sequencing_platform_short", "sex", "single_cell_platform", "smoking_status_num", "subject_type", "subject_ID", "anatomical_region_ccf_score", "nose", "age". You need to use the results of the third level clustering for analysis. Note that your analysis should be condected from sample-level observations. 
/data/yangliu/llm-agent/data/final_adata_with_subclusters.h5ad: integrated single-cell data after multi-level clustering. Results of multi-level clustering are stored in "cluster_level_1", "cluster_level_2", "cluster_level_3" and "cluster_level_4". Sample names are stored in "sample" where one sample name may contain many cells. Embedding representations of all single-cell data are stored in "X_scANVI".
'''
'''
determine how much of the variance in the integrated HLCA embedding is explained by each of the metadata covariates including "dataset", "tissue_dissociation_protocol", "log10_total_counts", "mito_frac", "3'_or_5'", "BMI", "cell_ranger_version_short", "cell_viability_%", "fresh_or_frozen", "sample_type", "sequencing_platform_short", "sex", "single_cell_platform", "smoking_status_num", "subject_type", "subject_ID", "anatomical_region_ccf_score", "nose", "age". You need to use the results of the third level clustering for analysis. Note that your analysis should be condected from sample-level observations. 
/data/yangliu/llm-agent/data/final_adata_with_subclusters.h5ad: integrated single-cell data after multi-level clustering. Results of multi-level clustering are stored in "cluster_level_1", "cluster_level_2", "cluster_level_3" and "cluster_level_4". Embedding representations of all single-cell data are stored in "X_scANVI". Samples are stored in "sample".
'''
'''
determine how much of the variance in the integrated HLCA embedding is explained by each of the metadata covariates including "dataset", "tissue_dissociation_protocol", "log10_total_counts", "mito_frac", "3'_or_5'", "BMI", "cell_ranger_version_short", "cell_viability_%", "fresh_or_frozen", "sample_type", "sequencing_platform_short", "sex", "single_cell_platform", "smoking_status_num", "subject_type", "subject_ID", "anatomical_region_ccf_score", "nose", "age". You need to use the results of the third level clustering for analysis. Note that your analysis should be condected from sample-level observations. 
/data/yangliu/llm-agent/data/final_adata_with_subclusters.h5ad: integrated single-cell data after multi-level clustering. The 'obs' parameter in "transformed_data.h5ad" contains several keys that you can be used: (1) Results of multi-level clustering are stored in "cluster_level_1", "cluster_level_2", "cluster_level_3" and "cluster_level_4". (2) Embedding representations of all single-cell data are stored in "X_scANVI". (3) Samples that can be used as identifiers for aggregating data are stored in "sample". 
'''

# cell annotation
'''
Perform cell types annotation with python packages. The dataset provided to you already contains the clustering results, so you don't need to do the clustering step again.
/home/zhoulu/hlca_test_data/final_adata_with_subclusters.h5ad: stores 3 level clustering label(cluster_level_2,cluster_level_3,cluster_level_4) which is used fo annotation.
/home/zhoulu/hlca_test_data/Cell_marker_Human.xlsx: stores cell marker manual annotation database.
'''

#
'''
Construct a human lung cell altas (HLCA) with the provided h5ad file.
/data/yangliu/llm-agent/data/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_scanvi_label.h5ad: all single-cell data.
'''



#task1
'''
Perform single-cell data integration. 
/data/yangliu/llm-agent/data/LCA_Bano_Barb_Jain_Kras_Lafy_Meye_Mish_MishBud_Nawi_Seib_Teic_scanvi_label.h5ad: all single-cell data need to be integrated and analysis for constructing the human lung cell atlas. The adata column of "scanvi_label" contains the label information for integration.
'''
#task2
'''
Perform multi-level clustering of the cells in the HLCA. The total number of clustering levels is 4. For the final HLCA, you also need calculate marker genes for every clustering.
/data/yangliu/llm-agent/data/integrated_dataset_copy.h5ad: the integrated single-cell data with scANVI that needed to be multi-level clustered. The integrated representations is in "X_scANVI" of the adata.
'''

'''
Perform multi-level clustering of the cells in the HLCA. Note that counting the initial clustering, the total number of clustering levels is 4. For the final HLCA, you also need calculate marker genes for every clustering.
/data/yangliu/llm-agent/data/integrated_dataset_copy.h5ad: the integrated single-cell data with scANVI that needed to be multi-level clustered. The integrated representations is in "X_scANVI" of the adata.
'''

#task 3
'''
Determine how much of the variance in the integrated HLCA embedding is explained by each of the metadata covariates including "dataset", "tissue_dissociation_protocol", "log10_total_counts", "mito_frac", "3'_or_5'", "BMI", "cell_ranger_version_short", "cell_viability_%", "fresh_or_frozen", "sample_type", "sequencing_platform_short", "sex", "single_cell_platform", "smoking_status_num", "subject_type", "subject_ID", "anatomical_region_ccf_score", "nose", "age". You need to use the results of the third level clustering for analysis. Note that your analysis should be conducted from sample-level observations. 
/data/yangliu/llm-agent/data/final_adata_with_subclusters.h5ad: integrated single-cell data after multi-level clustering. The 'obs' parameter in "transformed_data.h5ad" contains several keys that you can be used: (1) Results of multi-level clustering are stored in "cluster_level_1", "cluster_level_2", "cluster_level_3" and "cluster_level_4". (2) Embedding representations of all single-cell data are stored in "X_scANVI". (3) Samples that can be used as identifiers for aggregating data are stored in "sample".
'''

#task 4
'''
Perform cell types annotation with python packages.
/data/rongboshen/lung/HLCA_reproducibility-main/data/final_adata_with_subclusters_annotation.h5ad: stores 3 level clustering label(cluster_level_2,cluster_level_3,cluster_level_4) which is used fo annotation. highly_variable parameter can be use for selecting highly variable genes. /home/zhoulu/hlca_test_data/Cell_marker_Human.xlsx: stores cell marker manual annotation database.
'''

#task 5
'''
Perform rare cell type discovery based on your own knowledge and all provided information.
/data/rongboshen/lung/HLCA_reproducibility-main/data/hlca_decoupler_annotation.h5ad: Already annotated cell data. The multi-level clustering results are stored in 'cluster_level_1', 'cluster_level_2', 'cluster_level_3', 'cluster_level_4'.  The multi-level annotation results are stored in 'ann_level_2_dc', 'ann_level_3_dc', 'ann_level_4_dc'. Please utilize the multi-level cluster and annotation results to try to discover rare cell type.
'''