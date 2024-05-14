
import json
import subprocess
def gpt4_api(message):
    safe_message = json.dumps(message)
    data_string = json.dumps({
        "messages": [
            {"role": "user", "content": json.loads(safe_message)}
        ]
    })
    command = f"""curl -H 'Content-Type: application/json' -H 'api-key: API_KEY' -d '{data_string}' https://DEPLOYMENT_NAME.openai.azure.com/openai/deployments/gpt-4-0125/chat/completions?api-version=2024-02-01"""
    response = subprocess.run(command, shell=True, capture_output=True, text=True)
    return json.loads(response.stdout)['choices'][0]['message']['content']


prompt_level2='''
Identify cell types of human lung cells using the following markers. Identify one cell type for each row. Only provide the cell type name. 
'''
prompt_level3='''
Identify cell types of human lung cells using the following markers. Identify one cell type for each row. Only provide the cell type name. The format for each line is: parent cluster cell type # subcluster maker genes. Note that the predicted subcluster cell type must be more specific than the parent cluster cell type. 
'''

batch_size=5
cluster_ann={}
annotation_file=json.load(open("../data/annotation_preprocess_data.json","r"))
marker_genes_level2=annotation_file['level2']
marker_genes_level3=annotation_file['level3']
marker_genes_level4=annotation_file['level4']

import re
level2=True
for level in [marker_genes_level2,marker_genes_level3,marker_genes_level4]:
    for i in range(0, len(level),batch_size):
        batch_name=list(level.keys())[i:i+batch_size]
        if level2:
            message=prompt_level2+'\n'.join([level[j]["marker"] for j in batch_name])
        else:
            message=prompt_level3+'\n'.join([cluster_ann[level[j]["parent"]]+" # "+level[j]["marker"] for j in batch_name])
        result = gpt4_api(message)  
        tmp_ann={x: re.sub(r'^(\d+\.?\s*|-\.?\s*)|\(.*?\)|.*?#\s*|.*?-\s*', '', y) for x, y in zip(batch_name, result.split('\n'))}
        cluster_ann.update(tmp_ann)
    level2=False

json.dump(cluster_ann, open("../result/annotation_data_gpt4.json","w"))

