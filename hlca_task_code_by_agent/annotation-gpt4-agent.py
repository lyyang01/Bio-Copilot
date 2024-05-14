
# The core code for cell types annotation using gpt4 api.

### Step 1: Parse the Annotation Preprocess Data
import json
import os
import subprocess

# Define the path to the JSON file
json_file_path = '../data/annotation_preprocess_data.json'

# Check if the file exists before proceeding
if os.path.exists(json_file_path):
    try:
        # Load the JSON data
        with open(json_file_path, 'r') as file:
            annotation_data = json.load(file)
        
        # Function to extract hierarchical structure
        def extract_hierarchy(data, levels=['level2', 'level3', 'level4']):
            """
            Extracts hierarchical structure from the annotation data.
            
            Parameters:
            - data: The loaded JSON data.
            - levels: The levels of annotation to extract.
            
            Returns:
            A dictionary with keys as levels and values as lists of dictionaries containing
            'cluster', 'marker', and 'parent' information.
            """
            extracted_data = {}
            for level in levels:
                level_data = data.get(level, {})
                extracted_data[level] = [
                    {'cluster': cluster, 'marker': details['marker'], 'parent': details['parent']}
                    for cluster, details in level_data.items()
                ]
            return extracted_data
        
        # Extract hierarchical structure
        hierarchical_data = extract_hierarchy(annotation_data)
        
        # Data structuring for further use
        # Define a simple data model as nested dictionaries
        structured_data = {}
        for level, clusters in hierarchical_data.items():
            structured_data[level] = {}
            for cluster_info in clusters:
                cluster_name = cluster_info['cluster']
                structured_data[level][cluster_name] = {
                    'marker': cluster_info['marker'],
                    'parent': cluster_info['parent']
                }
        
        # Documentation and comments are included within the code as inline comments and docstrings
        
    except Exception as e:
        print(f"An error occurred while processing the file: {e}")
else:
    print(f"File does not exist at the specified path: {json_file_path}")


### Step 2: Level 2 Cell Type Annotation
import requests
import json

# Assuming the structured_data dictionary from the previous step is available and contains the 'level2' data
level2_data = structured_data['level2']

# Initialize an empty dictionary to store the cluster name and its predicted cell type
predicted_cell_types = {}

# Prepare batches of marker genes for batch processing
marker_genes_batches = [list(level2_data.values())[i:i+10] for i in range(0, len(level2_data), 10)]

# Set the prompt for GPT-4 in level 2
prompt_base = "Identify cell types of human lung cells using the following markers. Identify one cell type for each row. Only provide the cell type name."

# API information
API_KEY = "YOUR_API_KEY" # change to your api key
RESOURCE_NAME = "DEPLOYMENT_NAME" # change to your deployment name on azure
DEPLOYMENT_NAME = "gpt-4-0125"
API_VERSION = "2024-02-01"

# Process each batch
for batch in marker_genes_batches:
    # Construct the prompt with marker genes from the current batch
    prompt = prompt_base + "\n" + "\n".join([cluster_info['marker'] for cluster_info in batch])
    
    data = {
        "messages": [
            {"role": "system", "content": prompt}
        ]
    }
    
    # Construct the curl command
    curl_command = f'curl https://{RESOURCE_NAME}.openai.azure.com/openai/deployments/{DEPLOYMENT_NAME}/chat/completions?api-version={API_VERSION} ' \
                   f'-H "Content-Type: application/json" ' \
                   f'-H "api-key: {API_KEY}" ' \
                   f'-d \'{json.dumps(data)}\''
    
    # Execute the curl command
    process = subprocess.Popen(curl_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    if process.returncode == 0:
        try:
            # Attempt to parse the response
            response_data = json.loads(stdout)
            if 'choices' in response_data:
                predicted_types = response_data['choices'][0]['message']['content'].split('\n')
                
                # Update the predicted_cell_types dictionary with the new predictions
                for cluster_info, predicted_type in zip(batch, predicted_types):
                    cluster_name = list(level2_data.keys())[list(level2_data.values()).index(cluster_info)]
                    predicted_cell_types[cluster_name] = predicted_type.strip()
            else:
                print("The 'choices' key is not present in the response. Please check the API response structure.")
        except json.JSONDecodeError:
            print("Failed to decode JSON from the subprocess output. Please check the output format.")
    else:
        print(f"Error executing curl command: {stderr.decode('utf-8')}")



### Step 4: Level 3 Cell Type Annotation
import json
import os
import subprocess
import matplotlib.pyplot as plt

# Assuming the structured_data and predicted_cell_types dictionaries from previous steps are available

# Step 4: Level 3 Cell Type Annotation

# Initialize an empty dictionary to store the level 3 cluster name and its predicted cell type
predicted_cell_types_level3 = {}

# Prepare batches of marker genes for batch processing, including the parent cell type
marker_genes_batches_level3 = []
for cluster_name, cluster_info in structured_data['level3'].items():
    parent_cluster_name = cluster_info['parent']
    if parent_cluster_name in predicted_cell_types:
        parent_cell_type = predicted_cell_types[parent_cluster_name]
        marker_genes_batches_level3.append((cluster_name, f"{parent_cell_type} # {cluster_info['marker']}"))

# Set the prompt for GPT-4 in level 3
prompt_base_level3 = "Identify cell types of human lung cells using the following markers. Identify one cell type for each row. Only provide the cell type name. The format for each line is: parent cluster cell type # subcluster maker genes. Note that the predicted subcluster cell type must be more specific than the parent cluster cell type."

# Process each batch
batch_size = 10
for i in range(0, len(marker_genes_batches_level3), batch_size):
    batch = marker_genes_batches_level3[i:i+batch_size]
    prompt_level3 = prompt_base_level3 + "\n" + "\n".join([marker_genes_with_parent[1] for marker_genes_with_parent in batch])
    
    data = {
        "messages": [
            {"role": "system", "content": prompt_level3}
        ]
    }
    
    API_KEY = os.getenv("AZURE_API_KEY")  # Assuming the API key is stored in an environment variable
    RESOURCE_NAME = "DEPLOYMENT_NAME"
    DEPLOYMENT_NAME = "gpt-4-0125"
    API_VERSION = "2024-02-01"
    curl_command = f'curl https://{RESOURCE_NAME}.openai.azure.com/openai/deployments/{DEPLOYMENT_NAME}/chat/completions?api-version={API_VERSION} ' \
                   f'-H "Content-Type: application/json" ' \
                   f'-H "api-key: {API_KEY}" ' \
                   f'-d \'{json.dumps(data)}\''
    
    process = subprocess.Popen(curl_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    if process.returncode == 0:
        try:
            response_data = json.loads(stdout)
            if 'choices' in response_data:
                predicted_types = response_data['choices'][0]['message']['content'].split('\n')
                for cluster_info, predicted_type in zip(batch, predicted_types):
                    cluster_name = cluster_info[0]
                    predicted_cell_types_level3[cluster_name] = predicted_type.strip()
            else:
                print("The 'choices' key is not present in the response. Please check the API response structure.")
        except json.JSONDecodeError:
            print("Failed to decode JSON from the subprocess output. Please check the output format.")
    else:
        print(f"Error executing curl command: {stderr.decode('utf-8')}")

# Visualization of predicted cell types for level 3
plt.figure(figsize=(10, 6))
plt.bar(predicted_cell_types_level3.keys(), [len(cell_type) for cell_type in predicted_cell_types_level3.values()], color='skyblue')
plt.xlabel('Cluster Name')
plt.ylabel('Length of Predicted Cell Type Name')
plt.title('Predicted Cell Types for Level 3 Clusters')
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.show()


### Step 6: Level 4 Cell Type Annotation
import json
import os
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Assuming the structured_data and predicted_cell_types dictionaries from previous steps are available

# Step 6: Level 4 Cell Type Annotation

# 6.1 Prepare Level 4 Data
# Assuming the structured_data dictionary has been loaded with level4 data as per the provided data structure

# 6.2 Retrieve Updated Dictionary
# Assuming predicted_cell_types_level3 dictionary is updated with level3 annotations

# 6.3 Annotate Level 4 Data
predicted_cell_types_level4 = {}
marker_genes_batches_level4 = []
for cluster_name, cluster_info in structured_data['level4'].items():
    parent_cluster_name = cluster_info['parent']
    if parent_cluster_name in predicted_cell_types_level3:
        parent_cell_type = predicted_cell_types_level3[parent_cluster_name]
        marker_genes_batches_level4.append((cluster_name, f"{parent_cell_type} # {cluster_info['marker']}"))

# Set the prompt for GPT-4 in level 4
prompt_base_level4 = "Identify cell types of human lung cells using the following markers. Identify one cell type for each row. Only provide the cell type name. Some could be a mixture of multiple cell types. Some could be unknown cell types. The format for each line is: parent cluster cell type # subcluster maker genes. Note that the predicted subcluster cell type must be more specific than the parent cluster cell type."

# Process each batch
batch_size = 10
for i in range(0, len(marker_genes_batches_level4), batch_size):
    batch = marker_genes_batches_level4[i:i+batch_size]
    prompt_level4 = prompt_base_level4 + "\n" + "\n".join([marker_genes_with_parent[1] for marker_genes_with_parent in batch])
    
    data = {
        "messages": [
            {"role": "system", "content": prompt_level4}
        ]
    }
    
    API_KEY = os.getenv("AZURE_API_KEY")  # Assuming the API key is stored in an environment variable
    RESOURCE_NAME = "xxx"
    DEPLOYMENT_NAME = "gpt-4-0125"
    API_VERSION = "2024-02-01"
    curl_command = f'curl https://{RESOURCE_NAME}.openai.azure.com/openai/deployments/{DEPLOYMENT_NAME}/chat/completions?api-version={API_VERSION} ' \
                   f'-H "Content-Type: application/json" ' \
                   f'-H "api-key: {API_KEY}" ' \
                   f'-d \'{json.dumps(data)}\''
    
    process = subprocess.Popen(curl_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    if process.returncode == 0:
        try:
            response_data = json.loads(stdout)
            if 'choices' in response_data:
                predicted_types = response_data['choices'][0]['message']['content'].split('\n')
                for cluster_info, predicted_type in zip(batch, predicted_types):
                    cluster_name = cluster_info[0]
                    predicted_cell_types_level4[cluster_name] = predicted_type.strip()
            else:
                print("The 'choices' key is not present in the response. Please check the API response structure.")
        except json.JSONDecodeError:
            print("Failed to decode JSON from the subprocess output. Please check the output format.")
    else:
        print(f"Error executing curl command: {stderr.decode('utf-8')}")

# 6.4 Update Dictionary with Level 4 Annotations
# The dictionary is updated in the loop above

# 6.5 Verification and Quality Control
# Manual verification and quality control steps would be performed here

# 6.6 Finalize and Save Updates
# Assuming the dictionary is stored in a JSON file for future use
updated_dictionary_path = '/path/to/updated_dictionary.json'
if not os.path.exists(updated_dictionary_path):
    with open(updated_dictionary_path, 'w') as file:
        json.dump(predicted_cell_types_level4, file)  # Save the level 4 updates

# Note: Replace '/path/to/updated_dictionary.json' with the actual path where you want to save the updated dictionary

### Step 8: Final Output Preparation
import json
import os
from datetime import datetime

def compile_annotations(structured_data, predicted_cell_types_level2, predicted_cell_types_level3, predicted_cell_types_level4):
    """
    Compiles annotations from levels 2, 3, and 4 into a comprehensive dictionary.
    
    Parameters:
    - structured_data: The structured data loaded from the JSON file.
    - predicted_cell_types_level2: Dictionary containing predicted cell types for level 2 clusters.
    - predicted_cell_types_level3: Dictionary containing predicted cell types for level 3 clusters.
    - predicted_cell_types_level4: Dictionary containing predicted cell types for level 4 clusters.
    
    Returns:
    A dictionary with levels as keys and dictionaries of cluster names and their predicted cell types as values.
    """
    compiled_annotations = {
        "Level 2": predicted_cell_types,
        "Level 3": predicted_cell_types_level3,
        "Level 4": predicted_cell_types_level4
    }
    return compiled_annotations

def save_to_json(dictionary, filename):
    """
    Saves a dictionary to a JSON file with pretty formatting.
    
    Parameters:
    - dictionary: The dictionary to save.
    - filename: The name of the file to save the dictionary to.
    """
    with open(filename, 'w') as file:
        json.dump(dictionary, file, indent=4)

# Assuming structured_data, predicted_cell_types_level2, predicted_cell_types_level3, and predicted_cell_types_level4 are available
# For demonstration, these would be loaded or defined here based on previous steps

# Compile the updated dictionary
compiled_annotations = compile_annotations(structured_data, predicted_cell_types, predicted_cell_types_level3, predicted_cell_types_level4)

# Include metadata
metadata = {
    "date_of_analysis": datetime.now().strftime("%Y-%m-%d"),
    "GPT-4_API_version": "4.x.x"
}

final_output = {
    "metadata": metadata,
    "annotations": compiled_annotations
}

# Save the compiled dictionary to a JSON file
output_filename = "final_cell_type_annotations.json"
save_to_json(final_output, output_filename)

# Final checks and execution would be performed manually or with additional scripts as needed