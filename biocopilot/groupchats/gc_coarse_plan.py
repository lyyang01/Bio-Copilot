from typing import Dict, Optional, Union, List

from autogen import Agent, AssistantAgent, UserProxyAgent, MyAssistantAgent
import chainlit as cl
import os
import autogen
from autogen.coding.jupyter import LocalJupyterServer
from autogen.coding.jupyter import JupyterCodeExecutor, JupyterConnectionInfo
from chainlitagent import ChainlitAssistantAgent, ChainlitCodeAssistantAgent, ChainlitUserProxyAgent, ChainlitUserProxyAgent_new, ChainlitUserProxyAgent_skip, ChainlitAssistantAgent_ex, ChainlitUserProxyAgent_plan, ChainlitUserProxyAgent_replan
from autogen.agentchat.contrib.retrieve_user_proxy_agent import RetrieveUserProxyAgent
from config import config_list
from typing_extensions import Annotated

#TODO
from .gc_read_data import data_agent, data_proxy
#set openai api parameters


#config
llm_config = {"config_list": config_list, "seed": 42, "timeout": 1000, "temperature": 0}

#build planner
planner_message = f"""
You are a planner with bioinformatics expertise and need to plan the bioinformatics task based on all messages.
"""
#if you think you do not have enough information, you can ask 'user_proxy' again.
planner = ChainlitAssistantAgent(
        name="Planner",
        llm_config=llm_config,
        system_message=planner_message,
        #code_execution_config={
        #    "executor": "ipython-embedded",
        #    "ipython-embedded": {"output_dir": "./coding","timeout": 10000},
        #},
    )

#bulid reviewer
reviewer_message = f"""
You are a reviewer with bioinformatics expertise and responsible for reviewing the plan from 'Planner'. Your anwser only include suggestions for revisions which can be in the form of a paragraph or key points. Do not include in your answer any revised plans based on the suggestions.
Specifically, you are responsible for two things:
(1) make suggestions on the structure or logic of the plan that are unreasonable. For example, some steps can be combined together, some steps need to be broken down, some steps need to be removed, or some steps need to be added.
(2) if 'Planner' proposes multiple methods in a certain step, you need to choose one of them based on your own knowledge and give the reason.
If no suggestions on the plan, just return your answer as "Great job!".
"""

"""
You are a reviewer with bioinformatics expertise and responsible for reviewing the plan from 'Planner'. Your anwser only include suggestions for revisions which can be in the form of a paragraph or key points. Do not include in your answer any revised plans based on the suggestions.
Specifically, you are responsible for two things:
(1) make suggestions on the structure or logic of the plan that are unreasonable. For example, some steps can be combined together, some steps need to be broken down, some steps need to be removed, or some steps need to be added.
(2) if 'Planner' proposes multiple methods in a certain step, you need to choose one of them based on your own knowledge and give the reason.
If no suggestions on the plan, just return your answer as "Great job!".
"""


"""
You are a reviewer with bioinformatics expertise and responsible for reviewing the plan from 'Planner'. Your anwser only include suggestions for revisions which can be in the form of a paragraph or key points. Do not include any revised plans in your answer!
Specifically, you are responsible for two things:
(1) make suggestions on combining steps that are closely related and easy to implement together.
(2) if 'Planner' proposes multiple methods in a certain step, you need to choose one of them based on your own knowledge and give the reason.
If no suggestions on the plan, just return your answer as "Great job!".
"""

"""
You are a reviewer with bioinformatics expertise and responsible for reviewing the plan from 'Planner'. Your anwser only include suggestions for revisions which can be in the form of a paragraph or key points. Do not include in your answer any revised plans based on the suggestions.
Specifically, you are responsible for two things:
(1) make suggestions on the structure or logic of the plan that are unreasonable. For example, some steps can be combined together, some steps need to be broken down, some steps need to be removed, or some steps need to be added.
(2) if 'Planner' proposes multiple methods in a certain step, you need to choose a suitable method that may not be among the methods mentioned in 'Planner' based on your own knowledge and data description and then give the reason.
If no suggestions on the plan, just return your answer as "Great job!".
"""
#Note that your advantage is to make good use of loops in planning. Therefore, improvements can be made to the steps of planning using cycles in the plan. So improvements can be made to those steps in the plan where loops can be used.

#if you think you do not have enough information, you can ask 'user_proxy' again.
reviewer = ChainlitAssistantAgent(
        name="Reviewer",
        llm_config=llm_config,
        system_message=reviewer_message
    )
    
#build planner proxy
user_message = f"""
You are a bioinformatics researcher who can help 'Planner' to finish the bioinformatics analysis process.
"""
plan_proxy = ChainlitUserProxyAgent_plan(
        name="plan_proxy",
        human_input_mode="ALWAYS",
        llm_config=llm_config,
        #max_consecutive_auto_reply=10,
        #is_termination_msg=lambda x: x.get("content", "").rstrip().endswith("TERMINATE"),
        code_execution_config=False,
        system_message=user_message
    )


#bulid reviewer
recommend_message = f"""
You are a recommender with bioinformatics expertise and responsible for recommending one best method for user task based on retrieved content.
If the retrieved content does not contain any useful information, please strictly reply with "No recommended methods. There is no relevant information for this task in the current knowledge base. Please update the knowledge base or use your own knowledge to make plan or write code."
Do not contain specific code blocks in your recommendations.
"""
#Note that your advantage is to make good use of loops in planning. Therefore, improvements can be made to the steps of planning using cycles in the plan. So improvements can be made to those steps in the plan where loops can be used.

#if you think you do not have enough information, you can ask 'user_proxy' again.
recommender = ChainlitAssistantAgent(
        name="Recommender",
        llm_config=llm_config,
        system_message=recommend_message
    )



recommend_proxy = ChainlitUserProxyAgent(
        name="recommend_proxy",
        human_input_mode="NEVER",
        llm_config=llm_config,
        code_execution_config=False,
        #max_consecutive_auto_reply=10,
        #is_termination_msg=lambda x: x.get("content", "").rstrip().endswith("TERMINATE"),
        system_message=""
    )

# == Begin Retriever Agent ==
retriever_call_plan = ChainlitAssistantAgent(
    name="retriever_call", 
    llm_config={"config_list": config_list, "seed": 42, "timeout": 1000, "temperature": 0},
    system_message="call the retriever",
)

retriever_executor_plan = ChainlitAssistantAgent_ex(
        name="retriever_executor", 
        llm_config={"config_list": config_list, "seed": 42, "timeout": 1000, "temperature": 0},
        system_message="execute the retriever",
    )

retriever_plan = RetrieveUserProxyAgent(
name="retriever",
human_input_mode="NEVER",
max_consecutive_auto_reply=5,
retrieve_config={
    "task": "code",
    "docs_path": [
        "./kb/data-integration.txt",
        "./kb/cell-annotation.txt",
    ],
},
code_execution_config=False,  # set to False if you don't want to execute the code
description="Assistant who has extra content retrieval power for solving difficult problems in coding.",
llm_config={"config_list": config_list, "seed": 42, "timeout": 1000, "temperature": 0},
)

def retrieve_content(
    message: Annotated[
        str,
        #"Refined message which keeps the original meaning of the current step and can be used to retrieve content for code generation.",
        "Extract the bioinformatics task from the user request.",
    ],
) -> str:
    data_path = "./kb/database.txt"
    with open(data_path, "r") as f:
        databse = f.readlines()
    PROMPT_CODE = f"""
You are a bioinformatics analysis expert. Your task is to recommend the most suitable method for a given bioinformatics analysis task based on a existing methods comparison information. Your recommendations should be accurate, practical, and tailored to the specific requirements of the task. Based on the provided existing methods comparison, analyze the input task description and recommend the best method(s) for the bioinformatics analysis. Your recommendation should include:
(1) The most suitable method(s) for the task.
(2) A brief explanation of why the method(s) are appropriate.

Bioinformatics Task is: 
{message}

Existing Methods Comparison is: 
{databse}

If the Existing Methods Comparison does not contain any useful information, please strictly reply with "No recommended methods. There is no relevant information for this task in the current knowledge base. Please update the knowledge base or use your own knowledge to make plan or write code."
Do not contain specific code blocks in your recommendations.
"""
    return PROMPT_CODE

for caller in [retriever_call_plan]:
    #import pdb
    #pdb.set_trace()
    d_retrieve_content = caller.register_for_llm(
        #description="help planner to make a plan on method selection.", api_style="function"
        description="make a method recommendation on user task.", api_style="function"
        #description="help coder to write code on cell annotation with decoupler, scanvi integration with scarches.", api_style="function"
        #description="help coder to write code. you need always call the function named retrieve_content.", api_style="function"
    )(retrieve_content)

for executor in [retriever_executor_plan]:
    executor.register_for_execution()(d_retrieve_content)

'''
def custom_speaker_selection_func_plan(last_speaker: Agent, groupchat: autogen.GroupChat):
    """Define a customized speaker selection function.
    A recommended way is to define a transition for each speaker in the groupchat.

    Returns:
        Return an `Agent` class or a string from ['auto', 'manual', 'random', 'round_robin'] to select a default method to use.
    """
    messages = groupchat.messages

    if len(messages) <= 1:
        #return planner
        return retriever_call_plan
    
    if last_speaker is retriever_executor_plan:
        return planner

    if last_speaker is retriever_call_plan:
        return retriever_executor_plan
        
    if last_speaker is plan_proxy:
        return planner

    if last_speaker is planner:
        #if messages[-1]["content"] == '':
        #    return retriever_executor_plan

        if len(messages) > 2+2:
            return plan_proxy
        else:
            return reviewer
    
    if last_speaker is reviewer:
        #import pdb
        #pdb.set_trace()
        if 'Great job' in messages[-1]['content']:
            return plan_proxy
        #else:
        return planner
'''
#

def custom_speaker_selection_func_plan(last_speaker: Agent, groupchat: autogen.GroupChat):
    """Define a customized speaker selection function.
    A recommended way is to define a transition for each speaker in the groupchat.

    Returns:
        Return an `Agent` class or a string from ['auto', 'manual', 'random', 'round_robin'] to select a default method to use.
    """
    messages = groupchat.messages

    if len(messages) <= 1:
        return planner
    
    if last_speaker is plan_proxy:
        return planner

    if last_speaker is planner:
        if len(messages) > 2:
            return plan_proxy
        else:
            return reviewer
    
    if last_speaker is reviewer:
        #import pdb
        #pdb.set_trace()
        if 'Great job' in messages[-1]['content']:
            return plan_proxy
        #else:
        return planner

def custom_speaker_selection_func_retrieve(last_speaker: Agent, groupchat: autogen.GroupChat):
    """Define a customized speaker selection function.
    A recommended way is to define a transition for each speaker in the groupchat.

    Returns:
        Return an `Agent` class or a string from ['auto', 'manual', 'random', 'round_robin'] to select a default method to use.
    """
    messages = groupchat.messages

    if len(messages) <= 1:
        #return planner
        return retriever_call_plan

    if last_speaker is retriever_call_plan:
        return retriever_executor_plan
    
    if last_speaker is retriever_executor_plan:
        return recommender

'''
def custom_speaker_selection_func_plan(last_speaker: Agent, groupchat: autogen.GroupChat):
    """Define a customized speaker selection function.
    A recommended way is to define a transition for each speaker in the groupchat.

    Returns:
        Return an `Agent` class or a string from ['auto', 'manual', 'random', 'round_robin'] to select a default method to use.
    """
    messages = groupchat.messages

    if len(messages) <= 1:
        return data_agent
    if last_speaker is data_agent:
        if '```python' in messages[-1]['content']:
            return planner
        else:
            return planner
    if last_speaker is data_proxy:
        return data_agent


    if last_speaker is plan_proxy:
        return planner

    if last_speaker is planner:
        if 'exitcode:' in messages[-1]['content']:
            return data_agent
        if len(messages) > 2+4:
            return plan_proxy
        else:
            return reviewer
    
    if last_speaker is reviewer:
        #import pdb
        #pdb.set_trace()
        if 'Great job' in messages[-1]['content']:
            return plan_proxy
        #else:
        return planner
'''

'''
def custom_speaker_selection_func_plan(last_speaker: Agent, groupchat: autogen.GroupChat):
        """Define a customized speaker selection function.
        A recommended way is to define a transition for each speaker in the groupchat.

        Returns:
            Return an `Agent` class or a string from ['auto', 'manual', 'random', 'round_robin'] to select a default method to use.
        """
        messages = groupchat.messages

        if len(messages) <= 1:
            return data_agent
        #    return toolselect_caller
        if last_speaker is data_agent:
            if '```python' in messages[-1]['content']:
                return data_proxy
            else:
                return planner
        if last_speaker is data_proxy:
            return data_agent

        if last_speaker is plan_proxy:
            return planner

        if last_speaker is planner:
            if len(messages) > 8:
                return plan_proxy
            else:
                return toolselect_caller
        
        if last_speaker is toolselect_caller:
            return toolselect_executor

        if last_speaker is toolselect_executor:
            return summary_agent
        
        if last_speaker is summary_agent:
            return reviewer
        
        if last_speaker is reviewer:
            if 'Great job' in messages[-1]['content']:
                return plan_proxy
            return planner
'''