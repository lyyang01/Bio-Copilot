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


llm_config = {"config_list": config_list, "seed": 42, "timeout": 1000, "temperature": 0}


planner_message = f"""
Your are responsible for providing a description of the input dataset. You can first write code to print out the required dataset information, let data_proxy execute it, and then describe the dataset based on the returned results.
When doing the second step, i.e., describing the dataset, do not include code blocks.
"""

'''
planner_message = f"""
Your are responsible for providing a description of the input dataset. You can first write code to print out the required dataset information, let data_proxy execute it, and describe the dataset based on the returned results.
"""
'''
#if you think you do not have enough information, you can ask 'user_proxy' again.
data_agent = ChainlitAssistantAgent(
        name="data_agent",
        llm_config=llm_config,
        system_message=planner_message,
)


user_message = f"""
You are responsible for code execution from 'data_agent' and return the results to 'data_agent'.
"""
data_proxy = ChainlitUserProxyAgent(
        name="data_proxy",
        human_input_mode="NEVER",
        llm_config=llm_config,
        #max_consecutive_auto_reply=10,
        #is_termination_msg=lambda x: x.get("content", "").rstrip().endswith("TERMINATE"),
        code_execution_config={
            "executor": "ipython-embedded",
            "ipython-embedded": {"output_dir": "./coding","timeout": 10000},
        },
        system_message=user_message
    )


def custom_speaker_selection_func_data(last_speaker: Agent, groupchat: autogen.GroupChat):
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
            return data_proxy
        if last_speaker is data_proxy:
            return data_agent

        