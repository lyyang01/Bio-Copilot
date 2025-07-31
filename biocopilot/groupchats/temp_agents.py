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
You need to determine whether the current task is suitable for three-stage planning. Three-stage planning means planning the main body of the task into three stages, including data preprocessing, method execution, and result saving. Please give your judgment and reasons.
"""
#if you think you do not have enough information, you can ask 'user_proxy' again.
temper = ChainlitAssistantAgent(
        name="decision_maker",
        llm_config=llm_config,
        system_message=planner_message,
        #code_execution_config={
        #    "executor": "ipython-embedded",
        #    "ipython-embedded": {"output_dir": "./coding","timeout": 10000},
        #},
    )

user_message = f"""
you are a proxy.
"""
temp_proxy = ChainlitUserProxyAgent_plan(
        name="decision_proxy",
        human_input_mode="ALWAYS",
        llm_config=llm_config,
        #max_consecutive_auto_reply=10,
        #is_termination_msg=lambda x: x.get("content", "").rstrip().endswith("TERMINATE"),
        code_execution_config=False,
        system_message=user_message
    )
