from typing import Dict, Optional, Union, List

from autogen import Agent, AssistantAgent, UserProxyAgent, MyAssistantAgent
import chainlit as cl
import os
import autogen
from autogen.coding.jupyter import LocalJupyterServer
from autogen.coding.jupyter import JupyterCodeExecutor, JupyterConnectionInfo
from chainlitagent import ChainlitUserProxyAgent_decision_1, ChainlitUserProxyAgent_decision_2
from autogen.agentchat.contrib.retrieve_user_proxy_agent import RetrieveUserProxyAgent
from config import config_list
from typing_extensions import Annotated


llm_config = {"config_list": config_list, "seed": 42, "timeout": 1000, "temperature": 0}


decisioner1 = ChainlitUserProxyAgent_decision_1(
        name="decision_proxy",
        human_input_mode="ALWAYS",
        llm_config=llm_config,
        #max_consecutive_auto_reply=10,
        #is_termination_msg=lambda x: x.get("content", "").rstrip().endswith("TERMINATE"),
        code_execution_config=False,
        system_message="you are a proxy.",
    )

decisioner2 = ChainlitUserProxyAgent_decision_2(
        name="decision_proxy",
        human_input_mode="ALWAYS",
        llm_config=llm_config,
        #max_consecutive_auto_reply=10,
        #is_termination_msg=lambda x: x.get("content", "").rstrip().endswith("TERMINATE"),
        code_execution_config=False,
        system_message="you are a proxy.",
    )
