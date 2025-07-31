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




jupytertoken = "8ec3a18eab7aa7d0deb068056ac4b8db4589c960549c9611"#"7fec73b094676d16ae4a85928f176a3bcf7aa7380e12f844"


llm_config = {"config_list": config_list, "seed": 42, "timeout": 1000, "temperature": 0}


# ===== Begin Skip =====
#skip_message = f"""
#You are a biometric analysis expert and receive a rough plan for the current step from 'skip_proxy'. You need to recommend whether to skip the current step directly or whether to skip the subsequent detailed planning execution.
#"""
#Note that your advantage is to make good use of loops in planning. Therefore, improvements can be made to the steps of planning using cycles in the plan. So improvements can be made to those steps in the plan where loops can be used.

#if you think you do not have enough information, you can ask 'user_proxy' again.
#skip = ChainlitAssistantAgent(
#        name="Skip",
#        llm_config=llm_config,
#        system_message=skip_message
#    )

skip_message = f"""
You provide an interface for human to decide whether to skip the current step or refined planner.
"""
skip_proxy = ChainlitUserProxyAgent_skip(
        name="skip_proxy",
        human_input_mode="ALWAYS",
        llm_config=llm_config,
        #max_consecutive_auto_reply=10,
        #is_termination_msg=lambda x: x.get("content", "").rstrip().endswith("TERMINATE"),
        code_execution_config=False,
        system_message=skip_message,
    )
# ===== End Skip =====


# ===== Begin Plan Refiner =====
#contains plan refiner, refiner reviewer, user_proxy_refine 
plan_refiner = ChainlitAssistantAgent(
        name="plan_refiner",
        llm_config=llm_config,
        # the default system message of the AssistantAgent is overwritten here
        system_message='''You are a plan refiner, responsible for incorporating more details into each step of the plan you receive.
        Assume all needed packages are installed and should not write install software as a separate step.
        '''
        #You are a plan refiner responsible for breaking down a complex task into steps that can be implemented through code.
        #Do not suggest concrete code. For any action beyond writing code or reasoning, convert it to a step that can be implemented by writing code. 
        #For example, browsing the web can be implemented by writing code that reads and prints the content of a web page.
        #Try to include recommended packages in your refinement step plan.
        #Note that you do not need to include environment setup in your refined plan.
        #You can include some necessary code blocks.
    
)


refiner_reviewer_msg = f"""
You are a reviewer with bioinformatics expertise and responsible for reviewing the plan from 'plan_refiner'. Your anwser only include suggestions for revisions which can be in the form of a paragraph or key points. Do not include in your answer any revised plans based on the suggestions.
Specifically, you are responsible for one thing: correct the wrong part, reduce unnecessary operations, and supplement the missing parts in the refined plan.
If no suggestions on the plan, just return your answer as "Great job!".
"""

plan_refiner_reviewer = ChainlitAssistantAgent(
        name="Re-Reviewer",
        llm_config=llm_config,
        # the default system message of the AssistantAgent is overwritten here
        system_message=refiner_reviewer_msg

#'''
#        You are a reviewer with bioinformatics expertise and need to review refined plans from 'plan_refiner'. Then you have to return your modification comments to the #'plan_refiner'.
#If you feel there is nothing that needs to be modified, you need to return your answer as "Great job!".
#What you need to do is to supplement the missing parts and reduce unnecessary operations. For example, unnecessary packages are imported or some functions are missing necessary #parameter settings.
#        '''
        #You are a plan refiner responsible for breaking down a complex task into steps that can be implemented through code.
        #Do not suggest concrete code. For any action beyond writing code or reasoning, convert it to a step that can be implemented by writing code. 
        #For example, browsing the web can be implemented by writing code that reads and prints the content of a web page.
        #Try to include recommended packages in your refinement step plan.
        #Note that you do not need to include environment setup in your refined plan.
    
)

user_proxy_refine = ChainlitUserProxyAgent_replan(
            name="planrefine_proxy",
            human_input_mode="ALWAYS",
            system_message="You are the terminator. Terminate the conversation promptly when there are no further modifications needed for planning.",
            max_consecutive_auto_reply=0,
            code_execution_config=False,
            llm_config=llm_config
)
# ===== End Plan Refiner =====

# ===== Begin coder =====
#contain coder, code_proxy_hard, code_proxy_hard2, checker
coder_message = f"""
    You are a coder with bioinformatical expertise and need to write code based on the plan description from 'code_proxy'.
    Then, you need to make changes to the code refer to the instructions from 'Checker' and your own knowledge.
    """
    
coder = ChainlitAssistantAgent(
            name="Coder",
            llm_config=llm_config,
            system_message=coder_message
        )

    #user proxy for coder
user_message3 = f"""
        You mainly return the execution result of the code to 'Checker' and 'Coder'.
        """
server = LocalJupyterServer()
#server = JupyterConnectionInfo(host="127.0.0.1", use_https=False, port=8888, token=jupytertoken)

code_proxy_hard = ChainlitUserProxyAgent(
        name="code_proxy",
        human_input_mode="TERMINATE",
        max_consecutive_auto_reply=5,
        is_termination_msg=lambda x: x.get("content", "").rstrip().endswith("TERMINATE"),
        #code_execution_config={
        #    "executor": "jupyter",
        #    "jupyter":{
        #        "output_dir": "/data/yangliu/llm-agent/coding",
        #        "jupyter_server": server,
        #    }
        #},
        code_execution_config={
            "executor": "ipython-embedded",
            "ipython-embedded": {"output_dir": "./coding","timeout": 10000},
        },
        system_message=user_message3
    )

code_proxy_hard2 = ChainlitUserProxyAgent_new(
        name="code_proxy2",
        human_input_mode="ALWAYS",
        max_consecutive_auto_reply=10,
        is_termination_msg=lambda x: x.get("content", "").rstrip().endswith("TERMINATE"),
        code_execution_config=False,
        system_message='''You provide a human-agent interaction interface to modify the code generated by 'Coder'.'''
    )

code_proxy_simple = ChainlitUserProxyAgent(
        name="code_proxy",
        human_input_mode="TERMINATE",
        max_consecutive_auto_reply=10,
        is_termination_msg=lambda x: x.get("content", "").rstrip().endswith("TERMINATE"),
        code_execution_config={
            "work_dir": "./coding",
            "use_docker": False,
        },
        system_message=user_message3
    )
    #import pdb
    #pdb.set_trace()

checker_message = f"""
    You are a code checker that receives code feedback from 'code_proxy'. 
    If the code feedback shows an error, you need to give 'Coder' suggestions on how to modify the code to correct the error.
    Do not suggest concrete code. For any action beyond writing code or reasoning, convert it to a step that can be implemented by writing code. For example, browsing the web can be implemented by writing code that reads and prints the content of a web page.
    when you receive exitcode: 0 from 'chat_manager', you should end your answer with 'Great job!'.
    Never end with "TERMINATE"!
    """
    #You do not need to return the full code to 'Coder', just where the changes were.
    #Never add "TERMINATE" at the end of your answer!.

checker = ChainlitAssistantAgent(
            name="Checker",
            llm_config=llm_config,
            system_message=checker_message
        )
    # ===== End coder =====

# == Begin Retriever Agent ==
retriever_call = ChainlitAssistantAgent(
    name="retriever_call", 
    llm_config={"config_list": config_list, "seed": 42, "timeout": 1000, "temperature": 0},
    system_message="call the retriever",
)

retriever_executor = ChainlitAssistantAgent_ex(
        name="retriever_executor", 
        llm_config={"config_list": config_list, "seed": 42, "timeout": 1000, "temperature": 0},
        system_message="execute the retriever",
    )

retriever_code = RetrieveUserProxyAgent(
name="retriever",
human_input_mode="NEVER",
max_consecutive_auto_reply=5,
retrieve_config={
    "task": "code",
    "docs_path": [
        "./tools/scanvi.md",
        "./tools/decoupler.md",
    ],
},
code_execution_config=False,  # set to False if you don't want to execute the code
description="Assistant who has extra content retrieval power for solving difficult problems in coding.",
llm_config={"config_list": config_list, "seed": 42, "timeout": 1000, "temperature": 0},
)


def retrieve_content(
    message: Annotated[
        str,
        "Refined message which keeps the original meaning of the current step and can be used to retrieve content for code generation.",
        #"Extract specific task from the top line 'Refined Step ...' or 'Step ...' or 'Plan for ...' that can be used to retrieve content for code generation.",
    ],
    #n_results: Annotated[int, "number of results"] = 3,
) -> str:
    #data_path = "./tools/cell2location.md"
    data_path = "./tools/multivi.md"
    with open(data_path, "r") as f:
        databse = f.readlines()
    PROMPT_CODE = f"""
You're a retrieve augmented coding assistant. You wrtie code for user's requirements or make modifications of codes based on your own knowledge and the context provided by the user. It should be noted that not all content in the Context must be included in your answer. If the current context does not contain the knowledge you need to answer, then ignore the context and answer using only your own knowledge. For code generation, you must obey the following rules:
Rule 1. You MUST NOT install any packages because all the packages needed are already installed.
Rule 2. You must follow the formats below to write your code:
```python
# your code
```

User's requirement is: {message}

Context is: {databse}
"""
    return PROMPT_CODE

for caller in [coder]:
    #import pdb
    #pdb.set_trace()
    d_retrieve_content = caller.register_for_llm(
        description="retrieve content for code generation and question answering.", api_style="function"
        #description="help coder to write code on cell annotation with decoupler, scanvi integration with scarches.", api_style="function"
        #description="help coder to write code. you need always call the function named retrieve_content.", api_style="function"
    )(retrieve_content)

for executor in [retriever_executor]:
    executor.register_for_execution()(d_retrieve_content)




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

def custom_speaker_selection_func_replan(last_speaker: Agent, groupchat: autogen.GroupChat):
    """Define a customized speaker selection function.
    A recommended way is to define a transition for each speaker in the groupchat.

    Returns:
        Return an `Agent` class or a string from ['auto', 'manual', 'random', 'round_robin'] to select a default method to use.
    """
    messages = groupchat.messages

    if len(messages) <= 1:
        return plan_refiner
    
    if last_speaker is user_proxy_refine:
        return plan_refiner

    if last_speaker is plan_refiner:
        if len(messages) > 2:
            return user_proxy_refine
        else:
            return plan_refiner_reviewer
    
    if last_speaker is plan_refiner_reviewer:
        #import pdb
        #pdb.set_trace()
        if 'Great job' in messages[-1]['content']:
            return user_proxy_refine
        #    return plan_proxy
        #else:
        return plan_refiner

def custom_speaker_selection_func_code(last_speaker: Agent, groupchat: autogen.GroupChat):
    """Define a customized speaker selection function.
    A recommended way is to define a transition for each speaker in the groupchat.

    Returns:
        Return an `Agent` class or a string from ['auto', 'manual', 'random', 'round_robin'] to select a default method to use.
    """
    messages = groupchat.messages

    if len(messages) <= 1:
        return coder
        #return retriever_call
    #TODO
    #if last_speaker is code_proxy_hard2:
    #    return coder
    if last_speaker is code_proxy_hard2:
        #import pdb
        #pdb.set_trace()
        #if len(code_proxy_hard2._human_input) > 0:
        #    if code_proxy_hard2._human_input[-1] == "retrieval":
        #        #import pdb
        #        #pdb.set_trace()
        #        return retriever_call
        #    else:
        #        return coder
        #else:
        #    return coder
        return coder

    #if last_speaker is retriever_call:
    #    return retriever_executor

    if last_speaker is retriever_executor:
        #import pdb
        #pdb.set_trace()
        return coder

    if last_speaker is coder:
        #import pdb
        #pdb.set_trace()
        if messages[-1]["content"] == '':
            return retriever_executor

        if "psudo-exit" in code_proxy_hard2._human_input:
            return code_proxy_hard
        #elif len(code_proxy_hard2._human_input) > 0:
        #    if code_proxy_hard2._human_input[-1] == "retrieval":
        #        return retriever_call
        else:
            return code_proxy_hard2

    elif last_speaker is code_proxy_hard:
        # Always let the user to speak after the planner
        if "USER-MARK" in messages[-2]["content"]:
            return coder
        else:
            return checker

    elif last_speaker is checker:
        ###
        #if "USER-MARK" in messages[-1]["content"]:
        #    return coder
        if 'TERMINATE' in messages[-1]["content"]:
            messages[-1]["content"] = messages[-1]["content"].replace("TERMINATE", "USER-MARK") + "\n" + "TERMINATE"
            return code_proxy_hard
        else:
            return coder
    else:
        return "random"


'''
def custom_speaker_selection_func_code(last_speaker: Agent, groupchat: autogen.GroupChat):
    """Define a customized speaker selection function.
    A recommended way is to define a transition for each speaker in the groupchat.

    Returns:
        Return an `Agent` class or a string from ['auto', 'manual', 'random', 'round_robin'] to select a default method to use.
    """
    messages = groupchat.messages

    if len(messages) <= 1:
        return coder
    #TODO
    #if last_speaker is code_proxy_hard2:
    #    return coder
    if last_speaker is code_proxy_hard2:
        #import pdb
        #pdb.set_trace()
        if len(code_proxy_hard2._human_input) > 0:
            if code_proxy_hard2._human_input[-1] == "retrieval":
                #import pdb
                #pdb.set_trace()
                return retriever_call
            else:
                return coder
        else:
            return coder
    
    #if last_speaker is retriever_call:
    #    return retriever_executor

    if last_speaker is retriever_executor:
        #import pdb
        #pdb.set_trace()
        return coder

    if last_speaker is coder:
        #import pdb
        #pdb.set_trace()
        if "psudo-exit" in code_proxy_hard2._human_input:
            return code_proxy_hard
        #elif len(code_proxy_hard2._human_input) > 0:
        #    if code_proxy_hard2._human_input[-1] == "retrieval":
        #        return retriever_call
        else:
            return code_proxy_hard2

    elif last_speaker is code_proxy_hard:
        # Always let the user to speak after the planner
        if "USER-MARK" in messages[-2]["content"]:
            return coder
        else:
            return checker

    elif last_speaker is checker:
        ###
        #if "USER-MARK" in messages[-1]["content"]:
        #    return coder
        if 'TERMINATE' in messages[-1]["content"]:
            messages[-1]["content"] = messages[-1]["content"].replace("TERMINATE", "USER-MARK") + "\n" + "TERMINATE"
            return code_proxy_hard
        else:
            return coder
    else:
        return "random"
'''