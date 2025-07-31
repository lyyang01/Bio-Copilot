
from typing import Dict, Optional, Union, List

from autogen import Agent, AssistantAgent, UserProxyAgent, MyAssistantAgent
import chainlit as cl
import os
import autogen
from autogen.coding.jupyter import LocalJupyterServer
from autogen.coding.jupyter import JupyterCodeExecutor, JupyterConnectionInfo
from autogen.agentchat import GroupChat, GroupChatManager
from get_prompt import get_prompt_plan, get_prompt_plan_nodata, get_prompt_coder_simple, get_prompt_coder_hard
from chainlitagent import ChainlitAssistantAgent, ChainlitCodeAssistantAgent, ChainlitUserProxyAgent, ChainlitUserProxyAgent_new, ChainlitUserProxyAgent_skip, ChainlitAssistantAgent_ex, ChainlitUserProxyAgent_plan, ChainlitUserProxyAgent_replan
import re
import math
from typing_extensions import Annotated

#gpt config
from config import config_list

#import selection function
#from agents import custom_speaker_selection_func_plan, custom_speaker_selection_func_replan, custom_speaker_selection_func_code
from agents import custom_speaker_selection_func_replan, custom_speaker_selection_func_code

#coarse-level-plan
from groupchats.gc_coarse_plan import planner, reviewer, plan_proxy, custom_speaker_selection_func_plan
from groupchats.gc_coarse_plan import retriever_executor_plan, retriever_call_plan, retriever_plan, recommender, recommend_proxy, custom_speaker_selection_func_retrieve

#decision groupchat
from groupchats.gc_decision_maker import decisioner1, decisioner2
#data groupchat
#from groupchats.gc_read_data import data_agent, data_proxy, custom_speaker_selection_func_data

#skip
from agents import skip_proxy

#fine-level-plan
from agents import plan_refiner, plan_refiner_reviewer, user_proxy_refine

#action-executor
from agents import coder, checker, code_proxy_hard, code_proxy_hard2, code_proxy_simple
from agents import retriever_call, retriever_executor, retriever_code


#
os.environ["CUDA_VISIBLE_DEVICES"]= '3'

async def ask_helper(func, **kwargs):
    res = await func(**kwargs).send()
    while not res:
        res = await func(**kwargs).send()
    return res



import math
@cl.on_chat_start
async def on_chat_start():
    #OPENAI_API_KEY = os.getenv('OPENAI_API_KEY')
    #try:
        #llm_config
    llm_config = {"config_list": config_list, "seed": 42, "timeout": 1000, "temperature": 0}

    ## Planner Groupchat
    #planner
    #reviewer
    #plan_proxy
    #config
    #build planner
    

    ## Skip
    #skip_proxy

    ## Refined Planner Groupchat
    #plan_refiner
    #plan_refiner_reviewer
    #user_proxy_refine

    ## Action Execution Groupchat
    #coder
    #checker
    #code_proxy_hard
    #code_proxy_hard2
    #
    #code_proxy_simple

    #========================================start chainlit server================================================================

    msg = cl.AskUserMessage(content=f"""Hello! I am a bioinformatics LLM-agent. What task would you like to get done? 
                        """, timeout=100000000, raise_on_timeout=False)
    message = await msg.send()
    
    
    TASK = message['output']
    
    temp = TASK.split('\n')

    goal = temp[0]
    data_list = []
    
    #if len(temp) > 1:
    for i in range(1, len(temp)):
        data_list.append(temp[i])
            #prompt = get_prompt_plan(data_list, goal)
    #else:
    #    prompt = get_prompt_plan_nodata(data_list, goal)



    await cl.Message(content=f"""Starting planning on task. If you think the plan is done, please to Continue!""").send()
    
    PROMPT_CODE = """Now, you're a retrieve augmented Recommender. You goal is to recommend the best method based on the user's task and context. Context generally contains the performance of different methods of a bioinformatics task in different dimensions. Thus, you need to extract useful information from the context.\nIf the context does not contain the knowledge you need, please reply with "No recommended methods. There is no relevant information for this task in the current knowledge base. Please update the knowledge base or use your own knowledge to make plan or write code.".\n\nUser's task is: {input_question}\n\nContext is:\n{input_context}"""
    #PROMPT_CODE = """Now, you're a retrieve augmented Planner. You goal is to plan the specific bioinformatics task with the best method based on the user's task and context. Context generally contains the performance of different methods of a bioinformatics task in different dimensions. Thus, you need to extract useful information from the context.\nIf the context does not contain the knowledge you need to plan the task, then ignore the context and make a plan using only your own knowledge.\n\nUser's task is: {input_question}\n\nContext is: {input_context}"""
    
    retriever_plan.customized_prompt ="user's task:\n" + goal + '\n\n' + PROMPT_CODE
    #read data groupchat
    MAX_ROUND = 4
    groupchat_retrieve = GroupChat(agents=[recommend_proxy, retriever_executor_plan, retriever_call_plan, recommender], messages=[], max_round=MAX_ROUND, speaker_selection_method=custom_speaker_selection_func_retrieve)
    manager_data = GroupChatManager(groupchat=groupchat_retrieve, llm_config=llm_config)
    if len(groupchat_retrieve.messages) == 0:
        chat_history = await cl.make_async(recommend_proxy.initiate_chat)(manager_data, message=TASK)
    elif len(groupchat_retrieve.messages) <  MAX_ROUND:
        chat_history = await cl.make_async(recommend_proxy.send)(manager_data, message=TASK,)
    elif len(groupchat_retrieve.messages) ==  MAX_ROUND:
        chat_history = await cl.make_async(recommend_proxy.send)(manager_data, message="exit",)

    recommend_description = groupchat_retrieve.messages[-1]['content']
    
    await cl.Message(content=f"""To plan or code?""").send()
    await cl.make_async(decisioner1.initiate_chat)(decisioner1, message="")
    #import pdb
    #pdb.set_trace()
    if decisioner1._human_input[0] == 'code':

        PROMPT_CODE = """Now, you're a retrieve augmented Coder. You goal is to write the code for the current step task based on your own knowledge and the context. Context generally contains how specific packages are used. Thus, you need to extract useful information from the context.\nNote that you do not need to follow the code in Context strictly. In other words, You can modify or discard the lines of code in the context that are not necessary for the task, or you can write new code yourself.\nIf Context does not contain the knowledge you need to write the code, then ignore Context and write the code using only your own knowledge.\n\nContext is: {input_context}
                """
        retriever_code.customized_prompt = TASK + '\n\n' + PROMPT_CODE
        code_prompt = TASK + '\n\n' + recommend_description
        # == Begin coder ==
        indicator = 0
        MAX_ROUND = 100
        groupchat = GroupChat(agents=[code_proxy_hard, code_proxy_hard2, coder, checker, retriever_call, retriever_executor], messages=[], max_round=MAX_ROUND, speaker_selection_method=custom_speaker_selection_func_code)
        manager = GroupChatManager(groupchat=groupchat, llm_config=llm_config)
        
        if len(groupchat.messages) == indicator:
            chat_history = await cl.make_async(code_proxy_hard2.initiate_chat)(manager, message=str(code_prompt)+ '\n\nretrieve content')
            #await user_proxy.initiate_chat( manager, message=message)
        elif len(groupchat.messages) <  MAX_ROUND + indicator:
            chat_history = await cl.make_async(code_proxy_hard.send)(manager, message=TASK,)
        elif len(groupchat.messages) ==  MAX_ROUND + indicator:  
            chat_history = await cl.make_async(code_proxy_hard.send)(manager, message="exit",)
        

    else:
        #== Begin planner groupchat ==
        #
        prompt = get_prompt_plan(TASK, recommend_method=recommend_description) 
        #prompt = TASK + "\n\n" + "I also provide you the data details as the following.\n" +agent_data_description
        #prompt = TASK

        
        groupchat_plan = GroupChat(agents=[plan_proxy, planner, reviewer], messages=[], max_round=100, speaker_selection_method=custom_speaker_selection_func_plan)
        manager_plan = GroupChatManager(groupchat=groupchat_plan, llm_config=llm_config)

        if len(groupchat_plan.messages) == 0:
            chat_history = await cl.make_async(plan_proxy.initiate_chat)(manager_plan, message=str(prompt))
                        #await user_proxy.initiate_chat( manager, message=message)
        elif len(groupchat_plan.messages) <  100:
            chat_history = await cl.make_async(plan_proxy.send)(manager_plan, message=TASK,)
        elif len(groupchat_plan.messages) ==  100:
            chat_history = await cl.make_async(plan_proxy.send)(manager_plan, message="exit",)
        
        #== End planner groupchat ==

        final_plan = groupchat_plan.messages[-1]['content']
        plan_str = final_plan

        await cl.Message(content=f"""To sub-steps or code?""").send()
        await cl.make_async(decisioner2.initiate_chat)(decisioner2, message="")
        #import pdb
        #pdb.set_trace()

        if decisioner2._human_input[0] == 'code':
            PROMPT_CODE = """Now, you're a retrieve augmented Coder. You goal is to write the code for the current step task based on your own knowledge and the context. Context generally contains how specific packages are used. Thus, you need to extract useful information from the context.\nNote that you do not need to follow the code in Context strictly. In other words, You can modify or discard the lines of code in the context that are not necessary for the task, or you can write new code yourself.\nIf Context does not contain the knowledge you need to write the code, then ignore Context and write the code using only your own knowledge.\n\nContext is: {input_context}
                """
            retriever_code.customized_prompt = TASK + '\n\n' + PROMPT_CODE
            code_prompt = plan_str
            # == Begin coder ==
            indicator = 0
            MAX_ROUND = 100
            groupchat = GroupChat(agents=[code_proxy_hard, code_proxy_hard2, coder, checker, retriever_call, retriever_executor], messages=[], max_round=MAX_ROUND, speaker_selection_method=custom_speaker_selection_func_code)
            manager = GroupChatManager(groupchat=groupchat, llm_config=llm_config)
            
            if len(groupchat.messages) == indicator:
                chat_history = await cl.make_async(code_proxy_hard2.initiate_chat)(manager, message=str(code_prompt)+ '\n\nretrieve content')
                #await user_proxy.initiate_chat( manager, message=message)
            elif len(groupchat.messages) <  MAX_ROUND + indicator:
                chat_history = await cl.make_async(code_proxy_hard.send)(manager, message=TASK,)
            elif len(groupchat.messages) ==  MAX_ROUND + indicator:  
                chat_history = await cl.make_async(code_proxy_hard.send)(manager, message="exit",)
        
        else:
            #import pdb
            #pdb.set_trace()
            sub_steps = plan_str.split('\n\n')#[1:-1]
            sub_steps_temp = []
            for i in sub_steps:
                pattern = r'Step \d'
                #if 'Step' in i:
                if re.search(pattern, i):
                    sub_steps_temp.append(i)
            sub_steps = sub_steps_temp

            
            indicator = 0
            previous_step = ""
            coder_messages = []
            tmp_coder_messages = ""
            for idx, sub_step in enumerate(sub_steps):
                ori_sub_step = sub_step
                
                # == Begin skip ==
                message2 = "Step {} is:\n".format(str(idx+1)) + sub_step + "\n\nYour goal: determine whether skip the current step or refined planning."
                await cl.Message(content=f"""Starting Skipping. Please decide whether to skip!""").send()
                await cl.make_async(skip_proxy.initiate_chat)(skip_proxy, message=message2)
                # == End skip ==
                
                indicator_ = skip_proxy._human_input

            
                if indicator_[0] == 'skip-step':
                    continue
                if indicator_[0] != 'skip-replan':
                    
                    # == Build Plan Refiner ==
                    message2 = "Step {} is:\n".format(str(idx+1)) + sub_step + "\n\nYour goal: refine the current step plan further."
                    
                    await cl.Message(content=f"""Starting refined planning on task. If you think the refined plan is done, please to Continue!""").send()
                    groupchat_replan = GroupChat(agents=[user_proxy_refine, plan_refiner, plan_refiner_reviewer], messages=[], max_round=100, speaker_selection_method=custom_speaker_selection_func_replan)
                    manager_replan = GroupChatManager(groupchat=groupchat_replan, llm_config=llm_config)
                    if len(groupchat_replan.messages) == 0:
                        chat_history = await cl.make_async(user_proxy_refine.initiate_chat)(manager_replan, message=message2)
                                    #await user_proxy.initiate_chat( manager, message=message)
                    elif len(groupchat_replan.messages) <  100:
                        chat_history = await cl.make_async(user_proxy_refine.send)(manager_replan, message=TASK,)
                    elif len(groupchat_replan.messages) ==  100:  
                        chat_history = await cl.make_async(user_proxy_refine.send)(manager_replan, message="exit",)

                    # == End Plan Refiner ==
                    for i in groupchat_replan.messages:
                        if 'Step' in i['content']:
                            sub_step = i['content']
                
                
                goals = "write the code for Step {}.".format(str(idx+1))
                #message2 ="the user task need: " + TASK + '.\n\n' + "The total plan:\n" + plan_str + "\n\n" + "The step {} is:\n".format(str(idx+1)) + sub_step# + "\n\n" + "Your goal: write the code for the step {}.'".format(str(idx+1))
                if previous_step == "":
                    message2 ="Step {} is:\n".format(str(idx+1)) + sub_step# + "\n\n" + "Your goal: write the code for the step {}.'".format(str(idx+1))
                else:
                    #message2 ="Previous steps are:\n"+ previous_step + "\n" + "Codes of previous steps are:\n"  + tmp  + "\n\n" +"Step {} is:\n".format(str(idx+1)) + sub_step
                    message2 ="Codes of previous steps are:\n"  + tmp  + "\n\n" +"Step {} is:\n".format(str(idx+1)) + sub_step

                code_prompt = get_prompt_coder_hard(message2, goals, None)
                
                
                PROMPT_CODE = """Now, you're a retrieve augmented Coder. You goal is to write the code for the current step task based on your own knowledge and the context. Context generally contains how specific packages are used. Thus, you need to extract useful information from the context.\nNote that you do not need to follow the code in Context strictly. In other words, You can modify or discard the lines of code in the context that are not necessary for the task, or you can write new code yourself.\nIf Context does not contain the knowledge you need to write the code, then ignore Context and write the code using only your own knowledge.\n\nContext is: {input_context}
                """

                retriever_code.customized_prompt = "the task of step {}:\n".format(str(idx+1)) + ori_sub_step + '\n\n' + PROMPT_CODE

                # == Begin coder ==
                MAX_ROUND = 100
                groupchat = GroupChat(agents=[code_proxy_hard, code_proxy_hard2, coder, checker, retriever_call, retriever_executor], messages=[], max_round=MAX_ROUND, speaker_selection_method=custom_speaker_selection_func_code)
                manager = GroupChatManager(groupchat=groupchat, llm_config=llm_config)
                
                if len(groupchat.messages) == indicator:
                    chat_history = await cl.make_async(code_proxy_hard2.initiate_chat)(manager, message=str(code_prompt)+ '\n\nretrieve content')
                    #await user_proxy.initiate_chat( manager, message=message)
                elif len(groupchat.messages) <  MAX_ROUND + indicator:
                    chat_history = await cl.make_async(code_proxy_hard.send)(manager, message=TASK,)
                elif len(groupchat.messages) ==  MAX_ROUND + indicator:  
                    chat_history = await cl.make_async(code_proxy_hard.send)(manager, message="exit",)
                # == End code ==
                
                #indicator = len(groupchat.messages)
                
                previous_step = previous_step + sub_step + "\n"

                
                
                for i in groupchat.messages:
                    if i['name'] == "Coder":
                        if "```python" in i['content']:
                            #pattern = r'```python(.*?)```'
                            #result = re.findall(pattern, i['content'], re.DOTALL)
                            #tmp_coder_messages = "#step {}\n".format(str(idx+1)) +"```python" +'\n'.join(result)+ "```" + '\n' #i['content']
                            tmp_coder_messages = i["content"]
                            
                #import pdb
                #pdb.set_trace()
                coder_messages.append(tmp_coder_messages)
                tmp = '\n'.join(coder_messages)

    task_ = "temp"
    tmp = f"/home/liuyang/llm-agent-finalcode-v8-finalversion/experiments/{task_}/"
    if os.path.exists(tmp):
        pass
    else:
        os.makedirs(tmp)

    f = open(f'/home/liuyang/llm-agent-finalcode-v8-finalversion/experiments/{task_}/{task_}.txt', 'w')
    f.write(str(groupchat.messages))
    f.close()

    

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

'''
Do a spatial transcription analysis.
/Users/liuyang/home/bioautogen/biodataset/output_merfish_m1s3.h5ad: the spatial transcription data of mice
'''

'''
use scanpy to find the differentially expressed genes
/Users/liuyang/home/autogen/Autogen-UI-main/dataset/filtered_gene_bc_matrices/hg19: path to 10x mtx data
'''

'''
/Users/liuyang/home/packages/TrimGalore-0.4.5/trim_galore \
    --quality 20 \
    --length 36 \
    --fastqc \
    ./dataset/data/SRR1374921.fastq.gz \
    ./dataset/data/TruSeq3-SE.fa
'''