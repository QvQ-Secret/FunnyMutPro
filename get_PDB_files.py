import json
import os
import pickle
import time

import requests
from tqdm import tqdm


def ding_msg(text_notice):
    dingTalkToken = '?' # use your own token
    url = f'https://oapi.dingtalk.com/robot/send?access_token={dingTalkToken}'
    txt = json.dumps(
        {
            'msgtype': 'text',
            'text': {
                'content': f'>_< 内容抓取失败，对应位点为{text_notice}'
            }
        }
    )
    headers = {
        'content-type': 'application/json',
        'User_Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/100.0.4896.127 Safari/537.36'
    }
    req = requests.post(url, data=txt, headers=headers)


def data_req(content):
    try:
        data = requests.get(f'https://www.rcsb.org/fasta/entry/{content}/download').text
    except:
        time.sleep(5)
        data = requests.get(f'https://www.rcsb.org/fasta/entry/{content}/download').text
    return data


os.chdir(r'G:\RCSB_PDB')

with open('PDB_files.pickle', 'rb') as f:
    currentPDB = pickle.load(f)

latestPDB_url = 'https://data.rcsb.org/rest/v1/holdings/current/entry_ids'
latestPDB = eval(requests.get(latestPDB_url).text)

additionalPDB = list(set(latestPDB) - set(currentPDB))
addition_index = len(currentPDB)
currentPDB.extend(additionalPDB)

with open('PDB_files.pickle', 'wb') as f:
    pickle.dump(currentPDB, f)

for index, content in tqdm(enumerate(currentPDB)):
    if index > addition_index:
        print(index)
        try:
            path = f'package_{str(index // 2000)}'
            folder = os.path.exists(path)
            if not folder:
                os.makedirs(path)
            with open(fr'{path}/{content}.fasta', 'w') as f:
                f.write(data_req(content))
            time.sleep(1)
            print(index, content, 'Success')
        except:
            text_notice = f'{index},{content}'
            ding_msg(text_notice)
            print(index, content, 'Fail')