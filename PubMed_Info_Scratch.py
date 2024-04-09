import re
import time

import pandas as pd
import requests
from bs4 import BeautifulSoup as bf
from tqdm import tqdm


headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/113.0.0.0 Safari/537.36",
}


def getInfo(pubMedID: str) -> pd.DataFrame:
    url_new = f'https://pubmed.ncbi.nlm.nih.gov/{pubMedID}'
    article_page = bf(requests.get(url_new, headers=headers).content, 'html.parser').find('div', id='article-page')
    try:
        journal_shortname = article_page.find('button', id='full-view-journal-trigger').text.strip()  # 期刊缩写
        journal_longname = article_page.find('div', id='full-view-journal').get('aria-label').replace(
            'Dropdown menu for journal ', '')  # 期刊全称
    except:
        journal_shortname = '/'
        journal_longname = '/'
    article_title = article_page.find('h1', class_='heading-title').text.strip()  # 文献题目
    try:
        authors_get = article_page.find('div', class_='authors-list').find_all('span')
        authors_join = ','.join(re.findall('data-ga-label=\"(.*?)\"', str(authors_get)))  # 作者列表
    except:
        authors_join = '/'
    try:
        abstract = article_page.find('div', class_='abstract-content').text.strip()  # 文献摘要
        # trans = youDao(abstract)
    except:
        abstract = '/'
        trans = '/'
    try:
        cit = article_page.find('span', class_='cit').text
    except:
        cit = '/'
    try:
        secondary_date = article_page.find('span', class_='secondary-date').text.strip()
    except:
        secondary_date = '/'
    try:
        doi = article_page.find('span', class_='identifier doi').find('a', class_='id-link').get('href')
    except:
        doi = '/'
    new_df = pd.DataFrame(
        {
            'PMID': [pubMedID],
            'Title': [article_title],
            'Journal Short': [journal_shortname],
            'Journal Long': [journal_longname],
            'Authors': [authors_join],
            'Abstract': [abstract],
            'Cit': [cit],
            'Secondary-date': [secondary_date],
            'Article Url': [url_new],
            'DOI': [doi],
        }
    )
    return new_df


df = pd.read_excel(r"PMID.xlsx", index_col=None)
pubmed_list = list(df['PMID'])
result_df = pd.DataFrame()
for index, pmid in tqdm(enumerate(pubmed_list)):
    if index > -1:
        try:
            result_df = pd.concat([result_df, getInfo(pmid)], ignore_index=True)
        except:
            continue
        time.sleep(3)
        if index != 0 and index % 100 == 0:
            result_df.to_excel(f"ref_result_{index}.xlsx")
# print(result_df)
result_df.to_excel(r"ref_result.xlsx")
