# -*- coding: UTF-8 -*-
# 有点搞，本来基准叫做query，匹配内容叫做match，但是好像别人都叫target和query，现在就是最后的结果改成了target和query，但是中间的变量命名全都没有变
# 太坏了
from Bio.Align import PairwiseAligner


class Align:
    def __init__(self, querySeq, matchSeq):
        self.querySeq = querySeq
        self.matchSeq = matchSeq
        self.alignment = self.alignment()
        self.query = self.alignment[0]
        self.match = self.alignment[1]
        self.alignRange = self.alignPointGetter
        self.alignS = self.alignRange[0]
        self.alignE = self.alignRange[1]

    def alignment(self):
        aligner = PairwiseAligner(scoring='blastp')
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -0.5
        aligner.open_gap_score = -3
        aligner.extend_gap_score = -0.1
        aligner.target_end_gap_score = 0
        aligner.query_end_gap_score = 0
        alignments = aligner.align(self.querySeq, self.matchSeq)[0]
        return alignments[0], alignments[1], alignments

    @property
    def alignPointGetter(self):
        startPoint = None
        endPoint = None
        for _, __ in enumerate(self.query):
            if self.query[_] == self.match[_]:
                startPoint = _
                break
        queryReverse = ''.join(reversed(self.query))
        matchReverse = ''.join(reversed(self.match))
        for _, __ in enumerate(queryReverse):
            if queryReverse[_] == matchReverse[_]:
                endPoint = _
                break
        endPoint = -endPoint
        points = [startPoint, endPoint]
        for _ in range(2):
            if points[_] == 0:
                points[_] = None
        return points

    @property
    def queryAlign(self):
        return self.query[self.alignS:self.alignE]

    @property
    def matchAlign(self):
        return self.match[self.alignS:self.alignE]

    @property
    def queryInit(self):
        if self.alignS is not None:
            return self.query[:self.alignS]
        else:
            return ''

    @property
    def matchInit(self):
        if self.alignS is not None:
            return self.match[:self.alignS]
        else:
            return ''

    @property
    def queryTer(self):
        if self.alignE is not None:
            return self.query[self.alignE:]
        else:
            return ''

    @property
    def matchTer(self):
        if self.alignE is not None:
            return self.match[self.alignE:]
        else:
            return ''

    @property
    def alignLen(self):
        return len(self.queryAlign)

    # map residue index
    @staticmethod
    def resiFill(srcResi: list, srcSeq) -> list:
        dstResi = []
        cnt = 0
        for i in srcSeq:
            if i != '-':
                dstResi.append(srcResi[cnt])
                cnt += 1
            else:
                dstResi.append('/')
        return dstResi

    def queryAlignResi(self, resi: list) -> list:
        """

        param resi: the source resi
        :return: the resi slice of query protein within the alignment range
        """
        filledResi = self.resiFill(resi, self.query)
        return filledResi[self.alignS:self.alignE]

    def matchAlignResi(self, resi: list) -> list:
        """

        param resi: the source resi
        :return: the resi slice of match protein within the alignment range
        """
        filledResi = self.resiFill(resi, self.match)
        return filledResi[self.alignS:self.alignE]

    # map residue name
    @staticmethod
    def resnFill(srcResn: list, srcSeq) -> list:
        dstResn = []
        cnt = 0
        for i in srcSeq:
            if i != '-':
                dstResn.append(srcResn[cnt])
                cnt += 1
            else:
                dstResn.append('/')
        return dstResn

    def queryAlignResn(self, resn: list) -> list:
        """

        param resn: the source resn
        :return: the resn slice of query protein within the alignment range
        """
        filledResn = self.resiFill(resn, self.query)
        return filledResn[self.alignS:self.alignE]

    def matchAlignResn(self, resn: list) -> list:
        """

        param resn: the source resn
        :return: the resn slice of match protein within the alignment range
        """
        filledResn = self.resiFill(resn, self.match)
        return filledResn[self.alignS:self.alignE]

    # map secondary structure
    @staticmethod
    def ssFill(ss: list, srcSeq) -> list:
        dstSs = []
        cnt = 0
        for i in srcSeq:
            if i != '-':
                dstSs.append(ss[cnt])
                cnt += 1
            else:
                dstSs.append('/')
        return dstSs

    @staticmethod
    def formatSs(ssList):
        ssSet = set(ssList)
        if 'L' in ssSet:
            ssSet.remove('L')
        ss = [i[0] for i in ssSet]
        helixCnt = ss.count('H')
        sheetCnt = ss.count('S')
        ssFmt = ''
        if helixCnt != 0:
            if helixCnt == 1:
                ssFmt += 'α'
            else:
                ssFmt += f'{helixCnt}α'
        if sheetCnt != 0:
            if ssFmt != '':
                ssFmt += '+'
            if sheetCnt == 1:
                ssFmt += 'β'
            else:
                ssFmt += f'{sheetCnt}β'
        if ssFmt == '':
            ssFmt = '/'
        return ssFmt

    def queryAlignSs(self, ss: list) -> str:
        """

        param resn: the source resn
        :return: the resn slice of query protein within the alignment range
        """
        filledSs = self.ssFill(ss, self.query)
        alignSs = filledSs[self.alignS:self.alignE]
        ssFmt = self.formatSs(alignSs)
        return ssFmt

    def matchAlignSs(self, ss: list) -> str:
        """

        param resn: the source resn
        :return: the resn slice of match protein within the alignment range
        """
        filledSs = self.ssFill(ss, self.match)
        alignSs = filledSs[self.alignS:self.alignE]
        ssFmt = self.formatSs(alignSs)
        return ssFmt

    def alignVari(self, queryResi, matchResi, queryResn: list, matchResn: list):
        """
        Align区域的序列变化
        
        由于对应的内容会导致产生的文件较大，故暂时注释掉Resn部分
        :return:
        """
        alignDict = {
            'targetSub': [],
            'targetDelins': [],
            'targetDel': [],
            'targetIns': [],
            'querySub': [],
            'queryDelins': [],
            'queryDel': [],
            'queryIns': [],
        }
        queryAlignResi = self.queryAlignResi(queryResi)
        matchAlignResi = self.matchAlignResi(matchResi)
        # queryAlignResn = self.queryAlignResn(queryResn)
        # matchAlignResn = self.matchAlignResn(matchResn)
        sliceDict = {}
        count = 0
        querySliceAmino, matchSliceAmino, querySliceResi, matchSliceResi = [], [], [], []
        # querySliceResn, matchSliceResn = [], []
        for _ in range(self.alignLen):
            queryAmino = self.queryAlign[_]
            matchAmino = self.matchAlign[_]
            if queryAmino != matchAmino:
                flag = False
            else:
                flag = True
            if not flag:
                querySliceAmino.append(queryAmino)
                matchSliceAmino.append(matchAmino)
                querySliceResi.append(queryAlignResi[_])
                matchSliceResi.append(matchAlignResi[_])
                # querySliceResn.append(queryAlignResn[_])
                # matchSliceResn.append(matchAlignResn[_])
            elif len(querySliceResi) > 0:
                sliceDict[count] = {
                    'targetAmino': querySliceAmino,
                    'queryAmino': matchSliceAmino,
                    'targetResi': querySliceResi,
                    'queryResi': matchSliceResi,
                    # 'queryResn': querySliceResn,
                    # 'matchResn': matchSliceResn,
                }
                querySliceAmino, matchSliceAmino, querySliceResi, matchSliceResi = [], [], [], []
                # querySliceResn, matchSliceResn = [], []
                count += 1
        for values in sliceDict.values():
            # 删掉比对中为空的部分
            if True:
                while '-' in values['targetAmino']:
                    values['targetAmino'].remove('-')
                while '-' in values['queryAmino']:
                    values['queryAmino'].remove('-')
                while '/' in values['targetResi']:
                    values['targetResi'].remove('/')
                while '/' in values['queryResi']:
                    values['queryResi'].remove('/')
                # while '/' in values['queryResn']:
                #     values['queryResn'].remove('/')
                # while '/' in values['queryResn']:
                #     values['queryResn'].remove('/')
            queryBool = len(values['targetAmino'])
            matchBool = len(values['queryAmino'])
            if queryBool and matchBool:
                if len(values['targetAmino']) == 1:
                    if len(values['queryAmino']) == 1:
                        # query, match的序列变化均为1，此时为Sub
                        alignDict['targetSub'].append(f"{values['targetAmino'][0]}{values['targetResi'][0]}{values['queryAmino'][0]}")
                        alignDict['querySub'].append(f"{values['queryAmino'][0]}{values['queryResi'][0]}{values['targetAmino'][0]}")
                        # alignDict['querySub'].append(f"{values['queryResn']}{values['queryResi']}{values['matchResn']}")
                        # alignDict['matchSub'].append(f"{values['matchResn']}{values['matchResi']}{values['queryResn']}")
                    else:
                        # query为1，match>1,此时为delins
                        alignDict['targetDelins'].append(
                            f"{values['targetAmino'][0]}{values['targetResi'][0]}delins{''.join(values['queryAmino'])}"
                        )
                        alignDict['queryDelins'].append(
                            f"{values['queryAmino'][0]}{values['queryResi'][0]}_{values['queryAmino'][-1]}{values['queryResi'][-1]}"
                            f"delins"
                            f"{values['targetAmino'][0]}"
                        )
                        # alignDict['queryDelins'].append(
                        #     f"{values['queryResn'][0]}{values['queryResi'][0]}delins{''.join(values['matchResn'])}"
                        # )
                        # alignDict['matchDelins'].append(
                        #     f"{values['matchResn'][0]}{values['matchResi'][0]}_{values['matchResn'][-1]}{values['matchResi'][-1]}"
                        #     f"delins"
                        #     f"{values['queryResn'][0]}"
                        # )
                else:
                    if len(values['targetAmino']) == 1:
                        # query>1, match = 1, 此时为delins
                        alignDict['targetDelins'].append(
                            f"{values['targetAmino'][0]}{values['targetResi'][0]}_{values['targetAmino'][-1]}{values['targetResi'][-1]}"
                            f"delins"
                            f"{values['queryAmino'][0]}"
                        )
                        alignDict['queryDelins'].append(
                            f"{values['queryAmino'][0]}{values['queryResi'][0]}delins{''.join(values['targetAmino'])}"
                        )
                        # alignDict['queryDelins'].append(
                        #     f"{values['queryResn'][0]}{values['queryResi'][0]}_{values['queryResn'][-1]}{values['queryResi'][-1]}"
                        #     f"delins"
                        #     f"{values['matchResn'][0]}"
                        # )
                        # alignDict['matchDelins'].append(
                        #     f"{values['matchResn'][0]}{values['matchResi'][0]}delins{''.join(values['queryResn'])}"
                        # )
                    else:
                        # query, match > 1, 此时为delins
                        alignDict['targetDelins'].append(
                            f"{values['targetAmino'][0]}{values['targetResi'][0]}_{values['targetAmino'][-1]}{values['targetResi'][-1]}"
                            f"delins"
                            f"{''.join(values['queryAmino'])}"
                        )
                        alignDict['queryDelins'].append(
                            f"{values['queryAmino'][0]}{values['queryResi'][0]}_{values['queryAmino'][-1]}{values['queryResi'][-1]}"
                            f"delins"
                            f"{''.join(values['targetAmino'])}"
                        )
                        # alignDict['queryDelins'].append(
                        #     f"{values['queryResn'][0]}{values['queryResi'][0]}_{values['queryResn'][-1]}{values['queryResi'][-1]}"
                        #     f"delins"
                        #     f"{''.join(values['matchResn'])}"
                        # )
                        # alignDict['matchDelins'].append(
                        #     f"{values['matchResn'][0]}{values['matchResi'][0]}_{values['matchResn'][-1]}{values['matchResi'][-1]}"
                        #     f"delins"
                        #     f"{''.join(values['queryResn'])}"
                        # )
            elif queryBool and not matchBool:
                matchResiInit = queryAlignResi.index(values['targetResi'][0]) - 1
                if len(values['targetAmino']) == 1:
                    alignDict['targetDel'].append(f"{values['targetAmino'][0]}{values['targetResi'][0]}del")
                    matchResiTer = queryAlignResi.index(values['targetResi'][0]) + 1
                else:
                    alignDict['targetDel'].append(f"{values['targetAmino'][0]}{values['targetResi'][0]}"
                                                  f"_"
                                                  f"{values['targetAmino'][-1]}{values['targetResi'][-1]}"
                                                  f"del")
                    matchResiTer = queryAlignResi.index(values['targetResi'][-1]) + 1
                alignDict['queryIns'].append(f"{self.matchAlign[matchResiInit]}{matchAlignResi[matchResiInit]}"
                                             f"_"
                                             f"{self.matchAlign[matchResiTer]}{matchAlignResi[matchResiTer]}"
                                             f"ins{''.join(values['queryAmino'])}")
            elif not queryBool and matchBool:
                queryResiInit = matchAlignResi.index(values['queryResi'][0]) - 1
                if len(values['queryAmino']) == 1:
                    alignDict['queryDel'].append(f"{values['queryAmino'][0]}{values['queryResi'][0]}del")
                    queryResiTer = matchAlignResi.index(values['queryResi'][0]) + 1
                else:
                    alignDict['queryDel'].append(f"{values['queryAmino'][0]}{values['queryResi'][0]}"
                                                 f"_"
                                                 f"{values['queryAmino'][-1]}{values['queryResi'][-1]}"
                                                 f"del")
                    queryResiTer = matchAlignResi.index(values['queryResi'][-1]) + 1
                alignDict['targetIns'].append(f"{self.queryAlign[queryResiInit]}{queryAlignResi[queryResiInit]}"
                                              f"_"
                                              f"{self.queryAlign[queryResiTer]}{queryAlignResi[queryResiTer]}"
                                              f"ins{''.join(values['queryAmino'])}")
        for keys, values in alignDict.items():
            if len(values) == 0:
                alignDict[keys] = '/'
            else:
                alignDict[keys] = '; '.join(alignDict[keys])
        # print(alignDict)
        return alignDict

    def azylVari(self):
        """
        N端的序列变化
        :return:
        """
        variDict = {
            'targetExtInit': '/',
            'targetDelInit': '/',
            'targetDelinsInit': '/',
            'targetSubInit': '/',
            'queryExtInit': '/',
            'queryDelInit': '/',
            'queryDelinsInit': '/',
            'querySubInit': '/'
        }

        queryInit = self.queryInit.replace('-', '')
        matchInit = self.matchInit.replace('-', '')
        # 如果此处的布尔值为真，则对应的字符串不为空。
        queryBool = bool(len(queryInit))
        matchBool = bool(len(matchInit))
        if queryBool and matchBool:
            if len(queryInit) == 1 and len(matchInit) == 1:
                variDict['targetSubInit'] = f'{queryInit}subInit{matchInit}'
                variDict['querySubInit'] = f'{matchInit}subInit{queryInit}'
            else:
                variDict['targetDelinsInit'] = f'{queryInit}delinsInit{matchInit}'
                variDict['queryDelinsInit'] = f'{matchInit}delinsInit{queryInit}'
        elif queryBool and not matchBool:
            variDict['targetDelInit'] = f'{queryInit}delInit'
            variDict['queryExtInit'] = f'{queryInit}extInit'
        elif not queryBool and matchBool:
            variDict['targetExtInit'] = f'{matchInit}extInit'
            variDict['queryDelInit'] = f'{matchInit}delInit'
        return variDict

    def carboxylVari(self):
        """
        C端的序列变化
        :return:
        """
        variDict = {
            'targetExtTer': '/',
            'targetDelTer': '/',
            'targetDelinsTer': '/',
            'targetSubTer': '/',
            'queryExtTer': '/',
            'queryDelTer': '/',
            'queryDelinsTer': '/',
            'querySubTer': '/'
        }

        queryTer = self.queryTer.replace('-', '')
        matchTer = self.matchTer.replace('-', '')
        # 如果此处的布尔值为真，则对应的字符串不为空。
        queryBool = bool(len(queryTer))
        matchBool = bool(len(matchTer))
        if queryBool and matchBool:
            if len(queryTer) == 1 and len(matchTer) == 1:
                variDict['targetSubTer'] = f'{queryTer}subTer{matchTer}'
                variDict['querySubTer'] = f'{matchTer}subTer{queryTer}'
            else:
                variDict['targetDelinsTer'] = f'{queryTer}delinsTer{matchTer}'
                variDict['queryDelinsTer'] = f'{matchTer}delinsTer{queryTer}'
        elif queryBool and not matchBool:
            variDict['targetDelTer'] = f'{queryTer}delTer'
            variDict['queryExtTer'] = f'{queryTer}extTer'
        elif not queryBool and matchBool:
            variDict['targetExtTer'] = f'{matchTer}extTer'
            variDict['queryDelTer'] = f'{matchTer}delTer'
        return variDict


if __name__ == '__main__':
    myAlign = Align(
        'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKLGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG',
        'VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQGG'
    )
    print('2MGC原始序列:\n', myAlign.query)
    print('5OJ9原始序列:\n', myAlign.match)
    print(myAlign.alignment[2])
    # print(myAlign.alignPointGetter)
    # print(myAlign.queryInit, myAlign.queryAlign, myAlign.queryTer)
    # print(myAlign.matchInit, myAlign.matchAlign, myAlign.matchTer)
    print('序列匹配结果：')
    print('起始段结果：')
    print(myAlign.azylVari())
    print('匹配段结果：')
    print(myAlign.alignVari(
        queryResi=['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21',
                   '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41',
                   '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61',
                   '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81',
                   '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99', '100', '101',
                   '102', '103', '104', '105', '106', '107', '108', '109', '110', '111', '112', '113', '114', '115', '116', '117', '118',
                   '119', '120', '121', '122', '123', '124', '125', '126', '127', '128', '129', '130', '131', '132', '133', '134', '135',
                   '136', '137', '138', '139', '140', '141', '142', '143', '144', '145', '146', '147', '148', '149', '150', '151', '152',
                   '153'],
        matchResi=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21',
                   '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41',
                   '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61',
                   '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81',
                   '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '94', '95', '96', '97', '98', '99', '100', '101',
                   '102', '103', '104', '105', '106', '107', '108', '109', '110', '111', '112', '113', '114', '115', '116', '117', '118',
                   '119', '120', '121', '122', '123', '124', '125', '126', '127', '128', '129', '130', '131', '132', '133', '134', '135',
                   '136', '137', '138', '139', '140', '141', '142', '143', '144', '145', '146', '147', '148', '149', '150', '151', '152',
                   '153', '154'],
        queryResn=['GLY', 'ASP', 'GLU', 'LYS', 'ARG', 'PRO', 'ARG', 'THR', 'ALA', 'PHE', 'SER', 'SER', 'GLU', 'GLN', 'LEU', 'ALA', 'ARG',
                   'ALA',
                   'LYS', 'ARG', 'GLU', 'PHE', 'ASN', 'GLU', 'ASN', 'ARG', 'TYR', 'LEU', 'THR', 'GLU', 'ARG', 'ARG', 'ARG', 'GLN', 'GLN',
                   'LEU',
                   'SER', 'SER', 'GLU', 'LEU', 'GLY', 'LEU', 'ASN', 'GLU', 'ALA', 'GLN', 'ILE', 'LYS', 'ILE', 'TRP', 'PHE', 'GLN', 'ASN',
                   'LYS',
                   'ARG', 'ALA', 'LYS', 'ILE', 'ARG', 'ARG', 'SER'],
        matchResn=['GLY', 'ALA', 'MET', 'GLU', 'LYS', 'ARG', 'PRO', 'ARG', 'THR', 'ALA', 'PHE', 'SER', 'SER', 'GLU', 'GLN', 'LEU', 'ALA',
                   'ARG',
                   'LEU', 'LYS', 'ARG', 'GLU', 'PHE', 'ASN', 'GLU', 'ASN', 'ARG', 'TYR', 'LEU', 'THR', 'GLU', 'ARG', 'ARG', 'ARG', 'GLN',
                   'GLN',
                   'LEU', 'SER', 'SER', 'GLU', 'LEU', 'GLY', 'LEU', 'ASN', 'GLU', 'ALA', 'GLN', 'ILE', 'LYS', 'ILE', 'TRP', 'PHE', 'GLN',
                   'ASN',
                   'GLU', 'ARG', 'ALA', 'LYS', 'ILE', 'LYS', 'LYS', 'SER', 'GLY', 'SER']
    ))
    print('终止段结果：')
    print(myAlign.carboxylVari())
