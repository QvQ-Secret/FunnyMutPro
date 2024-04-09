from pymol import cmd
import numpy as np
from decimal import Decimal
import matplotlib.pyplot as plt
from Bio.Align import PairwiseAligner
from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap


def decimal_two(number:float) -> Decimal:
    return Decimal(number).quantize(Decimal("0.00"))


def get_alignment(querySeq,matchSeq):
        aligner = PairwiseAligner(scoring='blastp')
        aligner.mode = 'global'
        aligner.match_score = 2
        if '+' not in querySeq:
            aligner.mismatch_score = -0.7
        else:
            aligner.mismatch_score = -0.5
        aligner.open_gap_score = -3
        aligner.extend_gap_score = -0.1
        aligner.target_end_gap_score = 0
        aligner.query_end_gap_score = 0
        alignment = aligner.align(querySeq, matchSeq)[0]
        return alignment

_1ENH = 'RPRTAFSSEQLARLKREFNENRYLTERRRQQLSSELGLNEAQIKIWFQNKRAKI' # 3, 56
_1ZTR = 'GDEKRPRTAFSSEQLARAKREFNENRYLTERRRQQLSSELGLNEAQIKIWFQNKRAKIRRS' # -1, 59
alignment = get_alignment(_1ENH,_1ZTR)
# print(alignment)
indexQuery,indexMatch = [],[]
length = len(alignment.query)
startPositionTarget = 3
startPositionQuery = -1
for index in range(length):
    t = alignment[0][index]
    q = alignment[1][index]
    if t != '-':
        indexQuery.append(str(startPositionTarget))
        startPositionTarget += 1
    else:
        indexTarget.append('-')
    if m != '-':
        indexQuery.append(str(startPositionQuery))
        startPositionQuery += 1
    else:
        indexQuery.append('-')
for _ in range(length):
    print(_,'\t\t',indexTarget[_],'\t\t',indexQuery[_])

a = indexTarget.index('16')
b = indexQuery[a]
print(a,b,alignment[0][a],alignment[1][a])

def heatMap(targetPath,queryPath,resiPos):
    cmd.delete('all')
    cmd.load(targetPath,'prot')
    cmd.select('subs',f'resi {resiPos}')
    cmd.select('pocket', 'byres subs around 5')
    srcDictTarget = get_src_dict('pocket')
    resi = srcDictTarget['resi']
    # resn = srcDictTarget['resn']
    resiT, resiQ = [], []
    xLabel = []
    for residue in resi:
        try:
            index = indexTarget.index(residue)
            if alignment[0][index] != '-' and alignment[1][index] != '-':
                xLabel.append(f"{residue} {alignment[0][index]}")
                resiT.append(indexTarget[index])
                resiQ.append(indexQuery[index])
        except ValueError:
            pass

    shape = len(resiT)
    disMatrixT = np.zeros((shape,shape))
    for i in range(shape):
        for j in range(i, shape):
            disMatrixT[i, j] = decimal_two(cmd.distance(f'{resiT[i]}/CA',f'{resiT[j]}/CA'))
            disMatrixT[j, i] = disMatrixT[i, j]
    cmd.delete('all')
    cmd.load(queryPath,'prot')
    disMatrixQ = np.zeros((shape,shape))
    for i in range(shape):
        for j in range(i, shape):
            disMatrixQ[i, j] = decimal_two(cmd.distance(f'{resi[i]}/CA',f'{resi[j]}/CA'))
            disMatrixQ[j, i] = disMatrixM[i, j]
    cmd.delete('all')
    np.set_printoptions(precision=2,suppress=True)
    data = disMatrixQ-disMatrixT
    norm = Normalize(vmin=-10, vmax=10)
    cmap = LinearSegmentedColormap.from_list('rgw', ['cornflowerblue','w','tomato'], 512)
    ax = plt.imshow(data, cmap=cmap,norm=norm)

    ax = plt.gca()

    ax.xaxis.set_ticks_position('top') # 设置x轴到上方
    ax.yaxis.set_ticks_position('left') # 设置y轴到左方
    ax.tick_params(axis='x', rotation=45)

    # plt.imshow(data, cmap='OrRd',norm=norm)
    plt.xticks([i for i in range(shape)],xLabel)
    plt.yticks([i for i in range(shape)],xLabel)
    plt.colorbar()
    plt.show()

heatMap(targetPath=r'1ZTR_A.pdb',queryPath=r'1ENH_A.pdb',resiPos='16')
