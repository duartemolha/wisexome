import _pickle as cPickle
import argparse
import glob
import time
import numpy
from tqdm import tqdm


def quickSelect(bigList, amount):
    bestList = [[1]] * amount
    tmpList = [[1]] * amount
    curV = 0
    worst = 1
    for iR, valR in enumerate(bigList):
        if curV == amount:
            curV = 0
            bestList = sorted(bestList + tmpList)[:amount]
            worst = max([x[0] for x in bestList])
        if valR < worst:
            tmpList[curV] = [valR, iR]
            curV += 1
    return sorted(bestList + tmpList[:curV])[:amount]


def notQuickSelect(bigList, amount):
    tmpList = [[val, i] for i, val in enumerate(bigList)]
    return sorted(tmpList)[:amount]


parser = argparse.ArgumentParser(
    description="TBA (or CBA perhaps)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("refdir", type=str, help="files containg read hits per probe")
parser.add_argument("dropfile", type=str, help="file containg not too different probes")
parser.add_argument("targetchr", type=str, help="chromosome to target bins on")
parser.add_argument("refchr", type=str, help="chromosome to find reference bins on")
args = parser.parse_args()

refData = args.refdir
fileDrop = args.dropfile
tChr = args.targetchr
rChr = args.refchr

distances = []

targets = []
references = []

referenceFiles = glob.glob(refData + "/*.hits")
for refFile in referenceFiles:
    print("\tLoading:\t" + refFile)
    hitFile = cPickle.load(open(refFile, "r"))
    sumCount = float(
        sum([sum(hitFile[x]) for x in hitFile.keys()])
    )  # -sum(hitFile['chrX']) # - X is newly added after ex12
    target = [x / sumCount for x in hitFile[tChr]]
    reference = [x / sumCount for x in hitFile[rChr]]

    if tChr == "chrX" and numpy.median(target) < 0.00000175:  # 0.000002:
        target = [x / sumCount * 2 for x in hitFile[tChr]]
        # Well nevermind but technically we could accept the male samples for training
    else:
        targets.append(target)
        references.append(reference)


tarsT = map(list, zip(*targets))
refsT = map(list, zip(*references))

print(len(tarsT), len(tarsT[0]))
print(len(refsT), len(refsT[0]))

t = time.time()
toDoLen = len(tarsT)
distances = [[]] * len(tarsT)
# diffLists = [[]] * len(tarsT)
for tarTi, tarTv in enumerate(tqdm(tarsT)):
    curDists = [0] * len(refsT)
    for refTi, refTv in enumerate(refsT):
        curDist = 0
        for i in range(len(tarTv)):
            diff = tarTv[i] - refTv[i]
            curDist += diff * diff
        curDists[refTi] = curDist
    distances[tarTi] = quickSelect(curDists, 100)

cPickle.dump(distances, open(fileDrop, "wb"))
