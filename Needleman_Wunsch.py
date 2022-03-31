from typing import Union, Tuple

def needleman_wunsch(string1:str,
                     string2:str,
                     match:Union[int, float]=1,
                     mismatch:Union[int, float]=-1,
                     indel:Union[int, float]=-1
                    ) -> Tuple[str, str]:

    scores = [[-i-j for i in range(len(string1)+1)] for j in range(len(string2)+1)]

    for idx in range(1, len(string1)+1):
        for jdx in range(1, len(string2)+1):
            candidates = []
            if string1[idx-1] == string2[jdx-1]:
                candidates.append(scores[jdx-1][idx-1] + match)
            else:
                candidates.append(scores[jdx-1][idx-1] + mismatch)
            candidates.append(scores[jdx-1][idx] + indel)
            candidates.append(scores[jdx][idx-1] + indel)
            scores[jdx][idx] = max(candidates)
    scores = [i[1:] for i in scores][1:]

    idx = len(string1)-1
    jdx = len(string2)-1

    newstr1 = string1[-1]
    newstr2 = string2[-1]

    while idx>0 and jdx>0:
        mv = [scores[jdx-1][idx], scores[jdx][idx-1], scores[jdx-1][idx-1]]
        f = lambda i: mv[i]
        index = max(range(len(mv)), key=f)
        assert index in [0,1,2]
        if index == 2:
            newstr1 = string1[idx-1] + newstr1
            newstr2 = string2[jdx-1] + newstr2
            idx -= 1
            jdx -= 1
        elif index == 1:
            newstr1 = string1[idx-1] + newstr1
            newstr2 = '-' + newstr2
            idx -= 1
        else:
            newstr1 = '-' + newstr1
            newstr2 = string2[jdx-1] + newstr2
            jdx -= 1
            
    return newstr1, newstr2

if __name__ == '__main__':
    sequence1 = input("Sequence 1 : ")
    sequence2 = input("Sequence 2 : ")
    aligned1, aligned2 = needleman_wunsch(sequence1, sequence2)
    print("Aligned 1 :", aligned1)
    print("Aligned 2 :", aligned2)