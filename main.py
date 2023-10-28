sequence1 = "ATC"
sequence2 = "ATCC"

missmatchPenalty=-1
gapPenalty=-2
match=1

def createArray():
    array = [[0 for x in range(len(sequence1)+1)] for y in range(len(sequence2)+1)] 
    return array


def prettyPrint(array):
    nrow=0
    print('     ' + ', '.join(sequence1))
    for row in array:
        if(nrow==0):
            print(' '+str(row))
        else:
            print((sequence2[nrow-1])+str(row))

        nrow+=1
        
    

prettyPrint(createArray())