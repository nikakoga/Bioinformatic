sequence1 = "ATGCT"
sequence2 = "AGCT"

mismatchPenalty=-1
gapPenalty=-2
match=1

def createArray():
    array = [[0 for x in range(len(sequence1)+1)] for y in range(len(sequence2)+1)] 
    return array

def fillArray(array):
    for row in range(1,len(sequence2)+1):
        for col in range (1,len(sequence1)+1):
            topValue = array[row-1][col]+gapPenalty
            leftValue=array[row][col-1]+gapPenalty
            diagonalValue=countMatch(array,row,col)
            array[row][col]=max(correctIfNegative(topValue,leftValue,diagonalValue))
    return array

def countMatch(array,row,column):
    if(sequence1[column-1]==sequence2[row-1]):
        return array[row-1][column-1]+match
    else:
        return array[row-1][column-1]+mismatchPenalty

def correctIfNegative(number1,number2,number3):
    numbers=[number1,number2,number3]
    result=[]
    for number in numbers:
        if(number<0):
            result.append(0)
        else:
            result.append(number)
    return result

def prettyPrint(array):
    nrow=0
    print('     ' + ', '.join(sequence1))
    for row in array:
        if(nrow==0):
            print(' '+str(row))
        else:
            print((sequence2[nrow-1])+str(row))

        nrow+=1
        
    

prettyPrint(fillArray(createArray()))