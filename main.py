import sys
from Bio import SeqIO

def getNumberFromString(x):
    if x[0] == '-':
        number = x[1:]
        isNegative = True
    else:
        number = x
        isNegative = False
    if number.isnumeric():
        number=int(number)
        if isNegative:
            number = -1*number
        return True, number
    else:
        return False, 0
    

def getFastaSequences():
    firstSequence=""
    secondSequence=""
    isAlignmentPossible = True
    if "--fasta" in sys.argv and (sys.argv.index("--fasta")+1)<len(sys.argv):
        sequenceFilepath = sys.argv[sys.argv.index("--fasta")+1]
        try:
            file = open(sequenceFilepath, "r")
        except FileNotFoundError:
            print("The file you provided doesn't exist. Make sure the path and filename are correct.")
        sequences = list(SeqIO.parse(file, "fasta"))
        if len(sequences)<2:
            print("Please provide a fasta file with two sequences to run the algorithm.")
            isAlignmentPossible = False
        else:
            if len(sequences)>2:
                print("The file you provided has more than two sequences. The alignment will be created for the first two.")
            firstSequence = sequences[0].seq
            secondSequence = sequences[1].seq
    
    else:
        print("Please provide a .fasta file to run the algorithm.")
        isAlignmentPossible = False

    return isAlignmentPossible, firstSequence, secondSequence

def getMatchScore():
    if "--match" in sys.argv and (sys.argv.index("--match")+1)<len(sys.argv):
            matchScore = sys.argv[sys.argv.index("--match")+1]
            isNumeric, matchScore = getNumberFromString(matchScore)
            if isNumeric:
                print("Setting match score: {}".format(matchScore))
            else:
                print("The value for match score you provided is not a number. The default value will be used instead.")
                print("Setting match score: 1")
                matchScore=1
    else:
        print("Setting match score: 1")
        matchScore=1

    return matchScore

def getGapScore():
    if "--gap" in sys.argv and (sys.argv.index("--gap")+1)<len(sys.argv):
            gapScore = sys.argv[sys.argv.index("--gap")+1]
            isNumeric, gapScore = getNumberFromString(gapScore)
            if isNumeric:
                print("Setting gap score: {}".format(gapScore))
            else:
                print("The value for gap score you provided is not a number. The default value will be used instead.")
                print("Setting gap score: -1")
                gapScore=-1
    else:
        print("Setting gap score: -1")
        gapScore=-1

    return gapScore

def getMismatchScore():
    if "--mismatch" in sys.argv and (sys.argv.index("--mismatch")+1)<len(sys.argv):
            mismatchScore = sys.argv[sys.argv.index("--mismatch")+1]
            isNumeric, mismatchScore = getNumberFromString(mismatchScore)
            if isNumeric:
                print("Setting mismatch score: {}".format(mismatchScore))
            else:
                print("The value for mismatch score you provided is not a number. The default value will be used instead.")
                print("Setting mismatch score: -1")
                mismatchScore=-1
    else:
        print("Setting mismatch score: -1")
        mismatchScore=-1

    return mismatchScore

def createArray(firstSeqLength, secondSeqLength):
    array = [[0 for x in range(firstSeqLength)] for y in range(secondSeqLength)] 
    return array

def fillArray(array, sequence1, sequence2, match, mismatchPenalty, gapPenalty):
    for row in range(1,len(sequence2)+1):
        for col in range (1,len(sequence1)+1):
            topValue = array[row-1][col]+gapPenalty
            leftValue=array[row][col-1]+gapPenalty
            diagonalValue=array[row-1][col-1] + countMatch(sequence1[col-1],sequence2[row-1],match,mismatchPenalty)
            array[row][col]=max(topValue,leftValue,diagonalValue,0)
    return array

def countMatch(firstSeqChar,secondSeqChar,match,mismatchPenalty):
    if(firstSeqChar==secondSeqChar):
        return match
    else:
        return mismatchPenalty

def prettyPrint(array, sequence1, sequence2):
    nrow=0
    print('     ' + ', '.join(sequence1))
    for row in array:
        if(nrow==0):
            print(' '+str(row))
        else:
            print((sequence2[nrow-1])+str(row))

        nrow+=1
        
    
if __name__ == "__main__":
    isAlignmentPossible, firstSequence, secondSequence = getFastaSequences()

    if isAlignmentPossible:
        matchScore = getMatchScore()
        gapScore = getGapScore()
        mismatchScore = getMismatchScore()
        prettyPrint(fillArray(createArray(len(firstSequence),len(secondSequence)),firstSequence,secondSequence,matchScore,mismatchScore,gapScore))