import sys
from Bio import SeqIO

class textColors:
    WARNING = "\033[0;33;40m"
    ERROR = "\033[0;31;40m"
    RESET = "\033[0;0m"

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
            print(f"{textColors.ERROR}Error! The file you provided doesn't exist. Make sure the path and filename are correct.{textColors.RESET}")
        sequences = list(SeqIO.parse(file, "fasta"))
        if len(sequences)<2:
            print("Error! Please provide a fasta file with two sequences to run the algorithm.")
            isAlignmentPossible = False
        else:
            if len(sequences)>2:
                print(f"{textColors.WARNING}Warning: The file you provided has more than two sequences. The alignment will be created for the first two.{textColors.RESET}")
            firstSequence = sequences[0].seq
            secondSequence = sequences[1].seq
    
    else:
        print(f"{textColors.ERROR}Error! Please provide a .fasta file to run the algorithm.{textColors.RESET}")
        isAlignmentPossible = False

    return isAlignmentPossible, firstSequence, secondSequence

def getMatchScore():
    if "--match" in sys.argv and (sys.argv.index("--match")+1)<len(sys.argv):
            matchScore = sys.argv[sys.argv.index("--match")+1]
            isNumeric, matchScore = getNumberFromString(matchScore)
            if isNumeric:
                print("Setting match score: {}".format(matchScore))
            else:
                print(f"{textColors.WARNING}Warning: The value for match score you provided is not a number. The default value will be used instead.{textColors.RESET}")
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
                print(f"{textColors.WARNING}Warning: The value for gap score you provided is not a number. The default value will be used instead.{textColors.RESET}")
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
            print(f"{textColors.WARNING}Warning: The value for mismatch score you provided is not a number. The default value will be used instead.{textColors.RESET}")
            print("Setting mismatch score: -1")
            mismatchScore=-1
    else:
        print("Setting mismatch score: -1")
        mismatchScore=-1

    return mismatchScore

def getOutputFile():
    if "--output" in sys.argv and (sys.argv.index("--output")+1)<len(sys.argv):
        outputFile = sys.argv[sys.argv.index("--output")+1]
        if not outputFile.endswith(".txt"):
            outputFile = outputFile + ".txt"
    else:
        outputFile = "local_alignment.txt"
    return outputFile

def createArray(firstSeqLength, secondSeqLength):
    array = [[0 for x in range(firstSeqLength+1)] for y in range(secondSeqLength+1)] 
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

def getMaxFromMatrix(matrix, firstSeqLength, secondSeqLength):
    row, column, maximum = 0,0,0
    for i in range(1, secondSeqLength+1):
        for j in range(1, firstSeqLength+1):
            if matrix[i][j]>maximum:
                row, column, maximum = i, j, matrix[i][j]
    return row, column, maximum

def oneStepTraceback(matrix, row, column, gap, firstSequence, secondSequence):
    if matrix[row][column]==matrix[row-1][column]+gap:
        return row-1, column, matrix[row-1][column], "_", secondSequence[row-1]
    elif matrix[row][column]==matrix[row][column-1]+gap:
        return row, column-1, matrix[row][column-1], firstSequence[column-1], "_"
    else:
        return row-1, column-1, matrix[row-1][column-1], firstSequence[column-1], secondSequence[row-1]

def getAlignmentFromMatrix(matrix, firstSequence, secondSequence, gap):
    firstSeqAlignment, secondSeqAlignment = "", ""
    row, column, currentVal = getMaxFromMatrix(matrix, len(firstSequence), len(secondSequence))
    while currentVal != 0:
        row, column, currentVal, firstSeqChar, secondSeqChar=oneStepTraceback(matrix, row, column, gap, firstSequence, secondSequence)
        secondSeqAlignment, firstSeqAlignment = secondSeqChar+secondSeqAlignment, firstSeqChar+firstSeqAlignment

    print("Writing to file: {}".format(getOutputFile()))
    try:
        outputFile = open(getOutputFile(),"w")
        outputFile.write("Best local alignment:\n")
        outputFile.write(firstSeqAlignment + "\n")
        outputFile.write(secondSeqAlignment + "\n")
        outputFile.write("Score: {}".format(getMaxFromMatrix(matrix, len(firstSequence), len(secondSequence))[2]))
    except FileNotFoundError:
        print(f"{textColors.ERROR}Error! Make sure that path to the output file is correct.{textColors.RESET}")
        
    
if __name__ == "__main__":
    isAlignmentPossible, firstSequence, secondSequence = getFastaSequences()

    if isAlignmentPossible:
        matchScore = getMatchScore()
        gapScore = getGapScore()
        mismatchScore = getMismatchScore()
        alignmentMatrix = createArray(len(firstSequence), len(secondSequence))
        alignmentMatrix = fillArray(alignmentMatrix, firstSequence, secondSequence, matchScore, mismatchScore, gapScore)
        getAlignmentFromMatrix(alignmentMatrix, firstSequence, secondSequence, gapScore)