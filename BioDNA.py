import collections
import re
from PIL import Image
import time
from ansicol import color
import sys

def inputExplan(explain, lineLength):
# param explain: text file of explanations
# param lineLength: length of each line for explanation
# type explain: opened text file in ereading mode
# type lineLength: integer
# return explainDict: all explanations
# rtype explainDict: dictionary

    # puts explanations into list (if not \n)
    explainList = []
    for line in explain:
        if line is not "\n":
            explainList.append(line)

    # split lines according to characters
    find = "." * lineLength
    explainDict = collections.OrderedDict()
    for item in explainList:
        x = explainList.index(item)

        # split into list (with replacement field with x in name)
        name = "explain{}".format(x)
        explainDict[name] = re.findall(find, item)

        remainder = len(item) % lineLength
        length = len(item)
        if remainder != 0:
            explainDict[name].append(item[length-remainder:])

    # put whole words into end of previous item in list
    for listDict in explainDict.values():
        for itemDict in listDict:
            itemIndex = listDict.index(itemDict)
            if itemIndex is not 0:
                addToItem = ""
                previous = listDict[itemIndex-1]
                while previous is not " " and itemDict[0] is not " ":
                    addToItem += itemDict[0]
                    itemDict = itemDict[1:]
                listDict[itemIndex-1] += addToItem
                listDict[itemIndex] = itemDict[1:]

    return explainDict


def printExplan(explainDict, num):
# param explainDict: all explanations
# param num: number of explanation to print
# type explainDict: dictionary
# type num: integer

    print("\n")
    for text in explainDict["explain{}".format(num)]:
        print(text)
    print("\n")


def printImage(img):
# param img: name of image file
# type img: string

    image = Image.open(img)
    image.show()

    time.sleep(5)

    # process = subprocess.Popen(["display", img])
    # process.kill()                                          # TODO: close image window


def inputDNA():
# return DNA: user input gene of interest
# rtype DNA: string

    DNA = input("Input gene of interest: ")
    # data validation and verification (must be A, T, C or G)
    restart = True
    while True:
        while restart:
            DNA = DNA.upper()
            for letter in DNA:
                if letter is "A" or letter is "T" or letter is "C" or letter is "G":
                    continue
                else:
                    DNA = input("Input invalid. Input gene of interest: ")
                    restart = False
                break
            break
        if restart is False:
            restart = True
            continue
        break

    return DNA


def makeAminoAcidTable():
# return AATable: RNA triples and corresponding amino acid abbreviations
# rtype AATable: dictionary

    AAtable = { "UUU":"F|Phe","UUC":"F|Phe","UUA":"L|Leu","UUG":"L|Leu","UCU":"S|Ser","UCC":"S|Ser",
                "UCA":"S|Ser","UCG":"S|Ser","UAU":"Y|Tyr","UAC":"Y|Tyr","UAA":"*|***","UAG":"*|***",
                "UGU":"C|Cys","UGC":"C|Cys","UGA":"*|***","UGG":"W|Trp","CUU":"L|Leu","CUC":"L|Leu",
                "CUA":"L|Leu","CUG":"L|Leu","CCU":"P|Pro","CCC":"P|Pro","CCA":"P|Pro","CCG":"P|Pro",
                "CAU":"H|His","CAC":"H|His","CAA":"Q|Gln","CAG":"Q|Gln","CGU":"R|Arg","CGC":"R|Arg",
                "CGA":"R|Arg","CGG":"R|Arg","AUU":"I|Ile","AUC":"I|Ile","AUA":"I|Ile","AUG":"M|Met",
                "ACU":"T|Thr","ACC":"T|Thr","ACA":"T|Thr","ACG":"T|Thr","AAU":"N|Asn","AAC":"N|Asn",
                "AAA":"K|Lys","AAG":"K|Lys","AGU":"S|Ser","AGC":"S|Ser","AGA":"R|Arg","AGG":"R|Arg",
                "GUU":"V|Val","GUC":"V|Val","GUA":"V|Val","GUG":"V|Val","GCU":"A|Ala","GCC":"A|Ala",
                "GCA":"A|Ala","GCG":"A|Ala","GAU":"D|Asp","GAC":"D|Asp","GAA":"E|Glu",
                "GAG":"E|Glu","GGU":"G|Gly","GGC":"G|Gly","GGA":"G|Gly","GGG":"G|Gly"}

    return AAtable


def complementDNA_to_DNA(DNA):
# param DNA: user input (gene of interest)
# type DNA: string
# return compStrandDNA: DNA template strand (matching gene of interest)
# rtype compStrandDNA: string

    compStrandDNAList = []

    # loop through DNA to create list of complement strand
    for letter in DNA:
        if letter == "A":
            compStrandDNAList.append("T")
        elif letter == "T":
            compStrandDNAList.append("A")
        elif letter == "C":
            compStrandDNAList.append("G")
        elif letter == "G":
            compStrandDNAList.append("C")

    compStrandDNA = "".join(compStrandDNAList)
    return compStrandDNA



def complementDNA_to_RNA(DNA):
# param DNA: user input (gene of interest)
# type DNA: String
# return compStrandmRNA: mRNA complement strand (matching gene of interest)
# rType compStrandRNA: string

    compStrandRNAList = list(DNA)

    while "T" in compStrandRNAList:
        index = compStrandRNAList.index("T")
        compStrandRNAList.remove("T")
        compStrandRNAList.insert(index, "U")


    compStrandRNA = "".join(compStrandRNAList)
    return compStrandRNA


def transcribe(DNA):
# param DNA: user input (gene of interest)
# type DNA: String
# return compStrandmRNA: DNA complement strand (matching gene of interest)
# rType compStrandRNA: string
# return compStrandmRNA: mRNA complement strand (matching gene of interest)
# rType compStrandRNA: string

    compStrandDNA = complementDNA_to_DNA(DNA)
    compStrandmRNA = complementDNA_to_RNA(DNA)

    return compStrandDNA, compStrandmRNA


def askSymbol():
# return symbol: user input of desired symbol
# rType symbol: integer

    symbol = input("Do you want 3- or 1-letter symbols for the amino acids? ")
    while True:
        try:
            symbol = int(symbol)
        except ValueError:
            symbol = input("Must enter number. Do you want 3- or 1-letter symbols? ")
            continue
        if symbol is not 1 and symbol is not 3:
            symbol = input("Must enter 1 or 3. Do you want 3- or 1-letter symbols? ")
            continue
        break

    return symbol


def translate(mRNA, AAtable, symbol):
# param mRNA: mRNA complement to DNA
# param AATable: RNA triples and corresponding amino acid abbreviations
# symbol: 1- or 3-letter symbol
# type mRNA: String
# type AATable: dictionary
# type symbol: integer
# return protein: amino acids of gene of interest
# return metStart: Met at the start of gene
# return metIn: Met in gene
# return stopEnd: stop codon at the end of gene
# return stopIn: stop codon in gene
# rType protein: string
# rType metStart: boolean
# rType metIn: boolean
# rType stopEnd: boolean
# rType stopIn: boolean

    startBP = 0
    endBP = len(mRNA)
    nextCodonStart = startBP  # the index of the first nucleotide in the codon to be read, increments by 3 each time
    nextCodonEnd = nextCodonStart + 3

    protein = ""
    AASymbol = ""

    # process of translation -- until there are less than 3 base pairs to translate
    while nextCodonEnd <= endBP:
        nextCodon = mRNA[nextCodonStart:nextCodonEnd]       # the three nucleotides to be read
        nextAA = AAtable[nextCodon]                         # amino acid produced by codon

        # parse symbols based on user input
        if symbol is 1:
            AASymbol = " " + nextAA[0] + " "                # amino acid symbol -- parsed from AATable dictionary
        elif symbol is 3:
            AASymbol = nextAA[2:]

        protein += AASymbol

        nextCodonStart += 3
        nextCodonEnd = nextCodonStart + 3

    # check for summary
    metStart = False
    metIn = False
    stopEnd = False
    stopIn = False

    if len(protein) >= 3:
        # check for Met
        if protein[:3] == "Met" or protein[0] == "M":
            metStart = True
        for index in range(0, len(protein), symbol):
            if index != 0:
                if protein[index:index+3] == "Met" or protein[index] == "M":
                    metIn = True
                    break
        # check for stop codon
        if protein[len(protein)-3:] == "***":
            stopEnd = True
        for index in range(0, len(protein)-symbol, symbol):
            if protein[index:index+3] == "***" or protein[index] == "*":
                stopIn = True
                break

    return protein, metStart, metIn, stopEnd, stopIn


def printReadingFrames(num, compStrandDNA, mRNA, protein, metStart, metIn, stopEnd, stopIn):
# param num: reading frame number
# param compStrandDNA: DNA template strand (matching gene of interest)
# param compStrandmRNA: mRNA complement strand (matching gene of interest)
# param protein: amino acids of gene of interest
# param metStart: Met at the start of gene
# param metIn: Met in gene
# param stopEnd: stop codon at the end of gene
# param stopIn: stop codon in gene
# type num: integer
# type compStrandDNA: string
# type mRNA: string
# type protein: string
# type metStart: boolean
# type metIn: boolean
# type stopEnd: boolean
# type stopIn: boolean

    read = ""

    # beginning of reading frame
    read += "Reading Frame #{}\n".format(num)
    read += "RNA length: %i\n" % len(compStrandDNA)

    # add color to mRNA
    mRNAList = list(mRNA)
    for letter in mRNAList:
        letterIndex = mRNAList.index(letter)
        if letter == "A":
            mRNAList[letterIndex] = color('cyan') + letter + color('reset')
        elif letter == "U":
            mRNAList[letterIndex] = color('blue') + letter + color('reset')
        elif letter == "C":
            mRNAList[letterIndex] = color('bold red') + letter + color('reset')
        elif letter == "G":
            mRNAList[letterIndex] = color('bold magenta') + letter + color('reset')
    mRNA = ''.join(mRNAList)

    # continued beginning
    read += "\tRNA:\t\t%s\n" % mRNA
    lines = "|  " * (len(compStrandDNA) // 3)
    read += "\t\t\t\t %s\n" % lines

    read += "\tProtein:\t%s\n\n" % protein

    # summary for Met at start
    if metStart == True:
        read += "There is a Met at the start site.\n"
    else:
        read += "There is no Met at the start site.\n"

    # summary for stop at end
    if stopEnd == True:
        read += "There is a stop at the end.\n"
    else:
        read += "There is no stop at the end.\n"

    # summary for Met inside strand
    if metIn == True:
        read += "There is at least one Met somewhere other than at the start.\n"
    else:
        read += "There are no Mets other than at the start.\n"

    # summary for stop inside strand
    if stopIn == True:
        read += "There is at least one stop interrupting the gene."
    else:
        read += "There are no interrupting stops."

    return read


def chooseFrame(read1, read2, read3):
# param read1: first reading frame
# param read2: second reading frame
# param read3: third reading frame
# type read1: string
# type read2: string
# type read3: string

    # ask user for input to choose best reading frame
    print("\n")
    choose = input("Which is the best reading frame? ")
    while True:
        try:
            choose = int(choose)
        except ValueError:
            choose = input("Must enter number. Which is the best reading frame? ")
            continue
        if choose is not 1 and choose is not 2 and choose is not 3:
            choose = input("Must enter 1, 2 or 3. Which is the best reading frame? ")
        break

    # save to text file
    file = input("Enter the desired file name: ") + ".txt"

    fileInsert = input("Do you want to write over (\"w\") or append (\"a\") to the file? ")
    fileInsert.lower()
    while True:
        if fileInsert is not "w" and fileInsert is not "a":
            fileInsert = input("Invalid input. Do you want to write over (\"w\") or append (\"a\") to the file? ")
        else:
            break

    fileOpen = open(file, fileInsert)

    findColor = ["[0m", "[34m", "[36m", "[1;35m", "[1;31m"]
    if choose == 1:
        for colorCode in findColor:
            while colorCode in read1:
                read1 = read1.replace(colorCode, "")
        if fileInsert == "a":
            read1 = "\n\n" + read1
        fileOpen.write(read1)
    elif choose == 2:
        for colorCode in findColor:
            while colorCode in read2:
                read2 = read2.replace(colorCode, "")
        if fileInsert == "a":
            read2 = "\n\n" + read2
        fileOpen.write(read2)
    elif choose == 3:
        for colorCode in findColor:
            while colorCode in read3:
                read3 = read3.replace(colorCode, "")
        if fileInsert == "a":
            read3 = "\n\n" + read3
        fileOpen.write(read3)

    fileOpen.close()


def programCont():
    # ask user for input to continue or quit program
    next = input("Do you want to continue or quit? ")
    while True:
        next = next.lower()
        if next == 'c':
            print("\n\n")
            break
        elif next == "q":
            sys.exit()
        else:
            next = input("Invalid input. Do you want to continue or quit? ")