import BioDNA
import time

def main():

    print ("\n ---+++ IngeneT v5.136 +++---")

    # opens text file containing explanations
    explain = open('explain.txt', 'r')
    explainDict = BioDNA.inputExplan(explain, 150)

    # print introduction
    BioDNA.printExplan(explainDict, 0)

    times = 0
    while True:
        # asks for user input (gene desired)
        DNA = BioDNA.inputDNA()

        # prepare DNA variables for three reading frames
        DNA1 = DNA
        DNA2 = DNA[1:]
        DNA3 = DNA[2:]

        # print transcription explanation
        if times == 0:
            BioDNA.printExplan(explainDict, 1)
            BioDNA.printImage("transcription.jpg")

        # transcription process
        compStrandDNA1, mRNA1 = BioDNA.transcribe(DNA1)
        compStrandDNA2, mRNA2 = BioDNA.transcribe(DNA2)
        compStrandDNA3, mRNA3 = BioDNA.transcribe(DNA3)

        # print translation explanation
        if times == 0:
            BioDNA.printExplan(explainDict, 2)
            BioDNA.printImage("translation.jpg")
            # TODO: add press enter to continue
            BioDNA.printExplan(explainDict, 3)
            BioDNA.printImage("table.png")

        # translation process
        # create dictionary of amino acids
        AATable = BioDNA.makeAminoAcidTable()
        # ask for 1- or 3- letter symbols
        symbol = BioDNA.askSymbol()
        # translation
        protein1, metStart1, metIn1, stopEnd1, stopIn1 = BioDNA.translate(mRNA1, AATable, symbol)
        protein2, metStart2, metIn2, stopEnd2, stopIn2 = BioDNA.translate(mRNA2, AATable, symbol)
        protein3, metStart3, metIn3, stopEnd3, stopIn3 = BioDNA.translate(mRNA3, AATable, symbol)

        # print reading frames explanation
        if times == 0:
            BioDNA.printExplan(explainDict, 4)

        # print reading frames with proper formatting
        # first part (original DNA)
        print("==============================================================\n")
        print("   DNA: 5' %s 3'\t gene of interest to be transcribed" % DNA)
        print("   DNA: 3' %s 5'\t the template strand" % compStrandDNA1)
        lines = "|" * len(DNA)
        print("\t\t   %s" % lines)
        print("  mRNA: 5' %s 3'" % mRNA1)
        # second part (three reading frames)
        print("\n==============================================================")
        time.sleep(1)
        read1 = BioDNA.printReadingFrames(1, compStrandDNA1, mRNA1, protein1, metStart1, metIn1, stopEnd1, stopIn1)
        print(read1)
        print("======================================")
        time.sleep(2)
        print("\n======================================")
        read2 = BioDNA.printReadingFrames(2, compStrandDNA2, mRNA2, protein2, metStart2, metIn2, stopEnd2, stopIn2)
        print(read2)
        print("======================================")
        time.sleep(2)
        print("\n======================================")
        read3 = BioDNA.printReadingFrames(3, compStrandDNA3, mRNA3, protein3, metStart3, metIn3, stopEnd3, stopIn3)
        print(read3)
        print("==============================================================")

        # ask user to choose best reading frame
        BioDNA.chooseFrame(read1, read2, read3)

        # print conclusion
        if times == 0:
            BioDNA.printExplan(explainDict, 5)

        # ask user to continue or quit program
        BioDNA.programCont()
        times += 1


main()