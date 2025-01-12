#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import sys
import textwrap
import argparse


def ID(seq1, seq2):
    counter = 0
    for i in range(len(seq1)):
        if "-" not in [seq1[i], seq2[i]]:
            if seq1[i] == seq2[i]:
                counter += 1
    return (counter/len(seq1))*100


def RemoveLeadingSpaces(line):
    counter = 0
    newLine = ''
    for i in line:
        if i == " ":
            if counter == 0:
                pass
            else:
                newLine += i
        else:
            counter += 1
            newLine += i
    return newLine


def Unique(ls):
    unqList = []
    for i in ls:
        if i not in unqList:
            unqList.append(i)
    return unqList


def reverseComplement(seq):
    out = []
    for i in range(len(seq)-1, -1, -1):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def Complement(seq):
    out = []
    for i in range(0, len(seq)):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def ribosome(seq):
    NTs = ['T', 'C', 'A', 'G']
    stopCodons = ['TAA', 'TAG', 'TGA']
    Codons = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codon = NTs[i] + NTs[j] + NTs[k]
                # if not codon in stopCodons:
                Codons.append(codon)

    CodonTable = {}
    AAz = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    AAs = list(AAz)
    k = 0
    for base1 in NTs:
        for base2 in NTs:
            for base3 in NTs:
                codon = base1 + base2 + base3
                CodonTable[codon] = AAs[k]
                k += 1

    prot = []
    for j in range(0, len(seq), 3):
        codon = seq[j:j + 3]
        try:
            prot.append(CodonTable[codon])
        except KeyError:
            prot.append("")
    protein = ("".join(prot))
    return protein


def SeqCoord(seq, start, end):
    return seq[start:end]


def howMany(ls, exclude):
    counter = 0
    for i in ls:
        if i != exclude:
            counter += 1
    return counter


def stabilityCounter5(int):
    if len(str(int)) == 1:
        string = (str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 2:
        string = (str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 3:
        string = (str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 4:
        string = (str(0) + str(int))
        return (string)
    if len(str(int)) > 4:
        string = str(int)
        return (string)


def stabilityCounter6(int):
    if len(str(int)) == 1:
        string = (str(0) + str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 2:
        string = (str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 3:
        string = (str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 4:
        string = str(0) + (str(0) + str(int))
        return (string)
    if len(str(int)) == 5:
        string = (str(0) + str(int))
        return (string)
    if len(str(int)) > 5:
        string = str(int)
        return (string)


def sum(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count


def ave(ls):
    count = 0
    for i in ls:
        try:
            count += float(i)
        except ValueError:
            pass
    return count/len(ls)


def derep(ls):
    outLS = []
    for i in ls:
        if i not in outLS:
            outLS.append(i)
    return outLS


def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        #->>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        #->>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    # data = sorted(data)
    data.sort(key=int)
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


def GCcalc(seq):
    count = 0
    for i in seq:
        if i == "G" or i == "C":
            count += 1
    return count/len(seq)


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
            x = i
    return x


def RemoveDuplicates(ls):
    empLS = []
    counter = 0
    for i in ls:
        if i not in empLS:
            empLS.append(i)
        else:
            pass
    return empLS


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length-1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)-1]


def secondToLastItem(ls):
    x = ''
    for i in ls[0:len(ls)-1]:
        x = i
    return x


def pull(item, one, two):
    ls = []
    counter = 0
    for i in item:
        if counter == 0:
            if i != one:
                pass
            else:
                counter += 1
                ls.append(i)
        else:
            if i != two:
                ls.append(i)
            else:
                ls.append(i)
                counter = 0
    outstr = "".join(ls)
    return outstr


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


def replaceLS(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    return emptyList


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def removeLS(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    return emptyList


def fasta(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def fasta2(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def allButTheFirst(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(1, length):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)-1]


def filter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def filterRe(list, regex):
    ls1 = []
    ls2 = []
    for i in list:
        if re.findall(regex, i):
            ls1.append(i)
        else:
            ls2.append(i)
    return ls1, ls2


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls = filter(ls, [""])
    return ls


parser = argparse.ArgumentParser(
    prog="breseq_parser.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************

    Developed by Arkadiy Garber; University of Montana, Biological Sciences
    Please send comments and inquiries to arkadiy.garber@mso.umt.edu
    ************************************************************************
    '''))

parser.add_argument('-b', type=str, help="breseq_output_directory", default="NA")

parser.add_argument('-gbk', type=str, help="input genbank file", default="NA")

parser.add_argument('-gff', type=str, help="input gff file", default="")

parser.add_argument('-o', type=str, help="outdir", default="NA")

parser.add_argument('-outdir', type=str, help="outdir as defined by the vc_looper.sh", default="NA")

parser.add_argument('-m', type=str, help="mode for breseq run: clone (default) or pop", default="clone")


if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]

os.system("mkdir -p %s" % args.o)

gbfDict = defaultdict(lambda: '-')
oldLocus = ""
locus = ""
if args.gff != "":
    gff = open(args.gff)
else:
    gbk = open(args.gbk)
    for i in gbk:
        if re.findall(r'/locus_tag', i):
            locus = (remove(i.rstrip(), [" ", "\""]).split("=")[1])
            oldLocus = ""
        if re.findall(r'old_locus_tag', i):
            oldLocus = (remove(i.rstrip(), [" ", "\""]).split("=")[1])

        if len(oldLocus) > 0:
            if len(locus) > 0:
                gbfDict[locus] = oldLocus
        else:
            if len(locus) > 0:
                gbfDict[locus] = locus

gffDict = defaultdict(lambda: defaultdict(lambda: 'intergenic'))
gffDict2 = defaultdict(lambda: defaultdict(lambda: 'intergenic'))
gff = open(args.b + "/data/reference.gff3")
for i in gff:
    if re.match('##FASTA', i):
        break
    else:
        if not re.match(r'#', i):
            ls = i.rstrip().split("\t")
            try:
                alias = ls[8].split(";")[0].split("=")[1]

                if args.gff == "":
                    if ls[2] == "CDS":
                        name = (ls[8].split("Name=")[1].split(";")[0])
                        if name != alias:
                            gene = name
                        else:
                            gene = name

                        product = (ls[8].split("Note=")[1].split(";")[0])

                        gffDict2[alias]["gene"] = gene
                        gffDict2[alias]["product"] = product

                        for j in range(int(ls[3]), int(ls[4])+1):
                            gffDict[ls[0]][j] = alias
                else:
                    if ls[2] == "CDS":
                        gene = "-"

                        product = "-"

                        gffDict2[alias]["gene"] = gene
                        gffDict2[alias]["product"] = product

                        for j in range(int(ls[3]), int(ls[4]) + 1):
                            gffDict[ls[0]][j] = alias
            except IndexError:
                pass


gDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
gd = open(args.b + "/data/annotated.gd")
for i in gd:
    if not re.match(r'#', i):
        ls = i.rstrip().split("\t")
        if len(ls[0]) == 3:
            if args.m == "polymorphism-prediction":
                freq = float(i.rstrip().split("frequency=")[1].split("\t")[0]) * 100
            else:
                freq = 'NA'
            coord = ls[3] + "|" + ls[4]
            gDict[coord]["mutation"] = replaceLS(ls[0:10], [","], ";")
            gDict[coord]["freq"] = freq
        elif len(ls[0]) == 2:
            if ls[0] != "UN":
                coord = ls[3] + "|" + ls[4]
                if coord not in gDict.keys() and ls[0] != "UN":
                    coord = ls[3] + "|" + ls[7]
                gDict[coord]["evidence"] = replaceLS(ls[0:10], [","], ";")

eviDict = {"RA": "read alignment", "MC": "missing coverage", "JC": "new junction", "UN": "unknown base"}

outfile = args.b
out = open(args.o + "/" + lastItem(outfile.split("/")) + ".csv", "w")
out.write("identifier,sequence,position,mutation_type,mutation,seq_change,evidence,locus,old_locus,gene/locus,product,freq\n")
for i in gDict.keys():
    if len(gDict[i]) > 1:
        mutLS = gDict[i]["mutation"]
        eviLS = gDict[i]["evidence"]
        freq = gDict[i]["freq"]
        evidence = eviLS[0]
        if evidence in eviDict.keys():
            evidenceLong = eviDict[evidence]
        else:
            evidenceLong = 'unknown evidence'
        mutType = mutLS[0]
        contig = mutLS[3]
        position = mutLS[4]
        locus = (gffDict[contig][int(position)])
        oldLocus = gbfDict[locus]
        if mutLS[0] == "SNP" and eviLS[6] == "N":
            pass

        else:
            if mutType == "SNP":
                mutTypeLong = "single-nucleotide polymorphism"
                product = replace(gffDict2[locus]["product"], [","], ";")
                mutation = eviLS[6] + "->" + eviLS[7]
                AAposition = mutLS[7].split("=")[1]
                if locus != "intergenic":
                    gene = gffDict2[locus]["gene"]
                    newAA = mutLS[6].split("=")[1]
                    refAA = mutLS[8].split("=")[1]
                    change = refAA + AAposition + newAA

                    out.write(args.b + "," + contig + "," + str(position) + "," +
                              str(mutTypeLong + " (" + mutType + ")") + "," + str(mutation) + "," +
                              str(change) + "," + str(evidenceLong + " (" + evidence + ")") + "," +
                              locus + "," + oldLocus + "," + str(gene) + "," + str(product) + "," + str(freq) + "\n")
                else:
                    geneAlias = mutLS[6].split("=")[1]
                    change = mutLS[7].split("=")[1]

                    out.write(args.b + "," + contig + "," + str(position) + "," +
                              str(mutTypeLong + " (" + mutType + ")") + "," + str(mutation) + "," +
                              str(change) + "," + str(evidenceLong + " (" + evidence + ")") + "," +
                              locus + "," + oldLocus + "," + str(geneAlias) + "," + str(product) + "," + str(freq) + "\n")

            elif mutType == "INS":
                mutTypeLong = "insertion"
                product = replace(mutLS[8].split("=")[1], [","], ";")
                if locus != "intergenic":
                    gene = mutLS[6].split("=")[1]
                    mutation = "->" + mutLS[5]
                    change = mutLS[7].split("=")[1]

                    out.write(args.b + "," + contig + "," + str(position) + "," +
                              str(mutTypeLong + " (" + mutType + ")") + "," + str(mutation) + "," +
                              str(change) + "," + str(evidenceLong + " (" + evidence + ")")
                              + "," + locus + "," + oldLocus + "," + str(gene) + "," + str(product) + "," + str(freq) + "\n")
                else:
                    geneAlias = mutLS[6].split("=")[1]
                    change = mutLS[7].split("=")[1]
                    mutation = "->" + mutLS[5]

                    out.write(args.b + "," + contig + "," + str(position) + "," +
                              str(mutTypeLong + " (" + mutType + ")") + "," + str(mutation) + "," +
                              str(change) + "," + str(evidenceLong + " (" + evidence + ")")
                              + "," + locus + "," + oldLocus + "," + str(geneAlias) + "," + str(product) + "," + str(freq) + "\n")

            elif mutType == "DEL":
                mutTypeLong = "deletion"
                product = replace(mutLS[8].split("=")[1], [","], ";")
                gene = mutLS[6].split("=")[1]
                if locus != "intergenic":
                    mutation = "-" + mutLS[5] + " bp"
                    change = mutLS[7].split("=")[1]

                    out.write(args.b + "," + contig + "," + str(position) + "," +
                              str(mutTypeLong + " (" + mutType + ")") + "," + str(mutation) + "," +
                              str(change) + "," + str(evidenceLong + " (" + evidence + ")")
                              + "," + locus + "," + oldLocus + "," + str(gene) + "," + str(product) + "," + str(freq) + "\n")
                else:
                    geneAlias = mutLS[6].split("=")[1]
                    mutation = "-" + mutLS[5] + " bp"
                    change = mutLS[7].split("=")[1]

                    out.write(args.b + "," + contig + "," + str(position) + "," +
                              str(mutTypeLong + " (" + mutType + ")") + "," + str(mutation) + "," +
                              str(change) + "," + str(evidenceLong + " (" + evidence + ")")
                              + "," + locus + "," + oldLocus + "," + str(geneAlias) + "," + str(product) + "," + str(freq) + "\n")

            elif mutType == "SUB":
                mutTypeLong = "multiple base substitution"
                if locus != "intergenic":
                    geneAlias = mutLS[7].split("=")[1]
                    mutation = "-" + mutLS[5] + " bp" + "(" + "->+" + mutLS[6] + ")"
                    change = "(" + "->+" + mutLS[6] + ")"
                    product = replace(mutLS[8].split("=")[1], [","], ";")
                    mutLS = gDict[i]["mutation"]
                    eviLS = gDict[i]["evidence"]
                    evidence = eviLS[0]
                    if evidence in eviDict.keys():
                        evidenceLong = eviDict[evidence]
                    else:
                        evidenceLong = 'unknown evidence'
                    mutType = mutLS[0]
                    contig = mutLS[3]
                    position = mutLS[4]
                    locus = (gffDict[contig][int(position)])
                    oldLocus = gbfDict[locus]
                    # print(mutLS)
                    # print(eviLS)
                    # print("")

                    out.write(args.b + "," + contig + "," + str(position) + "," +
                              str(mutTypeLong + " (" + mutType + ")") + "," + str(mutation) + "," +
                              str(change) + "," + str(evidenceLong + " (" + evidence + ")")
                              + "," + locus + "," + oldLocus + "," + str(geneAlias) + "," + str(product) + "," + str(freq) + "\n")

                    # print(args.b + "," + contig + "," + str(position) + "," +
                    #           str(mutTypeLong + " (" + mutType + ")") + "," + str(mutation) + "," +
                    #           str(change) + "," + str(evidenceLong + " (" + evidence + ")")
                    #           + "," + locus + "," + oldLocus + "," + str(geneAlias) + "," + str(product) + "\n")
                else:
                    mutation = "-" + mutLS[5] + " bp" + "(" + mutLS[6] + ")"
                    change = mutLS[8]

                    geneAlias = mutLS[8].split("=")[1]

                    out.write(args.b + "," + contig + "," + str(position) + "," +
                              str(mutTypeLong + " (" + mutType + ")") + "," + str(mutation) + "," +
                               str(change) + "," + str(evidenceLong + " (" + evidence + ")")
                              + "," + locus + "," + oldLocus + "," + str(geneAlias) + ",-" + "," + str(freq) + "\n")
            else:
                pass
                # print("mutation" + "\t" + str(mutLS))
                # print("evidence" + "\t" + str(eviLS))
                # print("")

out.close()