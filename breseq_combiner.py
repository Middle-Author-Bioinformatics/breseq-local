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
    prog="breseq_combiner.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************

    Developed by Arkadiy Garber; University of Montana, Biological Sciences
    Please send comments and inquiries to arkadiy.garber@mso.umt.edu
    ************************************************************************
    '''))

parser.add_argument('-i', type=str, help="input dir", default="NA")

parser.add_argument('-m', type=str, help="mode for breseq run: clone (default) or pop", default="clone")


if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]

if args.m != "polymorphism-prediction":
    locusDict = defaultdict(lambda: defaultdict(lambda: '-'))
    counter = 0
    out = open(args.i + "/detailed_summary.csv", "w")
    summaryDict = defaultdict(lambda: defaultdict(lambda: '-'))
    summaries = os.listdir(args.i)
    for i in summaries:
        if i not in ["mutations_profile.csv", "detailed_summary.csv", "mutations_profile-melt.csv", "mutations_profile-filtered.csv", "mutations_profile.pdf"]:
            if re.findall(r'csv', i):
                file = open("%s/%s" % (args.i, i))
                for j in file:
                    ls = j.rstrip().split(",")
                    if ls[7] != "locus":
                        locus = remove(ls[7], ["[", "]"])
                        oldlocus = remove(ls[8], ["[", "]"])
                        if "missing coverage (MC)" not in [locus, oldlocus]:
                            locusDict[locus]["gene"] = remove(ls[9], ["[", "]"])
                            locusDict[locus]["product"] = remove(ls[10], ["[", "]"])
                            locusDict[locus]["oldlocus"] = oldlocus
                            out.write(j.rstrip() + "\n")

                            if locus not in ["intergenic"]:
                                summaryDict[locus][i.split(".")[0]] = ls[4].split(" (")[0]
                    else:
                        if counter == 0:
                            out.write(j.rstrip() + "\n")
                            counter += 1
    out.close()

    out = open(args.i + "/mutations_profile.csv", "w")
    out.write("locus,old_locus,gene,product")
    for i in summaries:
        if i not in ["mutations_profile.csv", "detailed_summary.csv",
                     "mutations_profile-melt.csv", "mutations_profile-filtered.csv", "mutations_profile.pdf"]:
            out.write("," + lastItem(i.split(".")[0].split("_")))
    out.write("\n")

    for i in summaryDict.keys():
        out.write(i + "," + str(locusDict[i]["oldlocus"]) + "," + str(locusDict[i]["gene"]) + "," + str(locusDict[i]["product"]))
        for j in summaries:
            if j not in ["mutations_profile.csv", "detailed_summary.csv", "mutations_profile-melt.csv",
                         "mutations_profile-filtered.csv", "mutations_profile.tiff"]:
                out.write("," + str(summaryDict[i][j.split(".")[0]]))
        out.write("\n")
    out.close()

else:
    locusDict = defaultdict(lambda: defaultdict(lambda: '-'))
    counter = 0
    out = open(args.i + "/detailed_summary.csv", "w")
    summaryDict = defaultdict(lambda: defaultdict(lambda: '-'))
    summaries = os.listdir(args.i)
    for i in summaries:
        if i not in ["mutations_profile.csv", "detailed_summary.csv", "mutations_profile-melt.csv",
                     "mutations_profile-filtered.csv", "mutations_profile.pdf"]:
            if re.findall(r'csv', i):
                file = open("%s/%s" % (args.i, i))
                for j in file:
                    ls = j.rstrip().split(",")
                    if ls[7] != "locus":
                        locus = remove(ls[7], ["[", "]"])
                        oldlocus = remove(ls[8], ["[", "]"])
                        freq = ls[11]
                        if "missing coverage (MC)" not in [locus, oldlocus]:
                            locusDict[locus]["gene"] = remove(ls[9], ["[", "]"])
                            locusDict[locus]["product"] = remove(ls[10], ["[", "]"])
                            locusDict[locus]["oldlocus"] = oldlocus
                            out.write(j.rstrip() + "\n")

                            if locus not in ["intergenic"]:
                                summaryDict[locus][i.split(".")[0]] = freq
                    else:
                        if counter == 0:
                            out.write(j.rstrip() + "\n")
                            counter += 1
    out.close()

    out = open(args.i + "/mutations_profile.csv", "w")
    out.write("locus,old_locus,gene,product")
    for i in summaries:
        if i not in ["mutations_profile.csv", "detailed_summary.csv",
                     "mutations_profile-melt.csv", "mutations_profile-filtered.csv", "mutations_profile.pdf"]:
            out.write("," + lastItem(i.split(".")[0].split("_")))
    out.write("\n")

    for i in summaryDict.keys():
        try:
            gene = (locusDict[i]["gene"])
            num = float(locusDict[i]["gene"])
            out.write(i + "," + str(locusDict[i]["oldlocus"]) + "," + str("-") + "," + str(locusDict[i]["product"]))
        except ValueError:
            out.write(i + "," + str(locusDict[i]["oldlocus"]) + "," + str(locusDict[i]["gene"]) + "," + str(locusDict[i]["product"]))
        for j in summaries:
            if j not in ["mutations_profile.csv", "detailed_summary.csv", "mutations_profile-melt.csv",
                         "mutations_profile-filtered.csv", "mutations_profile.tiff"]:
                out.write("," + str(summaryDict[i][j.split(".")[0]]))
        out.write("\n")
    out.close()



# UNDER CONSTRUCTION
#
# samples = []
# geneDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
# profileDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
# mutDict = defaultdict(list)
# profile = open(args.o + "/mutations_profile.csv")
# for i in profile:
#     ls = i.rstrip().split(",")
#     if ls[0] != "locus":
#         for j in range(4, len(samples)):
#             print(j)
#             print("+")
#             print(samples[j])
#             print(ls[j])
#             profileDict[ls[0]][int(samples[j])] = ls[j]
#             mutDict[ls[0]].append(ls[j])
#             geneDict[ls[0]]["gene"] = ls[1]
#             geneDict[ls[0]]["product"] = ls[2]
#         # print("")
#
#     else:
#         samples = ls
#         # print(samples)
#         header = i.rstrip()
#
# '''
# mutationsToRemove = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
# genesToRemove = []
# for i in profileDict.keys():
#     dupDict = defaultdict(list)
#     # print(i)
#     mutations = (mutDict[i])
#     # print(mutations)
#     sampleSize = len(mutations)
#     for j in mutations:
#         if j != "-":
#             dupDict[j].append(i)
#     counter = 0
#     for j in dupDict.keys():
#         prop = len(dupDict[j]) / sampleSize
#         if prop < 1:
#             # print(j + '\t\t' + str( prop ))
#             counter += 1
#         else:
#             mutationsToRemove[i][j] = "remove"
#     if counter == 0:
#         genesToRemove.append(i)
#     # print("")
#
# for i in genesToRemove:
#     profileDict.pop(i)
# '''
#
# out2 = open(args.o + "/mutations_profile-melt.csv", "w")
# out2.write("locus,sample,mutation\n")
# for i in sorted(profileDict.keys(), reverse=True):
#
#     for j in sorted(profileDict[i]):
#         mutation = (profileDict[i][j])
#         if mutation == "-":
#             TYPE = "none"
#         elif re.match(r'->', mutation):
#             TYPE = "insertion"
#         elif re.match(r'-', mutation):
#             TYPE = "deletion"
#         else:
#             ls = (mutation.split("->"))
#             if ls[0] in ["A", "C"] and ls[1] in ["A", "C"]:
#                 TYPE = "transversion"
#             if ls[0] in ["T", "G"] and ls[1] in ["T", "G"]:
#                 TYPE = "transversion"
#
#             if ls[0] in ["A", "T"] and ls[1] in ["A", "T"]:
#                 TYPE = "transversion"
#             if ls[0] in ["C", "G"] and ls[1] in ["C", "G"]:
#                 TYPE = "transversion"
#
#             if ls[0] in ["A", "G"] and ls[1] in ["A", "G"]:
#                 TYPE = "transition"
#             if ls[0] in ["T", "C"] and ls[1] in ["T", "C"]:
#                 TYPE = "transition"
#
#         out2.write(str(i) + "," + str(j) + "," + TYPE + "\n")
#
#
# out2.close()
# # out.close()
#
#














