import re, math
from collections import Counter
import csv
import numpy as np
import Bio
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd


def CountCTDC(seq1, seq2):
    sum = 0
    for aa in seq1:
        sum = sum + seq2.count(aa)
    return sum


def CountCTDD(aaSet, sequence):
	number = 0
	for aa in sequence:
		if aa in aaSet:
			number = number + 1
	cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
	cutoffNums = [i if i >=1 else 1 for i in cutoffNums]

	code = []
	for cutoff in cutoffNums:
		myCount = 0
		for i in range(len(sequence)):
			if sequence[i] in aaSet:
				myCount += 1
				if myCount == cutoff:
					code.append((i + 1) / len(sequence) * 100)
					break
		if myCount == 0:
			code.append(0)
	return code


def CalculateCTriad(sequence, gap, features, AADict):
    res = []
    for g in range(gap+1):
        myDict = {}
        for f in features:
            myDict[f] = 0
        for i in range(len(sequence)):
            if i+gap+1 < len(sequence) and i+2*gap+2<len(sequence):
                fea = AADict[sequence[i]] + '.' + AADict[sequence[i+gap+1]]+'.'+AADict[sequence[i+2*gap+2]]
                myDict[fea] = myDict[fea] + 1
        maxValue, minValue = max(myDict.values()), min(myDict.values())
        for f in features:
            res.append((myDict[f] - minValue) / maxValue)
    return res

def Rvalue(aa1, aa2, AADict, Matrix):
    return sum([(Matrix[i][AADict[aa1]] - Matrix[i][AADict[aa2]]) ** 2 for i in range(len(Matrix))]) / len(Matrix)

def count_motif(seq, motif):
    motif = motif.replace("-x",".").replace("-","").replace("{","[^").replace("}","]").replace("(","{").replace(")","}").replace("<","^").replace(">","$")
    return len(re.findall(motif,seq))



AA = 'ACDEFGHIKLMNPQRSTVWY'
lambdaValue = 29
w=0.05
nlag = 29

encodings = []
header = ['Entry','Entry name','Sequence']
#Constructing AAC, DPC and TPC Headers
diPeptides = [aa1 + aa2 for aa1 in AA for aa2 in AA]
triPeptides = [aa1 + aa2 + aa3 for aa1 in AA for aa2 in AA for aa3 in AA]
for i in AA:
    header.append(i)
header = header + diPeptides + triPeptides

#Constructing GAAC header
gaac_group = {
		'aliphatic': 'GAVLMI',
		'aromatic': 'FYW',
		'postivecharge': 'KRH',
		'negativecharge': 'DE',
		'uncharge': 'STCPNQ'
	}

gaac_groupKey = gaac_group.keys()
for key in gaac_groupKey:
    header.append(key)

#Constructing Moran Header
props=['CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102', 'CHOC760101', 'BIGC670101', 'CHAM810101', 'DAYM780201']
moran_head = []

for p in props:
    for n in range(1, nlag+1):
        moran_head.append('Moran_' + p + '.lag' + str(n))

header = header + moran_head


#Constructing CTD Header
group1 = {
		'hydrophobicity_PRAM900101': 'RKEDQN',
		'hydrophobicity_ARGP820101': 'QSTNGDE',
		'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
		'hydrophobicity_PONP930101': 'KPDESNQT',
		'hydrophobicity_CASG920101': 'KDEQPSRNTG',
		'hydrophobicity_ENGD860101': 'RDKENQHYP',
		'hydrophobicity_FASG890101': 'KERSQD',
		'normwaalsvolume': 'GASTPDC',
		'polarity':        'LIFWCMVY',
		'polarizability':  'GASDT',
		'charge':          'KR',
		'secondarystruct': 'EALMQKRH',
		'solventaccess':   'ALFCGIVW'
	}
group2 = {
		'hydrophobicity_PRAM900101': 'GASTPHY',
		'hydrophobicity_ARGP820101': 'RAHCKMV',
		'hydrophobicity_ZIMJ680101': 'HMCKV',
		'hydrophobicity_PONP930101': 'GRHA',
		'hydrophobicity_CASG920101': 'AHYMLV',
		'hydrophobicity_ENGD860101': 'SGTAW',
		'hydrophobicity_FASG890101': 'NTPG',
		'normwaalsvolume': 'NVEQIL',
		'polarity':        'PATGS',
		'polarizability':  'CPNVEQIL',
		'charge':          'ANCQGHILMFPSTWYV',
		'secondarystruct': 'VIYCWFT',
		'solventaccess':   'RKQEND'
        }
group3 = {
        	'hydrophobicity_PRAM900101': 'CLVIMFW',
		'hydrophobicity_ARGP820101': 'LYPFIW',
		'hydrophobicity_ZIMJ680101': 'LPFYI',
		'hydrophobicity_PONP930101': 'YMFWLCVI',
		'hydrophobicity_CASG920101': 'FIWC',
		'hydrophobicity_ENGD860101': 'CVLIMF',
		'hydrophobicity_FASG890101': 'AYHWVMFLIC',
		'normwaalsvolume': 'MHKFRYW',
		'polarity':        'HQRKNED',
		'polarizability':  'KMHFRYW',
		'charge':          'DE',
		'secondarystruct': 'GNPSD',
		'solventaccess':   'MSPTHY'
        }

groups = [group1, group2, group3]
property = (
	'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')


for p in property:
    for g in range(1, len(groups) + 1):
        header.append(p + '.G' + str(g))

for p in property:
    for g in ('1', '2', '3'):
        for d in ['0', '25', '50', '75', '100']:
            header.append(p + '.' + g + '.residue' + d)

for p in property:
    for tr in ('Tr1221', 'Tr1331', 'Tr2332'):
        header.append(p + '.' + tr)


#Constructing CTriad Header
AAGroup = {
	'g1': 'AGV',
	'g2': 'ILFP',
	'g3': 'YMTS',
	'g4': 'HNQW',
	'g5': 'RK',
	'g6': 'DE',
	'g7': 'C'
}

myGroups = sorted(AAGroup.keys())
AADict = {}
for g in myGroups:
    for aa in AAGroup[g]:
        AADict[aa] = g
features = [f1 + '.'+ f2 + '.' + f3 for f1 in myGroups for f2 in myGroups for f3 in myGroups]
for f in features:
    header.append(f)


#PAAC Header
with open('./data/PAAC.txt') as f:
    records = f.readlines()

AA2 = ''.join(records[0].rstrip().split()[1:])
AADict2 = {}
for i in range(len(AA2)):
    AADict2[AA2[i]] = i
AAProperty = []
AAPropertyNames = []
for i in range(1, len(records)):
    array = records[i].rstrip().split() if records[i].rstrip() != '' else None
    AAProperty.append([float(j) for j in array[1:]])
    AAPropertyNames.append(array[0])

AAProperty1 = []
for i in AAProperty:
    meanI = sum(i) / 20
    fenmu = math.sqrt(sum([(j-meanI)**2 for j in i])/20)
    AAProperty1.append([(j-meanI)/fenmu for j in i])

for aa in AA2:
    header.append('Xc1.' + aa)
for n in range(1, lambdaValue + 1):
    header.append('Xc2.lambda' + str(n))



#Constructing Header for other features
header.append("mol_wt")
header.append("isoelectric_pt")
header.append("gravy")
header.append("instability_index")
struct = ['helix','sheet','turn']
for i in struct:
    header.append(i)
ext = ['cysteines','cystines']
for i in ext:
    header.append(i)

header.append("Length")
header.append("Transmembrane")
header.append("DNA Binding")
header.append("Metal Binding")
header.append("Nucleotide Binding")


#Constructing Subcellular Location
subcell_list = ['cytoplasm','membrane','cell wall','secreted',
                'periplasm','cell surface','cell envelope','chlorosome',
                'cellular thylakoid membrane','cellular cromatopore membrane',
                'single-pass membrane protein','multi-pass membrane protein','peripheral membrane protein']

for s in subcell_list:
    header.append(s)

#Constructing Motifs Header
motif_list = []
with open("./data/unique_motifs.csv","r") as f:
    reader = csv.reader(f)
    for row in reader:
        motif_list.append(row[0])
        header.append(row[0])

#Constructing GoData Header
GoData_list = []
with open("./data/GoData.csv", "r") as g:
    reader = csv.reader(g)
    for row in reader:
        GoData_list.append(row[0])
        header.append(row[0])

#fields = ['Sequence','Gene ontology (molecular function)','Transmembrane']
data = pd.read_csv("modifiedPro.csv")

#name of the output dataset file to be generated
with open("Final_Dataset.csv","w+",newline='') as feature_file:
    #writing header to the file
    writer = csv.writer(feature_file)
    writer.writerow(header)
    for i in range(0,len(data)):
        sequence = data['Sequence'][i]
        go_id = data['Gene ontology (molecular function)'][i]
        trans = data['Transmembrane'][i]
        dna = data['DNA binding'][i]
        metal = data['Metal binding'][i]
        ncl = data['Nucleotide binding'][i]
        sub = data['Subcellular location [CC]'][i]
        entry_name = data['Entry name'][i]
        entry = data['Entry'][i]
        length = len(sequence)
        #skipping sequences that have any other than 20 Amino Acid
        if 'B' in sequence or 'J' in sequence or 'O' in sequence or 'U' in sequence or 'X' in sequence or 'Z' in sequence:
            continue
        #skipping sequences that have length less than 30
        if length<30:
            continue
        code = []
        code = [entry, entry_name, sequence]
        #Calculating AAC
        count = Counter(sequence)
        for key in count:
            count[key] = count[key]/len(sequence)
        for aa in AA:
            #appending AAC
            code.append(count[aa])

        #Calculating DPC
        AADict_AC = {}
        for i in range(len(AA)):
            AADict_AC[AA[i]] = i
        tmpCode = [0] * 400
        for j in range(len(sequence) - 2 + 1):
            tmpCode[AADict_AC[sequence[j]] * 20 + AADict_AC[sequence[j+1]]] = tmpCode[AADict_AC[sequence[j]] * 20 + AADict_AC[sequence[j+1]]] +1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        #appending DPC
        code = code + tmpCode


        #Calculating TPC
        tmpCode = [0] * 8000
        for j in range(len(sequence) - 3 + 1):
            tmpCode[AADict_AC[sequence[j]] * 400 + AADict_AC[sequence[j+1]]*20 + AADict_AC[sequence[j+2]]] = tmpCode[AADict_AC[sequence[j]] * 400 + AADict_AC[sequence[j+1]]*20 + AADict_AC[sequence[j+2]]] +1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        #appending TPC
        code = code + tmpCode

        #Calculating GAAC
        count = Counter(sequence)
        myDict_gaac = {}
        for key in gaac_groupKey:
            for aa in gaac_group[key]:
                myDict_gaac[key] = myDict_gaac.get(key, 0) + count[aa]

        for key in gaac_groupKey:
            #Appending GAAC
            code.append(myDict_gaac[key]/len(sequence))

        #Calculating Moran
        with open("./data/AAidx2.txt") as f:
            records = f.readlines()[1:]
        myDict = {}
        for i in records:
            array = i.rstrip().split('\t')
            myDict[array[0]] = array[1:]

        AAidx = []
        AAidxName = []
        for i in props:
            if i in myDict:
                AAidx.append(myDict[i])
                AAidxName.append(i)

        AAidx1 = np.array([float(j) for i in AAidx for j in i])
        AAidx = AAidx1.reshape((len(AAidx), 20))

        propMean = np.mean(AAidx,axis=1)
        propStd = np.std(AAidx, axis=1)

        for i in range(len(AAidx)):
            for j in range(len(AAidx[i])):
                AAidx[i][j] = (AAidx[i][j] - propMean[i]) / propStd[i]

        index = {}
        for i in range(len(AA)):
            index[AA[i]] = i

        moran_code = []
        geary_code = []
                        
        N = len(sequence)
        for prop in range(len(props)):
            xmean = sum([AAidx[prop][index[aa]] for aa in sequence]) / N
            for n in range(1, nlag + 1):
                if len(sequence) > nlag:
                    fenzi = sum([(AAidx[prop][index.get(sequence[j], 0)] - xmean) * (AAidx[prop][index.get(sequence[j + n], 0)] - xmean) for j in range(len(sequence) - n)]) / (N - n)
                    fenmu = sum([(AAidx[prop][index.get(sequence[j], 0)] - xmean) ** 2 for j in range(len(sequence))]) / N
                    rn = fenzi / fenmu
                    rn2 = (N-1)/(2*(N-n)) * ((sum([(AAidx[prop][index.get(sequence[j], 0)] - AAidx[prop][index.get(sequence[j + n], 0)])**2 for j in range(len(sequence)-n)])) / (sum([(AAidx[prop][index.get(sequence[j], 0)] - xmean) ** 2 for j in range(len(sequence))])))
                else:
                    rn = 'NA'
                    rn2 = 'NA'
                moran_code.append(rn)
                geary_code.append(rn2)
        #Appending Moran
        code = code + moran_code

        #Calculating CTD
        code_ctdc = []
        code_ctdd = []
        code_ctdt = []

        aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]

        for p in property:
            c1 = CountCTDC(group1[p], sequence) / len(sequence)
            c2 = CountCTDC(group2[p], sequence) / len(sequence)
            c3 = 1 - c1 - c2
            code_ctdc = code_ctdc + [c1, c2, c3]
            code_ctdd = code_ctdd + CountCTDD(group1[p], sequence) + CountCTDD(group2[p], sequence) + CountCTDD(group3[p], sequence)
            c1221, c1331, c2332 = 0, 0, 0
            for pair in aaPair:
                if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
                    c1221 = c1221 + 1
                    continue
                if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
                    c1331 = c1331 + 1
                    continue
                if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
                    c2332 = c2332 + 1
            code_ctdt = code_ctdt + [c1221/len(aaPair), c1331/len(aaPair), c2332/len(aaPair)]

        #appending CTDC
        code = code + code_ctdc + code_ctdd + code_ctdt




        #Calculating and appending CTriad
        code = code + CalculateCTriad(sequence, 0, features, AADict)



        #Calculating PAAC
        theta_PAAC = []
        theta_APAAC = []


        for n in range(1, lambdaValue + 1):
            theta_PAAC.append(sum([Rvalue(sequence[j], sequence[j + n], AADict2, AAProperty1) for j in range(len(sequence) - n)]) / (len(sequence) - n))
            
        myDict = {}
        for aa in AA:
            myDict[aa] = sequence.count(aa)

        #appending PAAC
        code = code + [myDict[aa] / (1 + w * sum(theta_PAAC)) for aa in AA]
        code = code + [(w * j) / (1 + w * sum(theta_PAAC)) for j in theta_PAAC]


        #Calculating other physiochemical features
        code.append(ProteinAnalysis(sequence).molecular_weight())
        code.append(ProteinAnalysis(sequence).isoelectric_point())
        code.append(ProteinAnalysis(sequence).gravy())
        code.append(ProteinAnalysis(sequence).instability_index())

        t = ProteinAnalysis(sequence).secondary_structure_fraction()
        for i in t:
            code.append(i)

        t = ProteinAnalysis(sequence).molar_extinction_coefficient()
        for i in t:
            code.append(i)

        code.append(length)
        
        #Calculating Transmembrane
        if type(trans) == str:
            code.append(1)
        else:
            code.append(0)

        #Calculating Bindings
        if type(dna) == str:
            code.append(1)
        else:
            code.append(0)

        if type(metal) == str:
            code.append(1)
        else:
            code.append(0)

        if type(ncl) == str:
            code.append(1)
        else:
            code.append(0)
            
        for i in range(0,len(subcell_list)):
            if type(sub) == str:
                sub = sub.lower()
                if subcell_list[i] in sub:
                    code.append(1)
                else:
                    code.append(0)
            else:
                code.append(0)
        
        #Calculating Motifs
        for i in range(0,len(motif_list)):
            c = count_motif(sequence,motif_list[i])
            code.append(c)


        #Appending GoData ClassLabels
        for i in range(0,len(GoData_list)):
            if GoData_list[i] in go_id:
                code.append(1)
            else:
                code.append(0)

        encodings.append(code)        
        writer = csv.writer(feature_file)
        writer.writerow(code)
    
    
