import os
import sys
import time
import math
import globalConstant
import numpy as np

def  make_Tri_BiPSDP_vector(seqline,posi_tri_fre,nega_tri_fre,posi_six_fre,nega_six_fre,lenSeq,methtype):
    vector_Tri_BIPSDP = []
    flag = True
    for i in range(lenSeq-5):
        for j in range(i+3,lenSeq-2):
            first_Tri=seqline[i:i+3]
            second_Tri=seqline[j:j+3]
            six_united=first_Tri+second_Tri
            if methtype == 'RNA':
                first_Tri_idx=globalConstant.tri_RNA.index(first_Tri)
                second_Tri_idx=globalConstant.tri_RNA.index(second_Tri)
                six_idx=globalConstant.six_RNA.index(six_united)
            else:
                first_Tri_idx = globalConstant.tri_DNA.index(first_Tri)
                second_Tri_idx = globalConstant.tri_DNA.index(second_Tri)
                six_idx=globalConstant.six_DNA.index(six_united)
            first_Tri_fr_posi=posi_tri_fre[first_Tri_idx][i]
            second_Tri_fr_posi=posi_tri_fre[second_Tri_idx][j]
            six_fr_posi=posi_six_fre[six_idx][int((2*lenSeq-9-i)*i/2+j-i-3)]
            first_Tri_fr_nega = nega_tri_fre[first_Tri_idx][i]
            second_Tri_fr_nega = nega_tri_fre[second_Tri_idx][j]
            six_fr_nega = nega_six_fre[six_idx][int((2*lenSeq-9-i)*i/2+j-i-3)]
            if first_Tri_fr_posi==0 or second_Tri_fr_posi==0 or six_fr_posi==0:
                posi_value=0
            else:
                posi_value=math.log(first_Tri_fr_posi)*math.log(second_Tri_fr_posi)/math.log(six_fr_posi)
            if first_Tri_fr_nega==0 or second_Tri_fr_nega==0 or six_fr_nega==0:
                nega_value=0
            else:
                nega_value = math.log(first_Tri_fr_nega) * math.log(second_Tri_fr_nega) / math.log(six_fr_nega)
            vector_Tri_BIPSDP=np.hstack((vector_Tri_BIPSDP, posi_value-nega_value))
    return vector_Tri_BIPSDP


def calculate_frequency_single(first, nucleStr, data):
    '''
    :param first: the first nucleotide
    :param nucleStr: the string of nucleotide
    :param data: positive or negative data
    :return: The frequency of nucleotides at the current position in positive or negative dataset
    '''
    singleNucle = []
    for line in data:
        singleNucle.append(line[first])
    fre = float(singleNucle.count(nucleStr)) / data.size
    return fre

def calculate_frequency_di(first, second, nucleStr, data):
    '''
    :param first: the first nucleotide
    :param second:  the second nucleotide
    :param nucleStr: the string of dinucleotide
    :param data: positive or negative data
    :return: The frequency of dinucleotides at the current position in positive or negative dataset
    '''
    diNucle = []
    for line in data:
        diNucleStr = line[first] + line[second]
        diNucle.append(diNucleStr)
    fre = float(diNucle.count(nucleStr)) / data.size
    return fre

def make_tri_fre(DataSet,seqlen,methtype):
    vector_fre=np.zeros((64,seqlen-2))
    for line in DataSet:
        for idx in range(seqlen-2):
            str = line[idx:idx+3]
            if methtype == "DNA":
                str_idx=globalConstant.tri_DNA.index(str)
            elif methtype ==  'RNA':
                str_idx=globalConstant.tri_RNA.index(str)
            vector_fre[str_idx][idx]+=1
    vector_fre = vector_fre/DataSet.size
    return vector_fre

def make_six_fre(DataSet, seqlen, methtype):
    vector_six=np.zeros((4096,int((seqlen-5)*(seqlen-4)/2)))
    id=0
    for line in DataSet:
        for idx in range(seqlen-5):
            for jdx in range(idx+3,seqlen-2):
                str = line[idx:idx+3]+line[jdx:jdx+3]
                if methtype == "DNA":
                    str_idx=globalConstant.six_DNA.index(str)
                elif methtype ==  'RNA':
                    str_idx=globalConstant.six_RNA.index(str)
                id+=1
                print(id)
                vector_six[str_idx][int((2*seqlen-9-idx)*idx/2+jdx-idx-3)]+=1
    vector_fre = vector_six/DataSet.size
    return vector_fre

def read_data(filePath, nucleType):
    '''
    :param filePath: the path of data
    :param nucleType: DNA or RNA
    :return: data with list type and label
    '''
    seqList = []
    label = []
    if nucleType == 'RNA':
        startTuple = tuple(globalConstant.single_RNA)
    elif nucleType == 'DNA':
        startTuple = tuple(globalConstant.single_DNA)
    with open(filePath) as files:
        for line in files:
            if line.startswith(startTuple):
                line = line.rstrip('\n')
                lenLine = len(line)
                seqList.append(line)
            else:
                seqName = line
                PosiOrNega = seqName[1]
                if PosiOrNega == 'P' or PosiOrNega == '+':
                    label.append(1)
                else:
                    label.append(2)
    return seqList, label, lenLine

def split_data_posi_nega(data, label):
    '''
    :param data: data with positive and negative
    :param label: positive and negative
    :return: positive and negative data
    '''
    posi_index = [i for i, x in enumerate(label) if x == 1]
    nega_index = [i for i, x in enumerate(label) if x == 2]
    posi_data = data[posi_index]
    nega_data = data[nega_index]
    return posi_data, nega_data

def save_result(filePath, methyType,resultpath,filename):
    '''
    :param filePath: the path of data
    :param alpha: the distance of between first nucleotide and second nucleotide
    :param beta: the distance of between second nucleotide and third nucleotide
    :return: None
    '''
    startTime = time.time()
    splitPath = filePath.split('/')
    nucleType = methyType.split('_')[0]
    filename_fre = filename.split('fasta')[0]
    seqList, label, lenSeq = read_data(filePath, nucleType)
    npSeqList = np.array(seqList)
    posiData, negaData = split_data_posi_nega(npSeqList, label)
    print(len(posiData))
    dataVector_Tri_BiPSDP = []
    seqID = 0
    posi_tri_fre=make_tri_fre(posiData,lenSeq,methyType)
    nega_tri_fre=make_tri_fre(negaData,lenSeq,methyType)
    posi_six_fre = make_six_fre(posiData, lenSeq, methyType)
    nega_six_fre = make_six_fre(negaData, lenSeq, methyType)
    print(1)
    for seqLine in seqList:
        print(seqID)
        seqID = seqID + 1
        lineVector_Tri_BiPSDP = make_Tri_BiPSDP_vector(seqLine,posi_tri_fre,nega_tri_fre,posi_six_fre,nega_six_fre,lenSeq,methyType)
        dataVector_Tri_BiPSDP.append(lineVector_Tri_BiPSDP)
    print('Done. Used time: %.2fs' % (time.time() - startTime))
    resultpaths = os.path.join(resultpath, 'data_'+filename_fre+'_Tri_BiPSDP.csv')
    print(resultpaths)
    np.savetxt(resultpaths, dataVector_Tri_BiPSDP, delimiter=" ")
    arrayLabel = np.array(label)
    resultpath_label = os.path.join(resultpath,  "label_"+filename_fre+ ".csv")
    np.savetxt(resultpath_label, arrayLabel, fmt='%d', delimiter=" ")

if __name__ == '__main__':
    '''
    methyType = 'DNA Sequence'
    if methyType == 'RNA_m6A':
        #fPath = './data/RNA_m6A/non-single-base'  # 非单分辨率
        fPath = './data/RNA_m6A/single-base'  # 单分辨率
        for i, j, k in os.walk(fPath):
            for fName in k:
                print(fName)
                file_path = fPath + '/' + fName
                save_result(file_path, 'RNA', './result', fName)
    elif methyType == 'DNA_4mC':
        fPath = './data/DNA_4mC'
        for i, j, k in os.walk(fPath):
            for fName in k:
                print(fName)
                file_path = fPath + '/' + fName
                save_result(file_path, 'DNA', './result', fName)
    elif methyType == 'DNA_6mA':
        fPath = './data/DNA_6mA'
        for i, j, k in os.walk(fPath):
            for fName in k:
                print(fName)
                file_path = fPath + '/' + fName
                save_result(file_path, 'DNA', './result', fName)
    elif methyType == 'DNA Sequence':
        fPath = './data/DNA Sequence'
        for i, j, k in os.walk(fPath):
            for fName in k:
                print(fName)
                file_path = fPath + '/' + fName
                save_result(file_path, 'DNA', './result', fName)
    #单个数据
    '''
    '''
    file_path = './data/6mA_data_MGF6mARice/Lv.fasta'
    save_result(file_path, 'DNA','./result' , 'Lv')
    '''
    data_names = ['Lv.fasta']
    for data_name in data_names:
        file_path = './data/6mA_data_MGF6mARice/'+data_name
        print(file_path)
        save_result(file_path, 'DNA', './result', data_name)
    print('Done all.')
