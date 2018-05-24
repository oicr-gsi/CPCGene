# -*- coding: utf-8 -*-
"""
Created on Wed May  9 15:32:05 2018

@author: rjovelin
"""

# use this script to check md5 and new headers on reheadered bams
# precondition: samtools need to be loaded
# module load samtools/1.5

import os
import sys
import subprocess
import argparse
import yaml


# use this function to match sample and bam name to full path
def MatchBamPath(folder):
    '''
    (str) -> tuple
    Look for the bam file in folder and return a tuple with the sample name,
    the  bam file name and the full path to the bam
    '''
    
    # make a list of files in folder
    L = [i for i in os.listdir(folder) if os.path.isfile(os.path.join(folder, i))]
    # loop over files
    sample, bam, path = '', '', ''
    if len(L) != 0:
        for filename in L:
            sample, bam, path = '', '', ''
            # check if file is bam
            if 'reheadered.bam' in filename:
                if filename[filename.index('reheadered.bam'):] == 'reheadered.bam':
                    sample = filename[:filename.index('_realigned_')]
                    bam = filename
                    path = os.path.join(folder, filename)
                    assert os.path.isfile(path)
                    break
    return (sample, bam, path)

# use this functions to match samples of bams with the bam full path
def GetBamPath(BamDir, Key):
    '''
    (str) -> dict
    Take the path to directories with bams for all samples and return a dictionary
    of sample: path to bam pairs for sample in Dataset if Key = sample, and a dictionary
    of reheadered bam names: path to reheadered bam for sample in Dataset if Key = file_final or file_temp
    '''
    # create a dict {sample: bam full path} pairs or {bam name: bam full path} pairs 
    BamFiles = {}
    assert os.path.isdir(BamDir)
    # make a list of CPCG sample dir
    L = [i for i in os.listdir(BamDir) if i.startswith('CPCG') and os.path.isdir(os.path.join(BamDir, i))]
    # loop over sample dirs    
    for i in L:
        # get the sample dir
        samplefolder = os.path.join(BamDir, i)
        assert os.path.isdir(samplefolder)
        # need to dig one more level to get the normal and tumor bams
        # make a list of subdirectories
        K = [item for item in os.listdir(samplefolder) if os.path.isdir(os.path.join(samplefolder, item))]
        # loop over subfolders if they exist
        if len(K) != 0:
            for m in range(len(K)):
                subfolder = os.path.join(samplefolder, K[m])
                assert os.path.isdir(subfolder)
                # check that subfolder is not empty
                if len(os.listdir(subfolder)) != 0:
                    if Key == 'sample':
                        sample, bam, path = MatchBamPath(subfolder)
                        assert sample != '' and bam != '' and path != ''
                        assert sample not in BamFiles
                        BamFiles[sample] = path
                    elif Key == 'file_final':
                        sample, bam, path = MatchBamPath(subfolder)
                        assert bam not in BamFiles
                        assert sample != '' and bam != '' and path != ''
                        BamFiles[bam] = path
                    elif Key == 'file_temp':
                        # get the bam in the reheader folder
                        # check that reheader folder is present
                        if 'reheader' in os.listdir(subfolder):
                            # get reheader folder name
                            reheaderfolder = os.path.join(subfolder, 'reheader')
                            assert os.path.isdir(reheaderfolder)
                            sample, bam, path = MatchBamPath(reheaderfolder)
                            assert bam not in BamFiles
                            # if bam has been moved already, directory is empty
                            # populate dict only if bam exists
                            if sample != '' and bam != '' and path != '':
                                BamFiles[bam] = path
                            else:
                                print('{0} is empty'.format(reheaderfolder))
    if Key == 'sample':
        print('matched samples to full path of original bams')    
    else:
        print('matched file names to full path of reheadered bams')
    return BamFiles

            
# use this function to extract the md5
def Getmd5(checksum):
    '''
    (file) -> str
    Take a md5 file and return the extracted md5 string
    '''
    infile = open(checksum)
    content = infile.read().rstrip().split()[0].strip()
    return content

# use this function to collect fields for bam header
def GrabRGBamHeader(filename):
    '''
    (file) -> dict 
    Take a bam file and return a dictionary with {pu: field list} pairs for given bam
    '''
    RG = {}
    header = subprocess.check_output('samtools view -H ' + filename, shell = True)
    header = header.decode('ascii')
    L = header.split('\n')    
    K = [i for i in L if i.startswith('@RG')]
    for s in K:
        name, cn, pu, ID, lb = '' ,'', '', '', '' 
        s = s.split()
        s = [i for i in s if ':' in i]
        for i in s:
            if i.startswith('SM') :
                i = i.split(':')
                name = i[1]
            elif i.startswith('CN'):
                i = i.split(':')
                cn = i[1]
            elif i.startswith('PU'):
                i = i.split(':')
                pu = i[1]
            elif i.startswith('ID'):
                i = i.split(':')
                ID = i[1]
            elif i.startswith('LB'):
                i = i.split(':')
                lb = i[1]
        assert name != '' and cn != '' and pu != '' and ID != '' and lb != ''
        # initialize inner dict
        assert pu not in RG
        RG[pu] = [pu, cn, ID, lb]
    return RG


# use this function to list of the 500PG files
def Grab500PGFiles(BamYaml):
    '''
    (file) -> list
    Take the yaml files with 500PG samples and return a list of files (full path)
    that are part of 500PG
    '''
    # make a list of bams from 500PG
    infile = open(BamYaml)
    Samples500PG = yaml.load(infile)
    infile.close()
    Bams500PG = []
    for sample in Samples500PG:
        Bams500PG.append(Samples500PG[sample]['normal']['path'])
        Bams500PG.append(Samples500PG[sample]['tumour']['path'])
    return Bams500PG


# use this function to extract the readg group info from the readgroup file
def ParseReadGroupFile(ReadGroupFile, FileType):
    '''
    (file) -> dict
    Take the file with read group information and return a dictionary with readgroup
    fields for each PU for each sample or file name depending on the Key value 
    '''
    
    # make a list of fields to extract from file
    Fields = ['PU', 'CN', 'LB', 'ID']
    if FileType == 'summary':
        for i in range(len(Fields)):
            Fields[i] = Fields[i] + '_rg'
        
    # create a dict with {sample: {PU:[CN, ID, LB]}} pairs from the readgroup file
    Rg = {}
    infile = open(ReadGroupFile)
    header = infile.readline().rstrip().split('\t')
    for line in infile:
        MasterKey = ''
        line = line.rstrip()
        if line.rstrip() != '':
            line = line.split()
            if FileType == 'readgroup':
                sample = line[header.index('SM')]
            elif FileType == 'summary':
                filename = line[1]
                filename = filename[filename.rfind('/')+1:]
            pu = line[header.index(Fields[0])]
            cn = line[header.index(Fields[1])]
            lb = line[header.index(Fields[2])]
            ID = line[header.index(Fields[3])]
            
            # initialize inner dict
            if FileType == 'readgroup':
                MasterKey = sample
            elif FileType == 'summary':
                MasterKey = filename
            assert MasterKey != ''
            if MasterKey not in Rg:
                Rg[MasterKey] = {}
            if pu in Rg[MasterKey]:
                if [pu, cn, ID, lb] not in Rg[MasterKey][pu]:
                    Rg[MasterKey][pu].append([pu, cn, ID, lb])
            else:
                Rg[MasterKey][pu] = [[pu, cn, ID, lb]]
                    
    infile.close()
    # check that 1 CN is recorded per PU
    for i in Rg:
        for pu in Rg[i]:
            assert len(Rg[i][pu]) ==1
            # reassign list to value
            Rg[i][pu] = Rg[i][pu][0]
    if FileType == 'readgroup':
        print('collected CN from readgroup file')
    elif FileType == 'summary':
        print('collected fields from summary file')
    return Rg

    
# use this function to review the readgroups in bams and file and write summary to files
def ReviewReadGroup(args):
    '''
    (list) -> files
    Take the list of argument from the command line, review readgroup info in bam
    headers and readgroup file and write a summary file
    '''
           
    # 1) make a list of bams from 500PG
    Bams500PG = Grab500PGFiles(args.yaml)
    
    # 2) collect fields for each samples from the readgroup file
    RgInfo = ParseReadGroupFile(args.readgroup, 'readgroup')
   
    # 3) Match bams to samples for each dataset
    path = '/.mounts/labs/cpcgene/private/samples/Analysis/wgs/bwa-mem/{0}/gatk/{1}/bam/realign-recal/'
    # make a list of bam directories
    BamDirs = [path.format('0.7.15', '3.7.0'), path.format('0.7.15', '3.6.0'), path.format('0.7.15', '3.5.0'),
               path.format('0.7.12', '3.5.0'), path.format('0.7.12', '3.4.0')]
    # make a list of dicts with bam: paths
    BamsDataset = [GetBamPath(BamDirs[i], 'sample') for i in range(len(BamDirs))]

    # 4) Collect fields for each sample from the bam header
    BamRGInfo = []
    for i in range(len(BamsDataset)):
        # create a dict {sample: {PU:[PU, CN, ID, LB]}}
        datasetRG = {}
        for sample in BamsDataset[i]:
            assert sample not in datasetRG
            datasetRG[sample] = GrabRGBamHeader(BamsDataset[i][sample])
        BamRGInfo.append(datasetRG)
    
    # 5) Compare cn, id, lb from readgtoup and bam headers
    # create lists of dicts for samples with matching pus and non-matching pus    
    MatchedHeaders, Absent = [], []
    # check bam headers and read group header for each dataset
    for i in range(len(BamRGInfo)):
        # create a dict to store bam and rg field values
        Matched, NotMatched = {}, {}
        # loop over samples for that dataset
        for sample in BamRGInfo[i]:
            # check that sample is in readgroup
            assert sample in RgInfo
            # initialize inner dict
            Matched[sample] = {}
            # compare the values of each PU for each sample
            for PU in BamRGInfo[i][sample]:
                # some sample have non-matching PUs between bams readgroup
                # create a key based on PU to match bamcenters[i][sample] and rgcenters[sample]
                newPu = ''
                if PU not in RgInfo[sample]:
                    #print(sample, PU, list(RgCenters[sample].keys()))
                    # loop over PU, check if PU from bam is contained in PU from rg
                    for j in list(RgInfo[sample].keys()):
                        if PU in j or j in PU:
                            newPu = j
                            break
                # check if matching PU is found
                else:
                    newPu = PU
                # check if a matching PU have been found
                if newPu == '':
                    # PUs do not match between bam and rg
                    NotMatched[sample] = BamsDataset[i][sample]
                else:
                    assert newPu in RgInfo[sample]
                    # record bam full path, bam fields and rg fields
                    Matched[sample][PU] = [BamsDataset[i][sample]]
                    Matched[sample][PU].extend(BamRGInfo[i][sample][PU])
                    Matched[sample][PU].extend(RgInfo[sample][newPu])
                    # compare the values of the same PU between bam header and readgroup
                    if BamRGInfo[i][sample][PU] != RgInfo[sample][newPu]:
                        # tag with no match
                        Matched[sample][PU].append('non-matching')
                    else:
                        Matched[sample][PU].append('matched')
    
        # remove sample from matched if at least 1 PU was never found
        to_remove = list(NotMatched.keys())
        if len(to_remove) != 0:
            for sample in to_remove:
                if sample in Matched:
                    del Matched[sample]
        # add dicts to respective lists
        MatchedHeaders.append(Matched)
        Absent.append(NotMatched)                
            
    # 6) check that all sample and PUs are recorded
    for i in range(len(MatchedHeaders)):
        assert len(MatchedHeaders[i]) + len(Absent[i]) == len(BamRGInfo[i]) == len(BamsDataset[i])    
        for sample in MatchedHeaders[i]:
            assert len(MatchedHeaders[i][sample]) == len(BamRGInfo[i][sample])
    
    # 7) Write results to files
    newfile = open(args.output, 'w')
    OutputHeader = ['Sample', 'Full_path', 'PU_Bam', 'CN_Bam', 'ID_Bam', 'LB_Bam', 'PU_rg', 'CN_rg', 'ID_rg', 'LB_rg', 'Match', '500PG'] 
    newfile.write('\t'.join(OutputHeader) + '\n')
    for i in range(len(MatchedHeaders)):
        for sample in MatchedHeaders[i]:
            for PU in MatchedHeaders[i][sample]:
                line = [sample]
                line.extend(MatchedHeaders[i][sample][PU])
                # check if bam is part of 500PG
                if MatchedHeaders[i][sample][PU][0] in Bams500PG:
                    line.append('Y')
                else:
                    line.append('N')
                newfile.write('\t'.join(line) + '\n')
    newfile.close()

    newfile = open(args.error, 'w')
    MissingHeader = ['Sample', 'Full_path', '500PG']
    newfile.write('\t'.join(MissingHeader) + '\n')
    for i in range(len(Absent)):
        if len(Absent[i]) != 0:
            for sample in Absent[i]:
                line = [sample, Absent[i][sample]]
                if Absent[i][sample] in Bams500PG:
                    line.append('Y')
                else:
                    line.append('N')
                newfile.write('\t'.join(line) + '\n')
    newfile.close()



# use this function to compare bams pre and post-reheadering
# compare md5sums and compare bam header and readgroup file
def ComparePrePost(args):
    '''    
    (list) -> file 
    Take the list of command line arguments and write a file with readgroup info
    from the readgroup file and post-reheadering bam and md5sums from the original
    and post reheadering bams
    '''
   
    # get the path the bams
    BamDir = '/.mounts/labs/cpcgene/private/samples/Analysis/wgs/bwa-mem/{0}/gatk/{1}/bam/realign-recal/'.format(args.bwa, args.gatk)
    
    # 1) make a dict new reheadered bams  {bam: full path} 
    if args.summary == 'temp':
        # write temporary summary to check if reheadered bams can replace original bams
        # make a dict of bams: full path with reheadered bams in the reheader folder
        NewBams = GetBamPath(BamDir, 'file_temp')
    elif args.summary == 'final':
        # write final summary once all bams have been replaced by reheadered bams
        # make a dict of bams: full path with reheadered bams 
        NewBams = GetBamPath(BamDir, 'file_final')
    
    # 2) match the new bams with their md5, check that md5 are the same
    Md5 = {}
    # get the checksum directory
    ChecksumDir = os.path.join('/.mounts/labs/gsiprojects/boutroslab/CPCGene/PIPELINES/REVIEW_BAM_HEADERS/CORRECT_HEADERS/bwa-mem_{0}_gatk_{1}/checksums'.format(args.bwa, args.gatk))
    assert os.path.isdir(ChecksumDir)
    L = os.listdir(ChecksumDir)
    if len(L) != 0:
        for bam in NewBams:
            Md5[bam] = ['', '']
            for filename in L:
                if bam in filename:
                    assert 'md5' in filename
                    if 'post_reheadering' in filename:
                        Md5[bam][1] = Getmd5(os.path.join(ChecksumDir, filename))
                    else:
                        Md5[bam][0] = Getmd5(os.path.join(ChecksumDir, filename))
                    
    # check that md5 have been collected for all bams
    assert set(Md5.keys()) == set(NewBams.keys())
    for bam in Md5:
        assert Md5[bam][0] != '' and Md5[bam][1] != ''
    print('extracted md5')

    # 3) collect fields for each file name from the readgroup file
    # create a dict with {bam: {PU:[CN, ID, LB]}} pairs from the summary file
    RGInfo = ParseReadGroupFile(args.review, 'summary')
        
    # 4) write summary file
    if args.summary == 'temp':
        Outputfile = 'TempStatus_bwa_{0}_gatk{1}.txt'.format(args.bwa, args.gatk)
    elif args.summary == 'final':
        Outputfile = 'StatusReheader_bwa_{0}_gatk{1}.txt'.format(args.bwa, args.gatk)
        
    newfile = open(Outputfile, 'w')
    Header = ['sample', 'filename', 'fullpath', 'md5_original', 'md5_post', 'md5_match', 'PU_rg', 'CN_rg', 'ID_rg', 'LB_rg', 'PU_bam', 'CN_bam', 'ID_bam', 'LB_bam', 'RG_match'] 
    newfile.write('\t'.join(Header) + '\n')

    for bam in NewBams:
        # collect readgroup info, compare to readgroup file
        RG = GrabRGBamHeader(NewBams[bam])
        # check that bam has info from RG file
        assert bam in RGInfo
        # compare all PU
        assert set(RGInfo[bam].keys()) == set(RG.keys())
        for pu in RG:
            line = []
            assert RG[pu] == RGInfo[bam][pu]
            # get sample name
            sample = bam[:bam.index('_')]
            # get full path
            path = NewBams[bam]
            # get md5 original, post-reheadering and md5 match
            md5_original, md5_post, md5_match = Md5[bam][0], Md5[bam][1], Md5[bam][0] == Md5[bam][1]
            # get RG file fields
            PU_rg, CN_rg, ID_rg, LB_rg = RGInfo[bam][pu]
            # get bam fields
            PU_bam, CN_bam, ID_bam, LB_bam = RG[pu]
            # check if rg fields match
            fields_match = RG[pu] == RGInfo[bam][pu]
            assert md5_match == True and fields_match == True
            line = [sample, bam, path, md5_original, md5_post, str(md5_match),
                    PU_rg, CN_rg, ID_rg, LB_rg, PU_bam, CN_bam, ID_bam, LB_bam, str(fields_match)]
            newfile.write('\t'.join(line) + '\n')
    newfile.close()


# use this fucntion to edit the bam header
def EditBamHeader(bam, BamHeader, BamRGHeaderFile):
    '''
    (str, str, file) -> str
    Take the full path string of a bam, the bam header string and the file
    with read group information
    '''
    
    # make a list strings from bamheader
    NewHeader = BamHeader.split('\n')
    
    infile = open(BamRGHeaderFile)
    header = infile.readline().rstrip().split('\t')
    for line in infile:
        if bam in line:
            line = line.rstrip().split('\t')
            # extract fields from bam and from RG
            Full_path = line[header.index('Full_path')]
            PU_Bam, CN_Bam = line[header.index('PU_Bam')], line[header.index('CN_Bam')]
            ID_Bam, LB_Bam = line[header.index('ID_Bam')], line[header.index('LB_Bam')]
            PU_rg, CN_rg = line[header.index('PU_rg')], line[header.index('CN_rg')]
            ID_rg, LB_rg = line[header.index('ID_rg')], line[header.index('LB_rg')]
            # create field with fields in bam and readgroup file
            BamField = [PU_Bam, CN_Bam, ID_Bam, LB_Bam]
            RGField = [PU_rg, CN_rg, ID_rg, LB_rg]
            # check that fields are present in bamheader
            for i in BamField:
                assert i in BamHeader
            # loop over strings in list, check if string contain @RG, 
            for i in range(len(NewHeader)):
                if '@RG' in NewHeader[i]:
                    if PU_Bam in NewHeader[i] and CN_Bam in NewHeader[i] and ID_Bam in NewHeader[i] and LB_Bam in NewHeader[i]:
                        for j in range(len(BamField)):
                            NewHeader[i] = NewHeader[i].replace(BamField[j], RGField[j])
    infile.close()
    NewHeader = '\n'.join(NewHeader)
    return NewHeader
  


# use this function to set scripts to reheader the bams
def ReheaderBams(args):
    '''
    (list) -> files
    Take the list of command line arguments and write scripts to reheader bams
    '''

    # make a list of bams to reheader
    AllBams = []        
    infile = open(args.review)
    infile.readline()
    for line in infile:
        if 'CPCG' in line:
            line = line.rstrip().split('\t')
            # check if fields are matching
            if line[-2] == 'non-matching':
                # collect full path of the bam
                AllBams.append(line[1])
    infile.close()

    # make a list of bams for versions of bwa and gatk
    MyBams = [bam for bam in AllBams if args.bwa in bam and args.gatk in bam]

    # get the working directory and output directory
    WorkDir = '/.mounts/labs/gsiprojects/boutroslab/CPCGene/PIPELINES/REVIEW_BAM_HEADERS/CORRECT_HEADERS/'
    DirOut = os.path.join(WorkDir, 'bwa-mem_{0}_gatk_{1}'.format(args.bwa, args.gatk))
    # check if DirOut exists
    if os.path.isdir(DirOut) == False:
        print('please verify that {0} exists or provide valid bwa and gatk versions'.format(DirOut))
        sys.exit(2)
    # get the log, edited bam headers and checksums directories
    ChecksumDir = os.path.join(DirOut, 'checksums')
    LogDir = os.path.join(DirOut, 'log')
    EditedHeaderDir = os.path.join(DirOut, 'edited_headers')
    # create dirs where logs and checksums are saved if they don't already exists
    if os.path.isdir(ChecksumDir) == False:
        os.mkdir(ChecksumDir)
    if os.path.isdir(LogDir) == False:
        os.mkdir(LogDir)
    if os.path.isdir(EditedHeaderDir) == False:
        os.mkdir(EditedHeaderDir)

    # loop over bams
    for bam in MyBams:
        BamHeader = ''
        # extract the file name and directory from the full path
        filename = bam[bam.rfind('/')+1:]
        bamDir = bam[:bam.rfind('/')]
        # create a qsub script
        QsubScript = os.path.join(DirOut, filename + '.reheader.sh')
        qsubFile = open(QsubScript, 'w')
        # check if reheader directory exists
        reheaderDir = os.path.join(bamDir, 'reheader')
        if os.path.isdir(reheaderDir) == False:
            qsubFile.write('mkdir {0}'.format(reheaderDir) + '\n')

        # get the header
        BamHeader = subprocess.check_output('samtools view -H {0}'.format(bam), shell=True).decode('utf-8')
        assert BamHeader != ''    
        # edit header, replace fields in header with fields from readgroup file
        NewBamHeader = EditBamHeader(bam, BamHeader, args.review)    
   
        # write edited bam header to file
        outheader = os.path.join(EditedHeaderDir, filename + '.edited_header.sam')
        newfile = open(outheader, 'w')
        newfile.write(NewBamHeader + '\n')
        newfile.close()
    
        # do a checksum on the body of the original bam
        original_md5 = os.path.join(ChecksumDir, filename + '.md5')
        line = 'qsub -cwd -b y -N md5.bwa-mem_{0}_gatk_{1}.{2} -e log -o log \"module load samtools;samtools view {3} | md5sum > {4}\"'.format(args.bwa, args.gatk, filename, bam, original_md5)
        qsubFile.write(line + '\n')
    
        # reheader the bam
        NewBamFile = os.path.join(reheaderDir, filename)
        line = 'qsub -cwd -b y -N reheader.bwa-mem_{0}_gatk_{1}.{2} -e log -o log \"module load samtools;samtools reheader -i {3} {4} > {5}\"'.format(args.bwa, args.gatk, filename, outheader, bam, NewBamFile)  
        qsubFile.write(line + '\n')
    
        # do a checksum on the body of the new bam
        reheadered_md5 = os.path.join(ChecksumDir, filename + '.post_reheadering.md5')
        line = 'qsub -cwd -b y -N md5.bwa-mem_{0}_gatk_{1}.{2} -hold_jid reheader.bwa-mem_{3}_gatk_{4}.{5} -e log -o log \"module load samtools;samtools view {6} | md5sum > {7}\"'.format(args.bwa, args.gatk, filename + '.post_reheadering', args.bwa, args.gatk, filename, NewBamFile, reheadered_md5)
        qsubFile.write(line + '\n')
        qsubFile.close()


# use this function to move reheadered bams to their 
def MoveBamsFinal(args):
    '''
    (list) -> None
    Take the list of command line arguments and move the reheadred bams to their
    final destination, replacing the original bams if the md5sums between original
    and reheadred bams match and if readgroup info match between readgroup file and reheadered bams
    '''
    
    # get the directory containing the original bams
    BamDir = '/.mounts/labs/cpcgene/private/samples/Analysis/wgs/bwa-mem/{0}/gatk/{1}/bam/realign-recal/'.format(args.bwa, args.gatk)
    
    # check that a temporary file with bam status exist
    StatusFile = 'TempStatus_bwa_{0}_gatk{1}.txt'.format(args.bwa, args.gatk)
    infile = open(StatusFile)
    header = infile.readline().rstrip().split('\t')
    # create a dict with bam : path to reheadered bam if md5sums and readgroup match
    Reheadered = {}
    NoMatch = []
    for line in infile:
        if 'CPCG' in line:
            line = line.rstrip().split('\t')
            # check that md5sums match
            md5_original = line[header.index('md5_original')]
            md5_post = line[header.index('md5_post')]
            md5_match = bool(line[header.index('md5_match')])
            RG = line[header.index('PU_rg'): header.index('PU_bam')]
            rg = line[header.index('PU_bam'):header.index('RG_match')]
            RG_match = bool(line[header.index('RG_match')])
            if md5_match == True and RG_match == True:
                assert md5_original == md5_post and RG == rg
                # get the bam name and its path
                bam, path = line[header.index('filename')], line[header.index('fullpath')]
                assert os.path.isfile(path)
                # bams may be recorded on multiple lines because PUs are the unique fields
                if bam in Reheadered:
                    assert Reheadered[bam] == path
                else:
                    Reheadered[bam] = path
            else:
                NoMatch.append(bam)
    infile.close()            
                
    # remove bams if md5 of RG did not match
    NoMatch = list(set(NoMatch))
    if len(NoMatch) != 0:
        for bam in NoMatch:
            print('no match for bam {0}'.format(bam))
            del Reheadered[bam]
    
    # make a dictionary of original bams: path pairs
    OriginalBams = GetBamPath(BamDir, 'file_final')
    
    # check that all reheadered bams can be match to their final destination
    for bam in Reheadered:
        print('moving bam {0}'.format(bam)) 
        OldFile = OriginalBams[bam]
        NewFile = Reheadered[bam]
        assert os.path.isfile(OldFile)
        assert os.path.isfile(NewFile)
        # replace original bam with reheadered bam
        os.system('mv {0} {1}'.format(NewFile, OldFile))
        # change group --> cpcgene
        os.system('chgrp cpcgene {0}'.format(OldFile))

if __name__ == '__main__':
    
    # create top-level parser
    parser = argparse.ArgumentParser(prog = 'ReviewBamHeaders.py', description='review read group fields in bams and readgroup file')
    subparsers = parser.add_subparsers(title='sub-commands', description='valid sub-commands', help = 'sub-commands help')
    
    # review sub-commands
    Review_parser = subparsers.add_parser('Review', help ='Review read group fields in bams and readgroup file')
    Review_parser.add_argument('readgroup', help='read group file')
    Review_parser.add_argument('-o', '--Out', dest='output', default='BamRGHeaders.txt', help='summary file with readgroup info. default is ./BamRGHeaders.txt')
    Review_parser.add_argument('-e', '--Err', dest='error', default='ProblematicBams.txt', help='bams that cannot be compared between headers and readgroup file. default is ./ProblematicBams.txt')
    Review_parser.add_argument('-y', '--Yaml', dest='yaml', default='/.mounts/labs/cpcgene/private/projects/500pg/data/samples/500pg_bam_all.yaml', help='yaml files with 500PG samples. default is /.mounts/labs/cpcgene/private/projects/500pg/data/samples/500pg_bam_all.yaml')
    Review_parser.set_defaults(func=ReviewReadGroup)

    # Compare bams post and pre-reheadering sub-commands
    Compare_parser = subparsers.add_parser('Compare', help='Compare bams post and pre-reheadering')               
    Compare_parser.add_argument('-r', '--Review', dest='review', default='/.mounts/labs/gsiprojects/boutroslab/CPCGene/PIPELINES/REVIEW_BAM_HEADERS/BamRGHeaders.txt', help='review file with bam and readgroup fields. default is /.mounts/labs/gsiprojects/boutroslab/CPCGene/PIPELINES/REVIEW_BAM_HEADERS/BamRGHeaders.txt')
    Compare_parser.add_argument('-b', '--Bwa', dest='bwa', default='0.7.15', help='Bwa version. default 0.7.15', choices=['0.7.12', '0.7.15'])
    Compare_parser.add_argument('-g', '--Gatk', dest='gatk', default='3.7.0', help='Gatk version. default 3.7.0', choices=['3.4.0', '3.5.0', '3.6.0', '3.7.0'])
    Compare_parser.add_argument('-s', '--Summary', dest='summary', help='Write temporary or final summary with md5sums and readgroup info', choices=['final', 'temp'])
    Compare_parser.set_defaults(func=ComparePrePost)
    
    # reheader bams sub-commands
    Reheader_parser = subparsers.add_parser('Reheader', help='Reheader bams with info from readgroup file')
    Reheader_parser.add_argument('-b', '--Bwa', dest='bwa', default='0.7.15', help='Bwa version. default 0.7.15', choices=['0.7.12', '0.7.15'])
    Reheader_parser.add_argument('-g', '--Gatk', dest='gatk', default='3.7.0', help='Gatk version. default 3.7.0', choices=['3.4.0', '3.5.0', '3.6.0', '3.7.0'])
    Reheader_parser.add_argument('-r', '--Review', dest='review', default='/.mounts/labs/gsiprojects/boutroslab/CPCGene/PIPELINES/REVIEW_BAM_HEADERS/BamRGHeaders.txt', help='review file with bam and readgroup fields. default is /.mounts/labs/gsiprojects/boutroslab/CPCGene/PIPELINES/REVIEW_BAM_HEADERS/BamRGHeaders.txt')
    Reheader_parser.set_defaults(func=ReheaderBams)
    
    # move bams sub-commands
    MoveBam_parser = subparsers.add_parser('MoveBam', help='Move reheadered bams to final directory')
    MoveBam_parser.add_argument('-b', '--Bwa', dest='bwa', default='0.7.15', help='Bwa version. default 0.7.15', choices=['0.7.12', '0.7.15'])
    MoveBam_parser.add_argument('-g', '--Gatk', dest='gatk', default='3.7.0', help='Gatk version. default 3.7.0', choices=['3.4.0', '3.5.0', '3.6.0', '3.7.0'])
    MoveBam_parser.set_defaults(func=MoveBamsFinal)
        
    # get arguments from the command line
    args = parser.parse_args()
    # pass the args to the default function
    args.func(args)

