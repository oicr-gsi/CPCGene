# -*- coding: utf-8 -*-
"""
Created on Thu May  3 14:33:51 2018

@author: rjovelin
"""

# use this script to generate missing files for the CPCG project

import argparse
import os
import sys
import subprocess
import copy


# Precondition bcftools module must be loaded prior
# module load bcftools/1.3.1
# module load samtools/1.5


# use this function to grab the normal and tumor bam for a single sample
def GetBamFiles(SampleFolder):
    '''
    (str) -> dict
    Take the full path of a sampe directory and return a dictionnary
    with normal and tumor bam if they exist
    Precondition: the sample folder contains 2 directories, 1 for tumor and 1 for normal    
    '''
    # create dict {'N': bam, 'T': bam}
    bams = {}
    
    # make a list of items 
    if os.path.isdir(SampleFolder) == True:
        L = os.listdir(SampleFolder)
        if len(L) != 0:
            # loop over normal and tumor subfolders, extract bams
            for i in range(len(L)):
                status = ''
                # check path
                subbamfolder = os.path.join(SampleFolder, L[i])
                assert os.path.isdir(subbamfolder) == True
                # get tumor/normal status
                if L[i][-2:] == 'B1' or L[i][-2:] == '_N': 
                    status = 'N'
                elif L[i][-2:] == 'F1' or L[i][-2:] == '_T':
                    status = 'T'
                if status != '':
                    if len(os.listdir(subbamfolder)) != 0:
                        for filename in os.listdir(subbamfolder):
                            if 'reheadered.bam' in filename:
                                if filename[filename.index('reheadered.bam'):] == 'reheadered.bam':
                                    assert os.path.isfile(os.path.join(subbamfolder, filename)) == True
                                    bams[status] = [os.path.join(subbamfolder, filename), filename]
                                    break
    return bams


# use this function to find the vcf of a given sample    
def FindSampleVCF(vcffolder):
    '''
    (str, str) -> str
    Take the full path of a sample directory and return the full path of the 
    vcf file if it exists
    '''
    # get the germline snv vcf file
    vcfFile = ''
    if os.path.isdir(vcffolder) == True:
        L = os.listdir(vcffolder)
        if len(L) != 0:
            for filename in L:
                if 'filtered_germline_snv_' in filename and '.vcf.gz' in filename:
                    if filename.startswith('filtered_germline_snv_') and filename[filename.index('.vcf.gz'):] == '.vcf.gz':
                        assert os.path.isfile(os.path.join(vcffolder, filename)) == True
                        vcfFile = os.path.join(vcffolder, filename)
                        break
    return vcfFile


# use this function to get the normal and tumor ID names of a given sample
def FindSampleNames(vcfFile):
    '''
    (str) -> dict
    Take the full path of a vcf file and return a dict with samples names
    contained in the vcf file header for each disease status
    '''
    # get the sample names
    assert os.path.isfile(vcfFile) == True
    mycomd = 'bcftools query -l ' + vcfFile 
    a = subprocess.check_output(mycomd, shell = True).decode('utf-8').split()
    sample_names = {}
    for i in a:
        if 'N_' in i or 'B1' in i or '_BD' in i or 'control' in i or i[-1] == 'b':
            sample_names['N'] = i
        elif 'T_' in i or 'F1' in i or '_PR' in i or 'tumor' in i or i[-1] == 'a':
            sample_names['T'] = i
    return sample_names

# use this function to look up normal and tumor ID for SU2C and TCGA samples
# normal/tumor status cannot easily be extracted from sample names
def GrabSampleNames(vcfFile, LookUp):
    '''
    (str, dict) -> dict
    Take the full path of a vcf file, a dictionary with sample : disease status
    and return a dict with sample names contained in the vcf file header for each disease status
    '''
        # get the sample names
    assert os.path.isfile(vcfFile) == True
    mycomd = 'bcftools query -l ' + vcfFile 
    a = subprocess.check_output(mycomd, shell = True).decode('utf-8').split()
    sample_names = {}
    for i in a:
        if i in LookUp:
            if LookUp[i] == 'normal':
                sample_names['N'] = i
            elif LookUp[i] == 'tumour':
                sample_names['T'] = i
    return sample_names

# use this function to extract disease information for SU2C dataset
def ParseClinical(Clinical):
    '''
    (file) -> dict
    Take a file with disease information for all samples and return a dictionary
    of sample: disease (normal/tumour) pairs
    '''
    DiseaseInfo = {}
    
    
    if 'phs000915' in Clinical:
        infile = open(Clinical)
        header = infile.readline().rstrip().split('\t')
        for line in infile:
            if line.startswith('WXS'):
                line = line.rstrip().split('\t')
                sample = line[header.index('Sample_Name')]
                if line[header.index('is_tumor')].upper() == 'YES':
                    disease = 'tumour'
                elif line[header.index('is_tumor')].upper() == 'NO':
                    disease = 'normal'
                assert sample not in DiseaseInfo
                DiseaseInfo[sample] = disease
        infile.close()
    return DiseaseInfo
    
# use this function to write a qsub script for individual contest
def WriteQsubContestInd(args, L, contaLevel, disease, Newfile):
    '''
    (list, list, str, str, file) -> None
    Take a list with command line arguments, a list of samples, a string specifying
    the contamination level, and the disease state, and write qsub commands to
    file NewFile open for writing    
    '''
    
    assert contaLevel in ['META', 'READGROUP']
    if contaLevel == 'META':
        fileExtension = '.per_bam.tsv'
    elif contaLevel == 'READGROUP':
        fileExtension = '.per_lane.tsv'
    DiseaseKey = ''
    if disease == 'normal':
        DiseaseKey = 'N'
    elif disease == 'tumour':
        DiseaseKey = 'T'
    
    # check if a look up table is provided to extract normal/tumour disease status
    if args.lookup:
        ### function to parse look up table
        LookUp = ParseClinical(args.lookup)
        
    # loop over sample lists
    if len(L) != 0:
        for sample in L:
            BamFolder = os.path.join(args.bamdir, sample)
            VcfFolder = os.path.join(args.vcfdir, sample)
            vcfFile = ''
            # check that folders exist            
            if os.path.isdir(BamFolder) and os.path.isdir(VcfFolder):
                # get the bams, vcf and normal, tumor sample names
                bams = GetBamFiles(BamFolder)
                vcfFile = FindSampleVCF(VcfFolder)
                if 'SU2C' in sample:
                    sample_names =  GrabSampleNames(vcfFile, LookUp)
                else:
                    sample_names = FindSampleNames(vcfFile)
            if DiseaseKey in bams and vcfFile != '' and DiseaseKey in sample_names:
                output_file = 'contest_' + sample_names[DiseaseKey] + fileExtension
                line = ['qsub -cwd -b y -N contest_{0} -l h_vmem={1}g'.format(contaLevel + '.' + bams[DiseaseKey][1], args.mem),
                        '\"module load Perl-BL; module load gatk/{0};'.format(args.gatk),
                        'java -Xmx{0}g'.format(args.java), 
                        '-jar \$GATKROOT/GenomeAnalysisTK.jar',
                        '-T ContEst -R /oicr/data/genomes/homo_sapiens_mc/1000g/phase2_reference_assembly/bwa/0.7.15/hs37d5.fa',
                        '-sn {0}'.format(sample_names[DiseaseKey]),
                        '--popfile /.mounts/labs/boutroslab/private/Resources/Sequencing/gatkBundle/2.8/b37/b37_population_stratified_af_hapmap_3.3.vcf.gz',
                        '--genotypes {0}'.format(vcfFile),
                        '-isr INTERSECTION',
                        '--lane_level_contamination {0}'.format(contaLevel),
                        '-I {0}'.format(bams[DiseaseKey][0]),
                        '-o {0}\"'.format(os.path.join(args.outdir, output_file))]
                Newfile.write(' '.join(line) + '\n')
            else:
                print('missing entry for {0} (bams:{1}, vcf:{2}, names:{3})'.format(sample, len(bams), vcfFile, len(sample_names)))


# use this function to run ContEst on individual samples
def RunContestInd(args):
    '''
    (list) -> file
    Tale the list of command line arguments and write a file with qsub commands
    to run ContEst on individual samples (as opposed to paired samples) 
    
    '''
    
    # get the file with missing info
    infile = open(args.status)
    header = infile.readline().rstrip().split()
    # create a list with sample names for which normal or tumor is missing contest files
    NB, NL, TB, TL = [], [], [] , []
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            name = line[0]
            # check file status, check also that a germline snv file is present
            if int(line[header.index('filtered_germline_snv:vcf.gz')]) == 0:
                print('missing germline snv calls for {0}'.format(name))
            else:
                if int(line[header.index('N:per_bam.tsv')]) == 0:
                    NB.append(name)
                if int(line[header.index('N:per_lane.tsv')]) == 0:
                    NL.append(name)
                if int(line[header.index('T:per_bam.tsv')]) == 0:
                    TB.append(name)
                if int(line[header.index('T:per_lane.tsv')]) == 0:
                    TL.append(name)
    infile.close()

    ### 2) write command to file for each sample for tumor and normal and for readgroup and meta level analyses
    newfile = open(args.qsub, 'w')
    
    # write commands for normal and tumor samples for META and READGROUP levels
    WriteQsubContestInd(args, NB, 'META', 'normal', newfile)
    WriteQsubContestInd(args, NL, 'READGROUP', 'normal', newfile)
    WriteQsubContestInd(args, TB, 'META', 'tumour', newfile)
    WriteQsubContestInd(args, TL, 'READGROUP', 'tumour', newfile)
        
    newfile.close()


# use this function to run ContEst on norma/tumor sample pairs    
def RunContestPairs(args):
    '''
    (list) -> file
    Take the list of command line arguments and write qsubs to outputfile
    to run ContEst on normal/tumor sample pairs
    '''
        
    # get the file with missing info
    infile = open(args.status) 
    header = infile.readline().rstrip().split()
    # create a list with sample names for which normal or tumor is missing contest files
    missing = []
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            name = line[0]
            assert name not in missing
            # check file status, check also that a germline snv file is present
            if 0 in [int(line[header.index('N:per_bam.tsv')]),
                     int(line[header.index('N:per_lane.tsv')]),
                     int(line[header.index('T:per_bam.tsv')]),
                     int(line[header.index('T:per_lane.tsv')])]:
                if int(line[header.index('filtered_germline_snv:vcf.gz')]) == 0:
                    print('missing germline snv calls for {0}'.format(name))
                else:
                    missing.append(name)
    infile.close()

    ### 2) write command to file for each sample for tumor and normal and for readgroup and meta level analyses
    newfile = open(args.qsub, 'w')

    for sample in missing:
        # get the bams
        bams = {}
        bamfolder = os.path.join(args.bamdir, sample)
        vcffolder = os.path.join(args.vcf, sample)    
        vcfFile = ''
        # check that folders exist
        if os.path.isdir(bamfolder) and os.path.isdir(vcffolder):
            vcfFile = FindSampleVCF(vcffolder)
            bams = GetBamFiles(bamfolder)
            sample_names = FindSampleNames(vcfFile)
            
        # check that bams, vcf, and samples exists
        if len(bams) != 0 and vcfFile != '' and len(sample_names) != 0:
            # run contest for for lane and meta for both normal and tumor
            for i in bams:
                line1 = ['ContEst.pl',
                         '--bam {0}'.format(bams[i][0]),
                         '--sample {0}'.format(sample_names[i]),
                         '--genotype {0}'.format(vcfFile),
                         '--analysis_level READGROUP',
                         '--ref /oicr/data/genomes/homo_sapiens_mc/1000g/phase2_reference_assembly/bwa/0.7.15/hs37d5.fa',
                         '--ref_type b37',                    
                         '--dir {0}'.format(args.outdir),
                         '--gatk_version {0}'.format(args.gatk),
                         '--ext .per_lane']
                line2 = ['ContEst.pl',
                         '--bam {0}'.format(bams[i][0]),
                         '--sample {0}'.format(sample_names[i]),
                         '--genotype {0}'.format(vcfFile),
                         '--analysis_level META',
                         '--ref /oicr/data/genomes/homo_sapiens_mc/1000g/phase2_reference_assembly/bwa/0.7.15/hs37d5.fa',
                         '--ref_type b37',
                         '--dir {0}'.format(args.outdir),
                         '--gatk_version {0}'.format(args.gatk),
                         '--ext .per_bam']
                line1 = ' '.join(line1)
                line2 = ' '.join(line2)
                BaseCom1 =  'qsub -cwd -b y -N {0} -e {1} -o {2} -l h_vmem=10g '.format(sample + '_L_' + i, os.path.join(args.outdir, 'logs'), os.path.join(args.outdir, 'logs'))
                BaseCom2 =  'qsub -cwd -b y -N {0} -e {1} -o {2} -l h_vmem=10g '.format(sample + '_B_' + i, os.path.join(args.outdir, 'logs'), os.path.join(args.outdir, 'logs'))
                newfile.write(BaseCom1 + '\"module load Perl-BL; ' + line1 + '\"' + '\n')
                newfile.write(BaseCom2 + '\"module load Perl-BL; ' + line2 + '\"' + '\n')
        else:
            print('missing entry for {0} (bams:{1}, vcf:{2}, names{3})'.format(sample, len(bams), len(vcfFile), len(sample_names)))
    
    newfile.close()


# use this function to run callable_bases through Gatk pipeline on pair of tumor/normal bams    
def RunGATK(args):
    '''
    (list) -> file
    Take the list of arguments from the command line and write a qsub script to run
    GATK on pair of normal/tumor sample
    '''
    ### 1) get the file with missing info
    infile = open(args.status)
    header = infile.readline().rstrip().split()
    # create a list with sample names missing callable-bases files
    missing = []
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            name = line[0]
            assert name not in missing
            if args.callable == 'Y':
                # check file status
                if 0 in [int(line[header.index('.tar.gz')]), int(line[header.index('CallableBases_exact')])]:
                    # check if bams are present
                    if int(line[header.index('N:reheadered.bam')]) == 0 or int(line[header.index('T:reheadered.bam')]) == 0:
                        print('missing bam for {0}'.format(name))
                    else:
                        missing.append(name)
            if args.snp == 'Y':
                # check vcf file status
                if 0 in [int(line[header.index('merged_raw_recalibrated.vcf.gz')]),
                         int(line[header.index('merged_raw_recalibrated.vcf.gz.tbi')]),
                         int(line[header.index('filtered_germline_indel:vcf.gz')]),
                         int(line[header.index('filtered_germline_indel:vcf.gz.tbi')]),
                         int(line[header.index('filtered_germline_snv:vcf.gz')]),
                         int(line[header.index('filtered_germline_snv:vcf.gz.tbi')])]:
                    if int(line[header.index('N:reheadered.bam')]) == 0 or int(line[header.index('T:reheadered.bam')]) == 0:
                        print('missing bam for {0}'.format(name))
                    else:
                        missing.append(name)

    missing = list(set(missing))
    infile.close()

    ### 2) write command to file for each sample with missing files; run tumor and normal together
    newfile = open(args.qsub, 'w')

    for sample in missing:
        # reinitialize variables for each new sample
        sampleName, tumorID, normalID = '', '', ''
        # get the bams
        bams = {}
        bamfolder = os.path.join(args.bamdir, sample)
        # get the normal and tumor bams
        bams = GetBamFiles(bamfolder)
        # check that bams for normal and tumor are present    
        if len(bams) == 2:
            # get sample name, tumorID and normalID
            sampleName, tumorID, normalID = '', '', '' 
            if 'Baca' in sample:
                sampleName = sample[sample.index('.')+1:]
                tumorID = 'T_' + sample[sample.index('.')+1:]
                normalID = 'N_' + sample[sample.index('.')+1:]
            elif 'Weischenfeldt' in sample:
                sampleName = sample[sample.index('.')+1:]        
                tumorID = bams['T'][1][bams['T'][1].index('Weischenfeldt')+len('Weischenfeldt')+1:bams['T'][1].index('_realigned_')]
                normalID = bams['N'][1][bams['N'][1].index('Weischenfeldt')+len('Weischenfeldt')+1:bams['N'][1].index('_realigned_')]
            elif 'PRAD-' in sample:
                sampleName = sample[sample.index('.')+1:]
                tumorID = bams['T'][1][:bams['T'][1].index('_realigned_')]
                normalID = bams['N'][1][:bams['N'][1].index('_realigned_')]
            elif 'Berger' in sample or 'TCGA-' in sample:
                sampleName = sample
                tumorID = bams['T'][1][:bams['T'][1].index('_realigned_')]
                normalID = bams['N'][1][:bams['N'][1].index('_realigned_')]
            elif 'CPCG' in sample:
                if '-' not in sample:
                    sampleName = sample
                    tumorID = bams['T'][1][:bams['T'][1].index('_realigned_')]
                    normalID = bams['N'][1][:bams['N'][1].index('_realigned_')]
                elif '-' in sample:
                    sampleName = sample
                    tumorID = bams['T'][1][:bams['T'][1].index('.merged')]
                    normalID = bams['N'][1][:bams['N'][1].index('.merged')]
            assert tumorID != '' and normalID != ''
        
        # check that bams and samples exist
        if len(bams) == 2 and sampleName != '' and tumorID != '' and normalID != '':
            # build command line
            line = ['GATK.pl',
                    '--ref /oicr/data/genomes/homo_sapiens_mc/1000g/phase2_reference_assembly/bwa/0.7.15/hs37d5.fa',
                    '--ref_type b37',
                    '--sample {0}'.format(sampleName),
                    '--normal {0}'.format(bams['N'][0]),
                    '--tumour {0}'.format(bams['T'][0]),
                    '--normalID {0}'.format(normalID),
                    '--tumourID {0}'.format(tumorID),
                    '--GATK_processed {0}'.format(args.gatk_processed),
                    '--Pre_processing {0}'.format(args.pre_processing),
                    '--callable_bases {0}'.format(args.callable),
                    '--SNP_calling {0}'.format(args.snp),
                    '--gatk_version {0}'.format(args.gatk),
                    '--output_dir {0}'.format(args.outdir)]
            line = ' '.join(line)
            BaseCom =  'qsub -cwd -b y -N {0} -e {1} -o {2} -l h_vmem={3}g '.format('gatk_' + sampleName, os.path.join(args.outdir, 'logs'), os.path.join(args.outdir, 'logs'), args.mem)
            newfile.write(BaseCom + '\"module load {0}; '.format(args.pbl) + line + '\"' + '\n')
        else:
            print('missing entry for {0} (smaple_name: {1}, normalID: {2}, tumorID: {3}, bams: {4})'.format(sample, sampleName, normalID, tumorID, len(bams)))
    
    newfile.close()

# use this function to run flagstat on missing samples
def RunFlagstat(args):
    '''
    (list) -> file
    Take the list of command line arguments and write qsubs to outputfile to run
    flagstat on indivial tumor or normal bams
    '''
    
    # get the file with missing info
    infile = open(args.status)
    header = infile.readline().rstrip().split()
    # create lists with sample names missing flagstat files for normal or tumor
    missingNormal, missingTumor = [], []
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            name = line[0]
            # check file status
            if int(line[header.index('N:bam.flagstat.txt')]) == 0:
                # check that bam is present
                if int(line[header.index('N:reheadered.bam')]) == 0:
                    print('missing normal bam for {0}'.format(name))
                else:
                    missingNormal.append(name)
            if int(line[header.index('T:bam.flagstat.txt')]) == 0:
                # check that bam is present
                if int(header.index('T:reheadered.bam')) == 0:
                    print('missing tumor bam for {0}'.format(name))
                else:
                    missingTumor.append(name)
    infile.close()

    ### 2) write command to file for each sample with missing files; run tumor and normal together
    newfile = open(args.qsub, 'w')

    # make a common list of missing samples
    MissingSamples = missingNormal + missingTumor
    MissingSamples = list(set(MissingSamples))
    # create a dict of bams per sample {sample: {'N': [path+bam, bam_name], 'T': [path+bam, nam_name]}}
    SampleBams = {}
    # grab the normal and tumor bams for each sample
    for sample in MissingSamples:
        SampleBams[sample] = {}
        bamfolder = os.path.join(args.bamdir, sample)
        # get bams for the given sample
        bams = GetBamFiles(bamfolder)
        for status in bams:
            SampleBams[sample][status] = bams[status]
        
    # loop over samples missing normal flagstat
    for sample in missingNormal:
        # check that sample has bam
        if sample in SampleBams:
            if 'N' in SampleBams[sample]:
                line = 'samtools flagstat {0} > {1}'.format(SampleBams[sample]['N'][0], os.path.join(args.outdir, SampleBams[sample]['N'][1] + '.flagstat.txt'))
                BaseCom =  'qsub -cwd -b y -N {0} -e {1} -o {2} -l h_vmem=10g '.format('Flag_' + sample + '_N', os.path.join(args.outdir, 'logs'), os.path.join(args.outdir, 'logs'))
                newfile.write(BaseCom + '\"module load samtools; ' + line + '\"' + '\n')
            else:
                print('missing normal bam for {0}'.format(sample))
    for sample in missingTumor:
        # check that sample has bam
        if sample in SampleBams:
            if 'T' in SampleBams[sample]:
                line = 'samtools flagstat {0} > {1}'.format(SampleBams[sample]['T'][0], os.path.join(args.outdir, SampleBams[sample]['T'][1] + '.flagstat.txt'))
                BaseCom =  'qsub -cwd -b y -N {0} -e {1} -o {2} -l h_vmem=10g '.format('Flag_' + sample + '_T', os.path.join(args.outdir, 'logs'), os.path.join(args.outdir, 'logs'))
                newfile.write(BaseCom + '\"module load samtools; ' + line + '\"' + '\n')
            else:
                print('missing tumour bam for {0}'.format(sample))
    newfile.close()


# use this script to run coverage stats on normal and tumor bams separately
def RunCoverageStats(args):
    '''
    (list) -> file
    Take the list of command line arguments and write qsubs to outputfile to run
    coverage stats on individual bams
    '''

    # get the file with missing info
    infile = open(args.status)
    header = infile.readline().rstrip().split()
    # create lists with sample names missing flagstat files for normal or tumor
    missingNormal, missingTumor = [], []
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            name = line[0]
            # check file status
            NormKeys = ['N:library_statistics', 'N:library_summary', 'N:read_group_statistics',
                       'N:read_group_summary', 'N:sample_statistics', 'N:sample_summary']
            TumKeys = ['T:library_statistics', 'T:library_summary', 'T:read_group_statistics',
                       'T:read_group_summary', 'T:sample_statistics', 'T:sample_summary']
            L1 = [line[header.index(NormKeys[i])] for i in range(len(NormKeys))]
            L2 = [line[header.index(TumKeys[i])] for i in range(len(TumKeys))]
            if 0 in list(map(lambda x: int(x), L1)):
                # check if normal bams are present
                if int(line[header.index('N:reheadered.bam')]) == 0:
                    print('missing normal bam for {0}'.format(name))
                else:
                    missingNormal.append(name)
            if 0 in list(map(lambda x: int(x), L2)):
                # check if tumor bams are present
                if int(line[header.index('T:reheadered.bam')]) == 0:
                    print('missing tumor bam for {0}'.format(name))
                else:
                    missingTumor.append(name)
    infile.close()

    ### 2) write command to file for each sample with missing files; run tumor and normal together
    newfile = open(args.qsub, 'w')

    # make a common list of missing samples
    MissingSamples = missingNormal + missingTumor
    MissingSamples = list(set(MissingSamples))
    # create a dict of bams per sample {sample: {'N': [path+bam, sample], 'T': [path+bam, sample]}}
    SampleBams = {}

    # grab the normal and tumor bams for each sample
    for sample in MissingSamples:
        SampleBams[sample] = {}
        bamfolder = os.path.join(args.outdir, sample)
        # grab bams 
        bams = GetBamFiles(bamfolder)
        for status in bams:
            SampleBams[sample][status] = bams[status]
        # add sample to list
        for sample in SampleBams:
            for status in SampleBams[sample]:
                filename = SampleBams[sample][status][1]
                sampleName = ''
                if '.merged' in filename:
                    sampleName = filename[:filename.index('.merged')]
                elif '_realigned_' in filename:
                    sampleName = filename[:filename.index('_realigned_')]
                assert sampleName != ''
                SampleBams[sample][status].append(sampleName)
          
    # loop over samples missing normal flagstat
    for sample in missingNormal:
        # check if sample has normal bam
        if sample in SampleBams:
            if 'N' in SampleBams[sample]:
                line = ['module load gatk/{0};'.format(args.gatk),
                        'java -Xmx10g -jar \$GATKROOT/GenomeAnalysisTK.jar',
                        '-T DepthOfCoverage',
                        '-R /oicr/data/genomes/homo_sapiens_mc/1000g/phase2_reference_assembly/bwa/0.7.15/hs37d5.fa',
                        '-o {0}'.format(os.path.join(args.outdir, SampleBams[sample]['N'][2])),
                        '-I {0}'.format(SampleBams[sample]['N'][0]),
                        '-omitBaseOutput -omitIntervals -omitLocusTable -nt 2 -pt sample -pt readgroup -pt library']
                line =  ' '.join(line)             
                BaseCom =  'qsub -cwd -b y -N {0} -e {1} -o {2} -l h_vmem=16g '.format('Cover_' + SampleBams[sample]['N'][1], os.path.join(args.outdir, 'logs'), os.path.join(args.outdir, 'logs'))
                newfile.write(BaseCom + '\"' + line + '\"' + '\n')
    for sample in missingTumor:
        # check if sample has tumor bam
        if sample in SampleBams:
            if 'T' in SampleBams[sample]:
                line = ['module load gatk/{0};'.format(args.gatk),
                        'java -Xmx10g -jar \$GATKROOT/GenomeAnalysisTK.jar',
                        '-T DepthOfCoverage',
                        '-R /oicr/data/genomes/homo_sapiens_mc/1000g/phase2_reference_assembly/bwa/0.7.15/hs37d5.fa',
                        '-o {0}'.format(os.path.join(args.outdir, SampleBams[sample]['T'][2])),
                        '-I {0}'.format(SampleBams[sample]['T'][0]),
                        '-omitBaseOutput -omitIntervals -omitLocusTable -nt 2 -pt sample -pt readgroup -pt library']
                line = ' '.join(line)             
                BaseCom =  'qsub -cwd -b y -N {0} -e {1} -o {2} -l h_vmem=16g '.format('Cov_' + sample + '_T', os.path.join(args.outdir, 'logs'), os.path.join(args.outdir, 'logs'))
                newfile.write(BaseCom + '\"' + line + '\"' + '\n')
    
    newfile.close()


  
# use this script to run filter_GATK_SNV_calls.pl to get germline genotype counts
def RunGermlineCount(args):
    '''
    (list) -> file
    Take the list of command line arguments and write qsubs to outputfile
    to run filter_GATK_SNV_calls.pl and generate germ line genotype counts on VCFs
    '''




    if args.status:
        # get the file with missing info
        infile = open(args.status)
        header = infile.readline().rstrip().split()
        # create a list with sample names for which normal or tumor is missing germline counts files
        missing = []
        for line in infile:
            if line.rstrip() != '':
                line = line.rstrip().split()
                name = line[0]
                assert name not in missing
                # check file status, check also that a germline snv file is present
                if 0 in [int(line[header.index('filtered_germline_genotype_count:tsv')]),
                         int(line[header.index('filtered_germline_variant_class_count:tsv')])]:
                    if int(line[header.index('merged_raw_recalibrated.vcf.gz')]) == 0:
                        print('missing germline snv calls for {0}'.format(name))
                    else:
                        missing.append(name)
        infile.close()
    elif args.samples:
        missing = []
        # get the list of samples
        infile = open(args.samples)
        for line in infile:
            line = line.rstrip()
            if line != '':
                missing.append(line)
        infile.close()

    ### 2) write command to file for each sample for tumor and normal and for readgroup and meta level analyses

    newfile = open(args.qsub, 'w')

    for sample in missing:
        # get the germline snv vcf file
        vcfFile = ''
        vcffolder = os.path.join(args.vcfdir, sample)
        assert os.path.isdir(vcffolder) == True
        L = os.listdir(vcffolder)
        if len(L) != 0:
            for filename in L:
                if 'merged_raw_' in filename and 'recalibrated.vcf.gz' in filename and 'tbi' not in filename:
                    assert sample in filename
                    assert os.path.isfile(os.path.join(vcffolder, filename)) == True
                    vcfFile = os.path.join(vcffolder, filename)
                    break

        # get the sample names
        assert os.path.isfile(os.path.join(vcffolder, vcfFile)) == True
        mycomd = 'bcftools query -l ' + os.path.join(vcffolder, vcfFile) 
        a = subprocess.check_output(mycomd, shell = True).decode('utf-8').split()
        sample_names = {}
        for i in a:
            if 'N_' in i or 'B1' in i or '_BD' in i or 'control' in i:
                sample_names['N'] = i
            elif 'T_' in i or 'F1' in i or '_PR' in i or 'tumor' in i:
                sample_names['T'] = i
        print(sample, sample_names, sep = '\t')
    
        # get the sample names
        if os.path.isfile(os.path.join(vcffolder, vcfFile)) == True and len(sample_names) == 2:
            line = ['filter_GATK_SNV_calls.pl',
                    '--input {0}'.format(vcfFile),
                    '--sample {0}'.format(sample),
                    '--ref /oicr/data/genomes/homo_sapiens_mc/1000g/phase2_reference_assembly/bwa/0.7.15/hs37d5.fa',
                    '--output_dir {0}'.format(args.outdir),
                    '--normal {0}'.format(sample_names['N']),
                    '--tumour {0}'.format(sample_names['T']),
                    '--filter_somatic Y',
                    '--filter_ambiguous Y',
                    '--split_calls Y']
            line = ' '.join(line)
            BaseCom =  'qsub -cwd -b y -N {0} -l h_vmem=14g '.format('germline_' + sample)
            newfile.write(BaseCom + '\"module load Perl-BL; ' + line + '\"' + '\n')
        else:
            print('missing entry for {0} (vcf:{1})'.format(sample, vcfFile))
    newfile.close()


if __name__ == '__main__':
    
    # create top-level parser
    parser = argparse.ArgumentParser(prog = 'GenerateMissingCPCG.py', description='generates missing files post dnaalign and gatk pipelines')
    subparsers = parser.add_subparsers(title='sub-commands', description='valid sub-commands', help = 'sub-commands help')
    
    # Callable_bases sub-commands
    GATK_parser = subparsers.add_parser('GATK', help ='Run phaseI (recalibrated_reheadered bams), phaseII (coverage or callable bases) or phaseIII (variant calling) analyses through GATK')
    GATK_parser.add_argument('-s', '--Status', dest='status', help='summary file with CPCG file status', required=True)
    GATK_parser.add_argument('-b', '--Bams', dest='bamdir', help='direcory with sample directories containing bams', required=True)
    GATK_parser.add_argument('-o', '--Out', dest='outdir', help='output directory where files are written', required=True)
    GATK_parser.add_argument('-g', '--Gatk', dest='gatk', help='gatk version', choices=['3.4.0', '3.5.0', '3.6.0', '3.7.0'], required=True)
    GATK_parser.add_argument('-p', '--PBL', dest='pbl', default='Perl-BL', help='Perl-BL version. eg: Per-BL/2017-12-01. default Perl-BL nightly build')
    GATK_parser.add_argument('-q', '--Qsub', dest='qsub', help='qsub script outputfile', required=True)
    GATK_parser.add_argument('-m', '--Mem', dest='mem', default='10', help='SGE memory. 10g by default')
    GATK_parser.add_argument('--GATK_processed', dest='gatk_processed', default='Y', choices =['Y', 'N'], help='if GATK pre-processing is needed. default = Y, assumes already gatk processed')
    GATK_parser.add_argument('--Pre_processing', dest='pre_processing', default='N', choices=['Y', 'N'], help='perform GATK pre-processing. default= N, assumes pre-processing is done')
    GATK_parser.add_argument('--callable_bases', dest='callable', choices=['Y','N'], help='perform callable bases through GATK')
    GATK_parser.add_argument('--SNP_calling', dest='snp', default = 'N', choices=['N', 'Y'], help = 'perform variant calling and generate the vcfs. default = N')
    GATK_parser.set_defaults(func=RunGATK)
 
    # Contest on individual samples sub-commands
    ContestInd_parser = subparsers.add_parser('contestInd', help='Run ContEst on individual bams')               
    ContestInd_parser.add_argument('-s', '--Status', dest='status', help='summary file with CPCG file status', required=True)
    ContestInd_parser.add_argument('-b', '--Bams', dest='bamdir', help='direcory with sample directories containing bams', required=True)
    ContestInd_parser.add_argument('-v', '--Vcf', dest='vcfdir', help='directory with sample directories containing vcf', required=True)    
    ContestInd_parser.add_argument('-o', '--Out', dest='outdir', help='output directory where files are written', required=True)
    ContestInd_parser.add_argument('-g', '--Gatk', dest='gatk', help='gatk version', choices=['3.4.0', '3.5.0', '3.6.0', '3.7.0'], required=True)
    ContestInd_parser.add_argument('-m', '--Mem', dest='mem', default='20', help='SGE memory. 20g by default')
    ContestInd_parser.add_argument('-j', '--java', dest='java', default='16', help='java memory. 16g by default')
    ContestInd_parser.add_argument('-q', '--Qsub', dest='qsub', help='qsub script outputfile', required=True)
    ContestInd_parser.add_argument('-l', '--Lookup', dest='lookup', help='file with normal/tumor info')
    ContestInd_parser.set_defaults(func=RunContestInd)
  
    # Contest on normal/tumour pairs sub-commands
    ContestPair_parser = subparsers.add_parser('contestPairs', help ='Run ContEst on normal/tumor sample pairs')
    ContestPair_parser.add_argument('-s', '--Status', dest='status', help='summary file with CPCG file status', required=True)
    ContestPair_parser.add_argument('-b', '--Bams', dest='bamdir', help='direcory with sample directories containing bams', required=True)
    ContestPair_parser.add_argument('-v', '--Vcf', dest='vcfdir', help='directory with sample directories containing vcf', required=True)
    ContestPair_parser.add_argument('-o', '--Out', dest='outdir', help='output directory where files are written', required=True)
    ContestPair_parser.add_argument('-g', '--Gatk', dest='gatk', help='gatk version', choices=['3.4.0', '3.5.0', '3.6.0', '3.7.0'], required=True)
    ContestPair_parser.add_argument('-q', '--Qsub', dest='qsub', help='qsub script outputfile', required=True)
    ContestPair_parser.set_defaults(func=RunContestPairs)

    # Flagstat sub-commands
    Flagstat_parser = subparsers.add_parser('flagstat', help='Run flagstat on individual bams')    
    Flagstat_parser.add_argument('-s', '--Status', dest='status', help='summary file with CPCG file status', required=True)
    Flagstat_parser.add_argument('-b', '--Bams', dest='bamdir', help='direcory with sample directories containing bams', required=True)
    Flagstat_parser.add_argument('-o', '--Out', dest='outdir', help='output directory where files are written', required=True)
    Flagstat_parser.add_argument('-q', '--Qsub', dest='qsub', help='qsub script outputfile', required=True)
    Flagstat_parser.set_defaults(func=RunFlagstat)
        
    # Coverage stats sub-commands
    CoverageStats_parser = subparsers.add_parser('coverage', help='Run coverage stats on individual bams')
    CoverageStats_parser.add_argument('-s', '--Status', dest='status', help='summary file with CPCG file status', required=True)
    CoverageStats_parser.add_argument('-b', '--Bams', dest='bamdir', help='direcory with sample directories containing bams', required=True)
    CoverageStats_parser.add_argument('-o', '--Out', dest='outdir', help='output directory where files are written', required=True)
    CoverageStats_parser.add_argument('-g', '--Gatk', dest='gatk', help='gatk version', choices=['3.4.0', '3.5.0', '3.6.0', '3.7.0'], required=True)
    CoverageStats_parser.add_argument('-q', '--Qsub', dest='qsub', help='qsub script outputfile', required=True)
    CoverageStats_parser.set_defaults(func=RunCoverageStats)

    # Germline count subcommands
    GermlineCount_parser = subparsers.add_parser('germlineCounts', help='run filter_GATK_SNV_calls to get the germline genotype counts')    
    GermlineCount_parser.add_argument('-s', '--Status', dest='status', help='summary file with CPCG file status')
    GermlineCount_parser.add_argument('-sp', '--Samples', dest='samples', help='file with sample names')
    GermlineCount_parser.add_argument('-v', '--Vcf', dest='vcfdir', help='directory with sample directories containing vcf', required=True)
    GermlineCount_parser.add_argument('-o', '--Out', dest='outdir', help='output directory where files are written', required=True)
    GermlineCount_parser.add_argument('-q', '--Qsub', dest='qsub', help='qsub script outputfile', required=True)
    GermlineCount_parser.set_defaults(func=RunGermlineCount)
    
    # get arguments from the command line
    args = parser.parse_args()
    # pass the args to the default function
    args.func(args)

