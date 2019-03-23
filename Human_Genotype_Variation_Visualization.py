""" This Python Code implements Visualization and Dimensionality Reduction of Genotype Data from the 
    Human Evolution and Genetic History Data based on the 1000 GENOMES PROJECT """

import time
extraction_start = time.time()
from collections import Counter
import gzip
import numpy as np

pca_n_components = [5, 15, 25]   # Combinations to iterate for optimalization of PCA Dimensionality Reduction
umap_neighbors = [5, 10, 15]     # UMAP iterates for Number of Neighbors and Minimum Distance
umap_min_dist = [0.1, 0.25, 0.5]

combined = [(a,c,d) for a in pca_n_components for c in umap_neighbors for d in umap_min_dist]  # List Comprehension for Combinations

genotype_file = 'ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz'  # File names for Data, VCF file for Genotype Data
population_file = '20131219.populations.tsv'                                      # Tsv file with Population Gene Classification
panel_file = 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/affy_samples.20141118.panel'


class snp(object):

    def __init__(self, line, select=False, autosome_only =True):
        """The initialization method takes in a line from the vcf file, as a string, 
        and records the relevant information. 
        line: a string from a vcf file
        select: a list of positions of individuals to be analyzed, where positions run from 0 to 
        nInd-1, the number of individuals
        """ 
        
        split_line = line.split()  #  First break down the line into a list of each field
        
        self.failed = False  # A label that we will set to True if something goes wrong.
        
        if line.startswith('#'):
            self.failed = True
            self.failure_cause = "line was a header line, not a snp"
            return
        
        if len(split_line)<=5:
            self.failed = True
            self.failure_cause = "incorrectly formatted line, should have at least 5 fields " + line
            return
          
        self.chrom = split_line[0]
        if autosome_only:
            if self.chrom not in ["%d" % (i,) for i in range(1,23)]:
                self.failed = True
                self.failure_cause = "not recognized as an autosome while autosome_only set to True"
                return
        
        self.chrom = int(split_line[0]) # Chromosome (numbered)
        self.position = int(split_line[1])  # The coordinates of the snp
        self.rid = split_line[2] # Name/Record ID
        self.ref_allele = split_line[3]
        self.alt_allele = split_line[4] # The alterate allele according to the vcf; also a string 
        # Only accept snps in ACGT. 
        if self.ref_allele not in ["A","C","G","T"] or self.alt_allele not in ["A","C","G","T"]:
            self.failed = True
            self.failure_cause = "ref or alt not in ACGT"
            return
        self.filter = split_line[6]  # See vcf format specifications for the interpretation of 
                                    # the filter field
        if self.filter not in ['PASS', '.'] :  # PASS indicates a SNP that passed all QC filters.
            self.failed = True
            self.failure_cause = self.filter
            return
              
        self.genotype_strings = split_line[9:]

        # Prepare a list that will contain the transformed genotypes. 
        # Since we already know how long the list will be, it makes sense 
        # to create an array of zeros of the same length as self.gtypes, 
        
        self.genotype_array = np.zeros(len(self.genotype_strings), dtype = np.int8)             

        # Count the number of each genotype. 
        # There may be different strings giving the same genotype so we increment the 
        # counts found so far for the genotype by the number of times the  
        # For example, "0/0" and "0\0" give homref, and "0|1" and "1|0" give het
        
        n_missing = 0
        for index,genotype_string in enumerate(self.genotype_strings):
            if genotype_string == './.':
                n_missing +=1 
                self.genotype_array[index]=-1
                continue # missing data will be left as 0
            allele_0 = genotype_string[0] # Get the first allele (as a string)
            allele_1 = genotype_string[2]
            if (allele_0=='1' and allele_1=='1'): # Use rstrip because windows machines will occasionally have extra \n
                self.genotype_array[index]=2
            elif ((allele_0=='0' and allele_1=='1') or (allele_0=='1' and allele_1=='0')):
                self.genotype_array[index]=1   
            elif (allele_0=='0' and allele_1=='0'):
                # The array was initialized to zero, so nothing to do here!
                continue
            else:
                print(("unknown genotype", genotype_string))
                self.failed=True
                self.failedreason="unknown genotype"
                return

# Specify the number of lines to skip to avoid storing every line in memory
number_of_lines_to_skip = 20


genotype_matrix = []  # Will contain our numerical genotype matrix. 
genotype_positions = []
genotype_names = []
x = 0
error_count = 0

with gzip.open(genotype_file,'rt') as f:
    count = 0
    for line in f:
        count+=1
        if count % number_of_lines_to_skip == 0:
            if line.startswith("#") or snp(line).failed:
                if snp(line).failure_cause != "line was a header line, not a snp":
                    error_count += 1
                    if x < 10:
                        print('Failed: ' + snp(line).failure_cause)
                        x+=1
                continue
            
            return_snp = snp(line)
            genotype_matrix.append(return_snp.genotype_array)
            genotype_names.append(return_snp.rid)
            genotype_positions.append([return_snp.chrom, return_snp.position])

extraction_end = time.time()

import logging
import sys


def log_file(logger_name, level=logging.INFO):
    """
    Method to return a custom logger with the given name and level
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(level)
    format_string = "%(levelname)s %(asctime)s - %(message)s"
    log_format = logging.Formatter(format_string)

    # Creating and adding the file handler
    file_handler = logging.FileHandler(logger_name, mode='w')
    file_handler.setFormatter(log_format)
    logger.addHandler(file_handler)
    return logger



for x in range(0,len(combined)):    # Iterating over possible combinations
    
    start = time.time()
    a = str(combined[x][0])
    b = str(combined[x][1])
    c = str(combined[x][2])
    
    #LOG_FORMAT = '%(levelname)s %(asctime)s - %(message)s'
    #filename = 'Log_Output_'+a+'_'+b+'_'+c+'.log'
    logger = log_file('Log_Output_'+a+'_'+b+'_'+c+'.log')
    #logger.basicConfig(filename = 'Log_Output_'+a+'_'+b+'_'+c+'.log', level = logger.INFO , format = LOG_FORMAT, filemode = 'w')
    logger.info('Logger Initiated')
    logger.info("Genotype Extraction from VCF file time consumed : {} seconds".format(extraction_end-extraction_start))
    logger.info("Genotype data lines skipped : {}".format(number_of_lines_to_skip))
    logger.info("Importing all the necessary Packages")

    import pandas as pd
    import sys
    import numpy as np
    import allel
    import seaborn as sns
    import matplotlib
    import matplotlib.pyplot as plt 
    import sklearn
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    import umap
    import gzip
    import warnings
    warnings.filterwarnings('ignore')


    logger.info("Packages Imported with Specified Versions")
    logger.info("Pandas - {} " .format(pd.__version__))
    logger.info("Scikit-Allel - {} " .format(allel.__version__))
    logger.info("Seaborn - {} " .format(sns.__version__))
    logger.info("Matplotlib - {} " .format(matplotlib.__version__))
    logger.info("Scikit-Learn - {} " .format(sklearn.__version__))
    logger.info("UMAP - {} " .format(umap.__version__))
    logger.info("Numpy - {} " .format(np.__version__))
    logger.info("Python - {} " .format(sys.version.split('|')[0]))

    logger.info("IMPORT COMPLETED")

    #genotype_file = 'ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz'
    #population_file = '20131219.populations.tsv'
    #panel_file = 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/affy_samples.20141118.panel'

    

    logger.info("Extraction Completed with data size {} rows".format(len(genotype_matrix)))
    logger.info("Scaling Data using Standard Scaler")

    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(genotype_matrix)

    logger.info("Data Scaled, Dimensionality Reduction using Principal Component Analysis")

    pca = PCA(n_components = combined[x][0],svd_solver = 'randomized',random_state=0)
    pca_df = pd.DataFrame(pca.fit_transform(scaled_data.transpose()))

    logger.info("Reduced Dimension with PCA and random state 0 and n_components : {}".format(combined[x][0]))

    reducer = umap.UMAP(n_neighbors = combined[x][1] , min_dist = combined[x][2] , metric = 'euclidean')

    logger.info("Initilaizing UMAP with neighbors :"+b+", distance :"+c+" and Euclidean metric")

    embedding = reducer.fit_transform(pca_df)


    final_df = pd.DataFrame(embedding,columns=['Factor1','Factor2'])
    final_df.to_csv('Reduced_Coordinates_'+a+'_'+b+'_'+c+'.csv',index=False)

    population_df = pd.read_csv(population_file,sep='\t')
    panel_df = pd.read_csv(panel_file,sep = '\t')
    panel_df['Population_Name'] = np.nan

    logger.info("Pre-Processing and Connecting Data from Poulations and Panel files for Plots")
    for i,row1 in panel_df.iterrows():
        for j,row2 in population_df.iterrows():
            if (row1['pop'] == row2['Population Code']):
                panel_df.loc[i,'Population_Name'] = row2['Population Description']
                break

    final_df['pop'] = panel_df['pop']
    final_df['desc'] = panel_df['Population_Name']

    logger.info("Data Processed, Initiating Plot")

    plt.figure(figsize=(50,50))


    for name in final_df['desc'].unique():


        plt.scatter(x=final_df.loc[final_df['desc']==name,'Factor1'], 
                    y=final_df.loc[final_df['desc']==name,'Factor2'], 
                    alpha=0.70,label=final_df.loc[final_df['desc']==name,'pop'].iloc[1])

        #add desc
        plt.annotate(name, 
                     final_df.loc[final_df['desc']==name,['Factor1','Factor2']].mean(),size=15)

    plt.legend(loc='upper right', prop={'size': 20},markerscale = 5)
    plt.savefig(fname = 'Image_'+a+'_'+b+'_'+c+'.pdf', format='pdf', dpi=1000)
    plt.title("Genetic History Visualization")
    logger.info("Plot file saved as : Image_"+a+'_'+b+'_'+c+".pdf")

    stop = time.time()
    logger.info(' TOTAL RUN TIME of Code : {} seconds'.format((stop-start)+(extraction_end-extraction_start)))