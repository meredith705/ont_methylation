## 08/10/2022
## Melissa Meredith
## Usage: python3 ../plot_methyl_fq_correlation.py -i 

## Takes in ONT and bisulfite methyl tags per location and outputs heatmap of methylation frequency

import argparse
import sys
import matplotlib.pyplot as plt
import gzip
import numpy as np
import pandas as pd


###########################################################
# Input File Explanation                                  #
###########################################################
'''
The bedMethyl file is a bed9+2 file containing the number of reads and the percent methylation.

Each column represents the following:

1.Reference chromosome or scaffold
2.Start position in chromosome
3.End position in chromosome
4.Name of item
5."Score" 1000 * (Nmod + Ncanon) / (Nmod + Ncanon + Nfilt + Nsub + Ndel). The quantity reflects the extent to which the calculated modification frequency in Column 11 is confounded by the alternative calls. The denominator here is the total read coverage as given in Column 10.
6.Strandedness, plus (+), minus (-), or unknown (.)
7.Start of where display should be thick (start codon)
8.End of where display should be thick (stop codon)
9.Color value (RGB)
10.Read coverage at reference position including all canonical, modified, undecided (filtered), substitutions from reference, and deletions. Nmod + Ncanon + Nfilt + Nsub + Ndel
11.Percentage of modified bases, as a proportion of canonical and modified (excluding filtered, substitutions, and deletions). 100 * Nmod / (Nmod + Ncanon)
12.Ncanon
13.Nmod
14.Nfilt those bases with a modification probability falling between given thresholds.

1,2,3,10,11

Bismark Coverage (*.bismark.cov.gz)
The MethylSeq v1.0 app provides a Bismark Coverage report in a GZIP compressed format (*.bismark.cov.gz). See Bismark.

Statistic           Definition
1.Chromosome          The chromosome name.
2.Start position      The genomic start position.
3.End position        The genomic end position.
4.Methylation %       The percentage of methylation at that position.
5.Count Methylated    The number of C bases that are methylated.
6.Count Unmethylated  The number of C bases that are unmethylated.

1,2,3,4,5

###########################################################
# Plotting Explanation                                    #
###########################################################

For each chromosome, start, and end position % bisulfite methylation and 
    % ONT mehtylation are incremented in an array[60][70] += 1. 
    Finally the methylation frequency for each technology is plotted as a heatmap.

'''


def increment_ont_counts(ont_input,score_cutoff,read_count_cutoff):
    ''' takes in file and creates, fills and returns positional Methyl% dictinary
        # ont: 1,2,3, 11%, 13 Nmod 100,000'''
    temp_positions = {}
    ont_scores = np.zeros([1001])
    ont_percents = np.zeros([101])
    ont_coverage = np.zeros([10000])
    ont_passFilter_coverage = np.zeros([8000])
    ont_passFilter_percents = np.array([])

    with open(ont_input, 'r') as ont_file:
        i =0
        o = 0
        low_coverage_ont = 0
        low_score = 0
        other_chrom = 0
        # for ont_line,bis_line in zip(ont_file,bisulfite_file):
        for ont_line in ont_file:
            o_line = ont_line.strip().split()
            chromosome = o_line[0]
            start_pos = int(o_line[1])
            end_pos = int(o_line[2])
            ont_score = int(o_line[4])
            total_raw_num_ont_reads = int(o_line[9])
            ont_line_percent = float(o_line[10])
            ont_connical_count = int(o_line[11])
            ont_mod_count = int(o_line[12])

            ont_scores[ont_score]+=1
            ont_percents[round(ont_line_percent)]+=1
            ont_coverage[ont_mod_count+ont_connical_count]+=1

            if ont_mod_count+ont_connical_count >= read_count_cutoff: 
                if ont_score >= score_cutoff:   
                    ont_passFilter_coverage[ont_mod_count+ont_connical_count]+=1
                    # ont_passFilter_percents=np.append(ont_passFilter_percents,np.array([ont_line_percent]))
                    # select only chr## aligned reads?  
                    if len(chromosome)<6:
                        if chromosome not in temp_positions.keys():
                            # positions['chr1']={2301:{ont_percent:90%,'bis_percent':0.0%}}
                            temp_positions[chromosome]={start_pos:{'ont_end_pos':end_pos,
                                                                    'ont_percent':ont_line_percent, 
                                                                    'ont_score':ont_score,
                                                                    'ont_nMod':ont_mod_count, 
                                                                    'ont_nConnical':ont_connical_count, 
                                                                    'ont_nTotal':total_raw_num_ont_reads}}
                            o+=1
                        else:
                            # if start_pos not in list(temp_positions[chromosome].keys())[-100]:
                            temp_positions[chromosome][start_pos]={'ont_end_pos':end_pos,
                                                                    'ont_percent':ont_line_percent, 
                                                                    'ont_nMod':ont_mod_count, 
                                                                    'ont_nConnical':ont_connical_count, 
                                                                    'ont_nTotal':total_raw_num_ont_reads}
                            o+=1
                    else:
                        other_chrom+=1
                else: low_score+=1
            else:
                low_coverage_ont += 1
            if i%10000000 == 0:
                print('ONT line:',i)
            i+=1
    plot_score_bar(ont_percents,"ONT %",'ONTpercents.svg')
    plot_score_bar(np.log(ont_percents), 'ONT %', 'ONTpercentsLog.svg')
    plot_score_bar(ont_scores,"ONT Scores",ont_input[:-4]+'ONTscores.svg')
    plot_score_bar(ont_coverage[1::],"ONT cov",'ONTcoverage.svg')
    plot_score_bar(ont_passFilter_coverage,'pass ONT cov',"ONT_passFilter_coverage.svg")
    # plot_list_bar(np.sort(ont_passFilter_percents), 'ont pass percents','ont_passFilter_percents.svg')
    print('ont lines:',i,'. ONT lines stored', o,'(',(o/i)*100, '%. # of not chr# positions',other_chrom,'(',(other_chrom/i)*100, '%). low score', low_score, '(',(low_score/i)*100,'%). Low coverage positions', low_coverage_ont,'(',(low_coverage_ont/i)*100,'%).')

    return temp_positions, low_coverage_ont

def increment_bis_counts(bfile,pos_dict,read_count_cutoff):
    ''' takes bisulfite file and ont filled dict and adds bis counts
        # bis: 1,2,3,  4%, 5 Nmod '''
    i=0
    b=0
    not_stored=0
    less_than_20_cov = 0
    print('chromosomes:', pos_dict.keys())
    with open(bfile,'r') as bisulfite_file:
        for line in bisulfite_file:
            b_line = line.strip().split()
            chromosome = b_line[0]
            start_pos = int(b_line[1])
            percent = float(b_line[3])
            nMod = int(b_line[4])
            nConnocal = int(b_line[5])

            # if chromosome not in pos_dict.keys():
            #     # positions['chr1']={2301:{ont_percent:90%,'bis_percent':0.0}}
            #     pos_dict[chromosome]={start_pos:{'bis_percent':percent}}#,'ont_percent':0.0}}
            #     b+=1
            # else:
            if nMod+nConnocal >= read_count_cutoff:
                if chromosome in pos_dict.keys():
                    if start_pos in pos_dict[chromosome].keys():
                        pos_dict[chromosome][start_pos]['bis_percent']=percent
                        pos_dict[chromosome][start_pos]['bis_nMod']=nMod
                        pos_dict[chromosome][start_pos]['bis_nConnocal']=nConnocal
                        pos_dict[chromosome][start_pos]['bis_nTotal']=nConnocal+nMod
                        b+=1
                            # pos_dict[chromosome][start_pos]={'bis_percent':percent,'ont_percent':0.0}
                    else:
                        not_stored+=1
                else:
                    not_stored+=1
            else:
                less_than_20_cov+=1
                    
                        # print('added bis:',chromosome,start_pos,pos_dict[chromosome][start_pos])
            # if b%1000000 ==0: 
            #     print('bis count:', b)
            if i%10000000 == 0:
                print('BiS line:',i)
            i+=1

    print('bis lines:',i,'bis lines stored', b,'(',(b/i)*100, '%). Not stored',not_stored, '(',(not_stored/i)*100, '%).low coverage', less_than_20_cov, '(',(less_than_20_cov/i)*100, '%).')

    return pos_dict, less_than_20_cov

def fill_heatmap_array(pos_dict,select_chromosome, plot_filename):
    ''' creates an array [101,101] and increments each square with (ont methyl %, bis methyl %) per site. 
        # positions['chr1']={2301:{ont_percent:90%,'bis_percent':0.0}}
        returns the array '''

    hm_arr = np.zeros([102,102], dtype=int)

    # 1D arrays for barchart
    bis_array = np.zeros([101], dtype=int)
    ont_array = np.zeros([101], dtype=int)
    with open(plot_filename+".bed",'w') as outBed:
        if select_chromosome == 'all':
            for chromosome in list(pos_dict.keys()):
                for position in pos_dict[chromosome].keys():
                    # store ont coverage
                    ont_array[round(pos_dict[chromosome][position]['ont_percent'])]+=pos_dict[chromosome][position]['ont_nTotal']
                    # check bis is in dictionary ( ont is filled first - it will always be there )
                    if 'bis_percent' in pos_dict[chromosome][position].keys():
                        row=round(pos_dict[chromosome][position]['ont_percent'])
                        col=round(pos_dict[chromosome][position]['bis_percent'])
                        # store bis field
                        bis_array[col]+=pos_dict[chromosome][position]['bis_nTotal']

                        if row > 100 or col > 100:
                            print('row',row,'col',col)
                        # arr[roman][catholics]
                        hm_arr[row][col]+=1
                        outstring='\t'.join([chromosome,str(position),str(pos_dict[chromosome][position]['ont_percent']),str(pos_dict[chromosome][position]['bis_percent'])])
                        outBed.write(outstring+'\n')


        else:
            for position in pos_dict[select_chromosome].keys():
                # store ont coverage
                ont_array[round(pos_dict[select_chromosome][position]['ont_percent'])]+=pos_dict[select_chromosome][position]['ont_nTotal']
                # check bis is in dictionary ( ont is filled first - it will always be there )
                if 'bis_percent' in pos_dict[select_chromosome][position].keys():
                    row=round(pos_dict[select_chromosome][position]['ont_percent'])
                    col=round(pos_dict[select_chromosome][position]['bis_percent'])
                    # store bis field
                    bis_array[col]+=pos_dict[select_chromosome][position]['bis_nTotal']

                    if row > 100 or col > 100:
                        print('row',row,'col',col)
                    # arr[roman][catholics]
                    hm_arr[row][col]+=1

    # add pseudo counts
    print('add pseudo counts')
    for i in np.arange(hm_arr.shape[0]):
        for j in np.arange(hm_arr.shape[1]):
            hm_arr[i][j]+=1

    # for a in np.arange(arr.shape[0]):
    #     print(arr[a])
    print('\tlog transform')
    log_data = np.log(hm_arr)
    print('\tcoeff calc')
    print(np.corrcoef(ont_array,bis_array))

    write_out_array('raw.'+plot_filename, hm_arr)

    half_log_data = half_the_array(log_data)

    return half_log_data, log_data, hm_arr, bis_array, ont_array

def half_the_array(big_array):
    ''' take the [101,101] array and make it [51,51] by averaging the [i][j],[i][j+1],[i+1][j],[i+1][j+1] positions
        hal_array[i][j] = avg(arr[i][j], arr[i][j+1], arr[i+1][j], arr[i+1][j+1]) '''
    
    size = big_array.shape[0]
    print('size',size,size//2, size//2*2*size//2*2)
    blocks = big_array.reshape(size//2,2,size//2, 2)
    block_mean = np.mean(blocks, axis=(1,-1))

    return block_mean

def write_out_array(plot_filename, hm_array):
    with open(plot_filename+'.hm_array.csv','w') as outfile:
        for i in hm_array:
            counter = 0
            for j in i:
                counter+=1
                if counter < hm_array.shape[1]:
                    outfile.write(str(j)+",")
                else:
                    outfile.write(str(j))
            outfile.write('\n')

def fill_fq_array(this_dict,select_field,select_chromosome):
    ''' creates an array [100] and increments each square with (ont methyl % or bis methyl %) per site. 
        # positions['chr1']={2301:{ont_percent:90%,'bis_percent':0.0}}
        returns the array '''

    bis_array = np.zeros([101], dtype=int)

    ont_array = np.zeros([101], dtype=int)

    if select_chromosome == 'all':
        # loop through all the chromosomes
        for chromosome in list(this_dict.keys()):
            for position in this_dict[chromosome].keys():
                # only include positions shared between the two technologies 
                # if len(this_dict[chromosome][position])==2:
                bis=round(this_dict[chromosome][position]['bis_'+select_field])
                bis_array[bis]+=1

                ont=round(this_dict[chromosome][position]['ont_'+select_field])
                ont_array[ont]+=1
    else: 
        for position in this_dict[select_chromosome].keys():
            bis=round(this_dict[select_chromosome][position][select_field])
            bis_array[bis]+=1

            ont=round(this_dict[chromosome][position]['ont_'+select_field])
            ont_array[ont]+=1
    # print(arr)
    return bis_array, ont_array

def plot_heatmap(meth_arr,fileName):
    ''' plots a heatmap from the positional methyl % array '''

    fig, ax = plt.subplots()
    im = ax.imshow(meth_arr,origin='lower')

    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('number of sites ('+fileName[0:3]+')', rotation=-90, va="bottom")

    ax.set_xticks(np.arange(0,meth_arr.shape[1]+1,20))
    ax.set_yticks(np.arange(0,meth_arr.shape[0]+1,20))

    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)

    ax.set_title("Metylation Comparison")
    ax.set_ylabel("ONT Methylation %")
    ax.set_xlabel("Bisulfite Methylation %")

    # fig.tight_layout()
    
    plt.savefig(fileName+".svg")

    # plt.show()

def plot_bar(ont_array,ont_label,bis_array,bis_label,outfileName):
    ''' plots bar chart of both technologies' arrays '''

    fig = plt.figure(figsize = (10,5))

    plt.bar(np.arange(0,101), ont_array, alpha=0.7, label=ont_label)
    plt.bar(np.arange(0,101), bis_array, alpha=0.7, label=bis_label)
    plt.xticks(np.arange(0,bis_array.shape[0]+1,5))
    plt.title(outfileName[:-4])
    plt.legend( loc ="upper right")

    # plt.show()#
    plt.savefig(outfileName)

def plot_score_bar(ont_array,ont_label,outfileName):
    ''' plots bar chart of both technologies' arrays '''

    fig = plt.figure(figsize = (10,5))

    plt.bar(np.arange(0,ont_array.shape[0]), ont_array, alpha=0.7, label=ont_label)
    # plt.xticks(np.arange(0,ont_array.shape[0]+1))
    plt.title(outfileName[:-4])
    # plt.legend( loc ="upper left")

    # plt.show()#
    plt.savefig(outfileName)

def plot_list_bar(ont_array,ont_label,outfileName):
    ''' plots bar chart of both technologies' arrays '''

    fig = plt.figure(figsize = (10,5))

    plt.bar(np.arange(0,len(ont_array)), ont_array, alpha=0.7, label=ont_label)
    # plt.xticks(np.arange(0,ont_array.shape[0]+1))
    plt.title(outfileName[:-4]+str(np.median(ont_array)))
    # plt.legend( loc ="upper left")

    # plt.show()#
    plt.savefig(outfileName)

def arInput(arr_file):
    ''' Reads in data from an array.csv file if it is provided as input '''
    file = open(arr_file, 'rb')
    data = np.loadtxt(file,delimiter = ",", usecols=range(101), dtype='int')

    print('input data shape',data.shape)

    # add pseudo counts
    for i in np.arange(data.shape[0]):
        for j in np.arange(data.shape[1]):
            data[i][j]+=1
    
    # log transform the array and make a correlation matrix    
    log_data = np.log(data)
    cceff = np.corrcoef(data)

    return log_data,data,cceff



def main(ont_input,bisulfite_input,arr_input,plot_filename,ont_score_cutoff,read_count_cutoff,select_chromosome):

    print("ONT bed",ont_input,"\nbisulfite cov.gz", bisulfite_input)
    print('ont_score_cutoff:', ont_score_cutoff, 'read_count_cutoff', read_count_cutoff)

    print('reading ont file')
    # fill dictionary with ont positions and percentages
    cpg_positional_dict, low_coverage_ont = increment_ont_counts(ont_input=ont_input,score_cutoff=ont_score_cutoff,read_count_cutoff=read_count_cutoff)
    # ont_dict = increment_ont_counts(ont_input,'ONT')
    

    print('reading bis file')
    # add bis positions and percentages
    cpg_positional_dict, low_coverage_bis = increment_bis_counts(bfile=bisulfite_input,pos_dict=cpg_positional_dict,read_count_cutoff=read_count_cutoff)
    # bis_dict={}
    # bis_dict = increment_bis_counts(bisulfite_gz_input,bis_dict)

    print('filling array')
    final_outfile_name = str(ont_score_cutoff)+"."+str(read_count_cutoff)+"."+plot_filename
    half_log_hm_array, log_hm_array, hm_array, bis_array, ont_array = fill_heatmap_array(cpg_positional_dict,select_chromosome, final_outfile_name)


    # ont_arr, bis_arr = fill_fq_array(positional_dict,'percent','all')
    print('ont array', ont_array, '\nbis array', bis_array)
    plot_bar(ont_array,"Nanopore",bis_array,"Bisulfite",final_outfile_name+'.'+select_chromosome+'.barchart.svg')

    # log_hm_array, hm_array, cceff = arInput(arr_input)
    print('plotting heatmap')
    plot_heatmap(log_hm_array,'log.'+plot_filename)
    plot_heatmap(half_log_hm_array,'log.bin2.'+plot_filename)
    plot_heatmap(hm_array,'raw.'+plot_filename)

    print('done')
    print('low coverage ONT:', low_coverage_ont, '# low coverage BIS:', low_coverage_bis)

    print(hm_array)
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="Input ONT methlyation bed"
    )

    parser.add_argument(
        "-b",
        required=True,
        type=str,
        help="Input bisulfite methlyation cov.gz"
    )

    parser.add_argument(
        "-n",
        required=True,
        type=str,
        help="plot filename"
    )

    parser.add_argument(
        "-c",
        required=False,
        type=str,
        default='all',
        help="chromosome to make plots from (eg. chr1 or all)"
    )

    parser.add_argument(
        "-a",
        required=False,
        type=str,
        help="Input arr"
    )

    parser.add_argument(
        "-s",
        required=False,
        type=int,
        default=1000,
        help="ONT bedMethyl min score cutoff (0-1000); default=1000"
    )

    parser.add_argument(
        "-r",
        required=False,
        type=int,
        default=20,
        help="Read coverage cutoff for any CpG site to be included in plot; default=20"
    )


    args = parser.parse_args()

    main(ont_input=args.i, bisulfite_input=args.b, arr_input=args.a, plot_filename=args.n, ont_score_cutoff=args.s, read_count_cutoff=args.r, select_chromosome=args.c)

