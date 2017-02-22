from sets import Set
from collections import defaultdict
from subprocess import Popen,PIPE
import argparse
from os import system,path


def extract_uniquely_mapped(bamin):
    query='samtools view -h -F 4 %s |grep -w "X0:i:1" > unique_%s.sam' % (bamin,bamin.split('.')[0]) 
    p1a=Popen(query,shell=True,stderr=PIPE,stdout=PIPE)
    stdout,stderr=p1a.communicate()
    grabhead='samtools view -H %s > header' % bamin
    p1b=Popen(grabhead,shell=True,stderr=PIPE,stdout=PIPE)
    stdout2,stderr2=p1b.communicate()
    system('cat header unique_%s.sam > wheader_unique_%s.sam' % (bamin.split('.')[0],bamin.split('.')[0]))
    system('rm header')
    system('samtools view -Sbh wheader_unique_%s.sam > unique_%s.bam' % (bamin.split('.')[0],bamin.split('.')[0]))
    system('rm wheader_unique_%s.sam' % bamin.split('.')[0]) 
    system('rm unique_%s.sam' % bamin.split('.')[0])
    if p1a.returncode==0 and p1b.returncode==0:
        return True
    else:
        return stderr

def bam2bed(bam):
    cmd='bamToBed -i %s > %s.bed' % (bam,bam.split('.')[0])
    p2=Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    stdout,stderr=p2.communicate()
    if p2.returncode==0:
        return True
    else:
        return stderr    
    
def extend_gene_boundaries(geneboundaries,numbases,outbed):
    fopen=open(geneboundaries,'r')
    fout=open(outbed,'w')
    for line in fopen:
        chrom,start,stop,strand,id=line.strip().split()
        if strand=="+":
            stop=int(stop)+numbases
        else:
            start=int(start)-numbases
        fout.write('%s\t%s\t%s\t%s\t%s\n' % (chrom,start,stop,strand,id))
    fout.close()


def intersect_reads_genes(readbed,genebed,outfile):
    cmd='intersectBed -wo -a %s -b %s > %s' % (readbed,genebed,outfile)
    p3=Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    stdout,stderr=p3.communicate()
    if p3.returncode==0:
        return True
    else:
        return False

def enforce_strandedness(intersectbed):
    """
    takes a bed file of intersections between
    read bed and gene boundary bed files
    and filters out entries where the read mapping
    and gene are on different strands.
    """ 
    strand_dict={'same':0,'diff':0}
    fopen=open(intersectbed,'r')
    fout=open('samestrand_'+intersectbed,'w')
    for line in fopen:
        linelist=line.strip().split()
        if linelist[5]==linelist[9]:
            fout.write(line)
            strand_dict['same']+=1
        else:
            strand_dict['diff']+=1
    
    fout.close()
    return strand_dict
   
    
def parse_intersect_to_dict(intersectbed):
    """
    convert the bed file resulting from intersectbed
    between uniquely mapping reads bed and gene
    boundaries into a dictionary with read id/umi
    combinations as keys, and values are dictionaries
    for each line in the bed file associated with that
    read id/umi combo
    """
    intersect_dict=defaultdict(list)
    bedin=open(intersectbed,'r')
    for line in bedin:
        linedict={}
        bedlist=line.strip().split()
        read_id=bedlist[3]
        read_instrument_data=read_id.split('#')[0]
        read_index=read_id.split('#')[1].split('/s')[0]
        read_umi=read_id.split('#')[1].split('/s')[1]

        linedict['index']=read_index
        linedict['chrom']=bedlist[0]
        linedict['read_zero_start']=int(bedlist[1])
        linedict['read_end']=int(bedlist[2])
        linedict['mapq']=int(bedlist[4])
        linedict['read_strand']=bedlist[5]
        linedict['gene_zero_start']=int(bedlist[7])
        linedict['gene_end']=int(bedlist[8])
        linedict['gene_strand']=bedlist[9]
        linedict['gene_id']=bedlist[10]
        linedict['overlap']=int(bedlist[11])
       
        intersect_dict[(read_instrument_data,read_umi)].append(linedict)

    return intersect_dict

def remove_multigene_reads(intersection_dictionary):
    """
    removes reads for which intersect bed operation with
    gene boundaries yields > 1 gene, to remove erroneous
    gene assignments in ambiguous cases
    """

    filtered_dict=defaultdict(list)
    for key in intersection_dictionary.keys():
        if len(intersection_dictionary[key])>1:
            genes=Set()
            for dictionary in intersection_dictionary[key]:
                genes.add(dictionary['gene_id'])
            if len(genes)==1:
                filtered_dict[key]=intersection_dictionary[key]    
        else:
            filtered_dict[key]=intersection_dictionary[key]

    return filtered_dict 

def build_gene_to_reads_dict(reads_dictionary):
    gene_dict=defaultdict(list)              
    for key in reads_dictionary.keys():
        for dictionary in reads_dictionary[key]:
            dictionary['read_id']=key
            gene_dict[dictionary['gene_id']].append(dictionary)

    return gene_dict
        
def extract_readcount_from_gene(gene_key,gene_dict):
    gene_set=Set()
    for read_dict in gene_dict[gene_key]:
        gene_set.add(read_dict['read_id'])
    return len(gene_set)
        


if  __name__=="__main__": 
    parser = argparse.ArgumentParser(description="single cell 3' end RNA=seq counts by feature and UMI collapse")
    parser.add_argument('-bam','--bamfile_in',dest='bam',type=str,help='in bamfile')
    parser.add_argument('-gbed','--gene_boundaries',dest='gbed',type=str,help='bed file of gene boundaries')
    parser.add_argument('-3ext','--three_prime_extend',dest='threeprime',type=int,default=0,help='# of bases to extend bedtools query in 3 prime direction')
    parser.add_argument('-intout','--read_gene_intersect_out',dest='interout',type=str,help='name of file that is intersection of unique reads and gene beds')
    parser.add_argument('-count','--gene_count_table',dest='gcount',type=str,help='name of gene count table outfile')
    parser.add_argument('-strand','--enforce_strand',dest='strand',type=bool,default=True,help='require mappings and genes be on same strand to intersect')
    opts = parser.parse_args()    


    # open stats log #
    log_handle=open('%s.log' % opts.bam.split('.')[0],'w3')

    # extend 3' utr boundaries if argument invoked
    if opts.threeprime>0:
        geneboundaries='extendby%sbp_%s' % (opts.threeprime,path.basename(opts.gbed))
        extend_gene_boundaries(opts.gbed,opts.threeprime,geneboundaries)
    else:
        geneboundaries=opts.gbed    
    

    # create uniquely mapped reads bam #
    unique_extract=extract_uniquely_mapped(opts.bam)
    if unique_extract==True:
        uniquebam='unique_%s' % opts.bam
        bed_from_unique=bam2bed(uniquebam)
        if bed_from_unique==True:
            uniquebed='unique_'+opts.bam.split('.')[0]+'.bed'
            intersect_reads_genes(uniquebed,geneboundaries,opts.interout)
            if opts.strand==True:
                strand_enforce=enforce_strandedness(opts.interout)
                intersect_dict=parse_intersect_to_dict('samestrand_%s' % opts.interout)
                strandlog=open('strandconcordance_%s.metrics' % opts.bam.split('.')[0],'w')
                strandlog.write('same:\t%s\ndiff:\t%s\n' % (strand_enforce['same'],strand_enforce['diff']))
                strandlog.close()
            else:
                intersect_dict=parse_intersect_to_dict(opts.interout)
            
            filter_dict=remove_multigene_reads(intersect_dict) # removes reads that overlap with more than one gene
          
            log_handle.write('Total reads intersecting genes: %s\nTotal reads overlapping only one gene: %s\n' % (len(intersect_dict),len(filter_dict))) 
  
            gene_to_read_dict=build_gene_to_reads_dict(filter_dict)
            counts_out=open(opts.gcount,'w')
            counts_out.write('geneid\t%s\n' % opts.bam.split('.')[0])
            genekeys=gene_to_read_dict.keys()
            genekeys.sort()
            for gene in genekeys:
                gene_count=extract_readcount_from_gene(gene,gene_to_read_dict)
                counts_out.write('%s\t%s\n' % (gene,gene_count))
            counts_out.close()
    
        else:
            print bed_from_unique 



    else:
        print 'ERROR:%\n%s' % unique_extract
        
    



   
        
