import argparse
from numpy import mean
from scipy.stats import pearsonr

def get_alt_allele_freq(ids,snpline,pop1,pop2):
    snplist = snpline.strip().split()
    chrom = snplist[0]
    position = snplist[1]
    snpdict = dict(zip(ids,snplist[9:]))
    numalleles = 1 + len(snplist[4].split(','))
    pop1_total = 0
    pop1_altcount = 0
    pop2_total = 0
    pop2_altcount = 0
    for pop in pop1.split(','):
        genotype = snpdict[pop][0]
        if len(genotype) > 3:
            raise ValueError('{} is not a valid format'.format(genotype))
        elif genotype !='./.' and numalleles <= 2:
            pop1_total += (genotype.count('0') + genotype.count('1'))
            pop1_altcount += genotype.count('1')
    for pop in pop2.split(','):
        genotype = snpdict[pop][0]
        if len(genotype) > 3:
            raise ValueError('{} is not a valid format'.format(genotype))
        elif genotype !='./.' and numalleles <= 2:
            pop2_total += (genotype.count('0') + genotype.count('1'))
            pop2_altcount += genotype.count('1')
    
    if pop1_total != 0:
        pop1_alt_frq = pop1_altcount/pop1_total
    else:
        pop1_alt_frq = 'NA'

    if pop2_total != 0:
        pop2_alt_frq = pop2_altcount/pop2_total
    else:
        pop2_alt_frq = 'NA'

    return chrom,position,pop1_alt_frq,pop2_alt_frq

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='generates allele frequencies for two populations from a vcf file')
    parser.add_argument('-v','--vcf',dest='vcf',type=str,help='vcf file used to calculate allele frequencies')
    parser.add_argument('-pop1','--sample-pop-1',dest='p1',type=str,help='first population for which to get allele freqs')
    parser.add_argument('-pop2','--sample-pop-2',dest='p2',type=str,help='second population for which to get allele freqs')
    parser.add_argument('-o','--frequency-outfile',dest='outfile',type=str,help='output file name to write freqs to')
    opts = parser.parse_args()

    fout = open(opts.outfile,'w')
    fout.write('pop1 is {}\n'.format(opts.p1))
    fout.write('pop2 is {}\n'.format(opts.p2))
    fout.write('chromosome\tposition\tpop1_altfrq\tpop2_altfrq\n')
    
    pop1_frq_list = []
    pop2_frq_list = []
    pop1_frq_avail = []
    pop2_frq_avail = []
    vcfin = open(opts.vcf,'r') 
    for line in vcfin:
        if line[0] == '#' and 'CHROM' in line:
            sample_ids = line.strip().split()[9:]
        elif line[0] != '#':
            chrom,position,pop1frq,pop2frq = get_alt_allele_freq(sample_ids,line,opts.p1,opts.p2)
            fout.write('{}\t{}\t{}\t{}\n'.format(chrom,position,pop1frq,pop2frq))
            if pop1frq != 'NA':
                pop1_frq_list.append(pop1frq)
            if pop2frq != 'NA':
                pop2_frq_list.append(pop2frq)
            if pop1frq != 'NA' and pop2frq != 'NA':
                pop1_frq_avail.append(pop1frq)
                pop2_frq_avail.append(pop2frq)
    fout.close()

    print('number of observable pop1 alt freqs = {}'.format(len(pop1_frq_list)))
    print('mean pop1 alt frq = {}\n'.format(mean(pop1_frq_list)))
    print('number of observable pop2 alt freqs = {}'.format(len(pop2_frq_list)))
    print('mean pop2 alt frq = {}\n'.format(mean(pop2_frq_list)))
    print('allele frequency correlation between callable sites = {}\n'.format(pearsonr(pop1_frq_avail,pop2_frq_avail)))
