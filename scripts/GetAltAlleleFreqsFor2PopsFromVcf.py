import argparse
import re

def pairwise_difference(alt_count, total):
    """Fraction of unordered pairs of called allele copies that differ at a site."""
    if total < 2:
        return 'NA'
    return 2 * alt_count * (total - alt_count) / (total * (total - 1))

def population_alt_counts(snpdict, population, genotype_index):
    total = 0
    alt_count = 0
    for sample_id in population.split(','):
        fields = snpdict[sample_id].split(':')
        if genotype_index >= len(fields):
            raise ValueError('missing GT field for {}'.format(sample_id))
        genotype = fields[genotype_index]
        alleles = re.split(r'[/|]', genotype)
        if any(allele not in ('0', '1', '.') for allele in alleles):
            raise ValueError('invalid biallelic genotype {!r} for {}'.format(genotype, sample_id))
        total += sum(allele != '.' for allele in alleles)
        alt_count += sum(allele == '1' for allele in alleles)
    return alt_count, total

def get_site_frequencies(ids, snpline, pop1, pop2):
    snplist = snpline.strip().split()
    if len(snplist) != 9 + len(ids):
        raise ValueError('VCF record has {} sample columns; expected {}'.format(
            len(snplist) - 9, len(ids)))
    chrom = snplist[0]
    position = snplist[1]
    if snplist[4] == '.' or ',' in snplist[4]:
        return chrom, position, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'
    format_fields = snplist[8].split(':')
    if 'GT' not in format_fields:
        raise ValueError('VCF record at {}:{} has no GT field'.format(chrom, position))
    genotype_index = format_fields.index('GT')
    snpdict = dict(zip(ids, snplist[9:]))
    pop1_alt_count, pop1_total = population_alt_counts(snpdict, pop1, genotype_index)
    pop2_alt_count, pop2_total = population_alt_counts(snpdict, pop2, genotype_index)
    pop1_phi = pairwise_difference(pop1_alt_count, pop1_total)
    pop2_phi = pairwise_difference(pop2_alt_count, pop2_total)
    phi_difference = (pop1_phi - pop2_phi
                      if pop1_phi != 'NA' and pop2_phi != 'NA' else 'NA')

    pooled_alt_count = pop1_alt_count + pop2_alt_count
    pooled_total = pop1_total + pop2_total
    if not pooled_total:
        return chrom, position, 'NA', 'NA', 'NA', pop1_phi, pop2_phi, phi_difference

    # At a pooled-frequency tie, choose ALT consistently. Differences are
    # invariant to this choice because both populations use the same allele.
    minor_is_alt = 2 * pooled_alt_count <= pooled_total
    pop1_minor_freq = ((pop1_alt_count if minor_is_alt else pop1_total - pop1_alt_count)
                       / pop1_total if pop1_total else 'NA')
    pop2_minor_freq = ((pop2_alt_count if minor_is_alt else pop2_total - pop2_alt_count)
                       / pop2_total if pop2_total else 'NA')
    difference = (abs(pop1_minor_freq - pop2_minor_freq)
                  if pop1_total and pop2_total else 'NA')
    return (chrom, position, pop1_minor_freq, pop2_minor_freq, difference,
            pop1_phi, pop2_phi, phi_difference)

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description='generates minor allele frequencies and pairwise diversity for two populations from a VCF')
    parser.add_argument('-v','--vcf',dest='vcf',type=str,required=True,help='vcf file used to calculate allele frequencies')
    parser.add_argument('-pop1','--sample-pop-1',dest='p1',type=str,required=True,help='first population for which to get allele freqs')
    parser.add_argument('-pop2','--sample-pop-2',dest='p2',type=str,required=True,help='second population for which to get allele freqs')
    parser.add_argument('-o','--frequency-outfile',dest='outfile',type=str,required=True,help='output file name to write freqs to')
    opts = parser.parse_args()

    with open(opts.vcf) as vcfin:
        sample_ids = None
        for line in vcfin:
            if line.startswith('#CHROM'):
                sample_ids = line.strip().split()[9:]
                break
        if sample_ids is None:
            raise ValueError('VCF is missing the #CHROM header')
        if len(sample_ids) != len(set(sample_ids)):
            raise ValueError('VCF has duplicate sample IDs')
        for population in (opts.p1, opts.p2):
            ids = population.split(',')
            missing = [sample_id for sample_id in ids if sample_id not in sample_ids]
            if missing:
                raise ValueError('sample IDs absent from VCF: {}'.format(', '.join(missing)))

    with open(opts.vcf) as vcfin, open(opts.outfile, 'w') as fout:
        fout.write('chromosome\tposition\tpop1_pooled_minor_frq\t'
                   'pop2_pooled_minor_frq\tabs_minor_frq_diff\t'
                   'pop1_phi\tpop2_phi\tphi_diff_pop1_minus_pop2\n')
        for line in vcfin:
            if not line.strip() or line.startswith('#'):
                continue
            result = get_site_frequencies(sample_ids, line, opts.p1, opts.p2)
            fout.write('\t'.join(map(str, result)) + '\n')
