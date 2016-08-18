import pysam
import pickle
import sys
import math
import itertools

from collections import defaultdict

from Bio import SeqIO

from codon_usage_memex import CodonUsageMemex
from codon_usage_memex import build_codon_usage_dict

# Amount to trim from read ends that is not included in the profile. Li 2012 set this @ 12
PROFILE_TRIM = 12

#List of bam file names
BAMFILE_NAMES = [
    "/scratch/apps/genome-designer/data/users/8fc1f831/projects/965a907d/genomes/43f65f2e/bowtie2_align.bam",
    "/scratch/apps/genome-designer/data/users/8fc1f831/projects/965a907d/genomes/b42eff81/bowtie2_align.bam",
    "/scratch/apps/genome-designer/data/users/8fc1f831/projects/965a907d/genomes/e18925d1/bowtie2_align.bam",
    "/scratch/apps/genome-designer/data/users/8fc1f831/projects/965a907d/genomes/e794cb3f/bowtie2_align.bam"]
     
PROFILE_NAME = "/scratch/dbg/genome-refactor/mds42_weissman_profile.pickled"

CTE_SCORE_FILE = '/scratch/dbg/genome-refactor/data/codon_cte_scores.txt'

MRNA_EXPRESSION_FILE = '/scratch/dbg/genome-refactor/data/GSM849371_genes.csv'

# The read position offset for the ribosomal A-site, used to identify 
# the per-codon pause average.
# Weissman found that it was between -8 and -11 (thats where
# the strong pause at stop codons was) and I found it similarly in three
# example sequences.

A_SITE_OFFSET = -8
A_SITE_WINDOW = 4

MASK_5PRIME = 10
MASK_3PRIME = 10

def get_per_gene_mrna_data(genome):
    '''
    from the MRNA_EXPRESION_FILE, extract a dictionary of per-gene mRNA
    expression levels for use in other functions.
    
    It returns a modified genome object with this RPKM and DPKM data added.
    RPKM: analog transcriptome read data: 
        (reads per kilobase per million mapped reads)
    DPKM: barcode based digital read data
        (digital unique read barcodes per kilobase per million mapped reads)
    '''
    
    rpkm = {}
    dpkm = {}
    missing = 0
    present = 0
    
    for line in open(MRNA_EXPRESSION_FILE):
        fields = line.split(',')
        rpkm[fields[0]] = fields[6]
        dpkm[fields[0]] = fields[7]
    
    for feature in genome.features:
                
        #only take CDS features:
        if feature.type != 'CDS' or 'gene' not in feature.qualifiers:
            continue
        
        gene_id = feature.qualifiers['gene'][0]
        
        if gene_id not in rpkm.keys():
            missing += 1
            continue
        
        present += 1
        feature.qualifiers['rpkm'] = float(rpkm[gene_id])
        feature.qualifiers['dpkm'] = float(dpkm[gene_id])
    
    print '%d genes missing, %d present in rna-seq data.' % (missing, present)

    return genome


def add_alignment_to_profile(profile, aligned_read, multi=1):
    '''
    Takes a pysam aligned read, trims it, and then adds it to the reads profile
    object. If this is one of multiple alignments for the same read, then only add 
    1 / multi to each base of the alignment. 
    '''
    
    # Get the strand
    strand = 'rev' if aligned_read.is_reverse else 'fwd'

    # Get the read start and end
    (r_start,r_end) = (aligned_read.aend - aligned_read.alen, aligned_read.aend)

    # Get profile start and end
    (p_start,p_end) = (r_start + PROFILE_TRIM, r_end - PROFILE_TRIM)

    # Get profile length, and divide by length to get score per nucleotide
    p_len = p_end - p_start
    p_score = float(1) / ( float(p_len) * float(multi) )

    for i in range(p_start, p_end):
        profile[strand][i] += p_score

def add_alignments_by_read(profile, bamfile):
    '''
    Go through each read adding bases to profiles. We need to consider multiple alignments 
    for the same read together, so we'll look ahead to see if the next alignment is the 
    same read.
    '''

    alignments_for_read = []
    current_read_name = None

    try:
        while True:
            
            # Go through all alignments, sorted by read name, skipping unmapped reads
            aligned_read = bamfile.next()
            if aligned_read.is_unmapped: continue
        
            # If this is a new read, add all alignments for previous read to profile, 
            # and clear the alignments list
            if current_read_name != None and current_read_name != aligned_read:
                for alignment in alignments_for_read:
                    add_alignment_to_profile(profile, alignment, len(alignments_for_read))
                alignments_for_read = []
        
            # Update the read name and append the aligned read to the alignments list
            current_read_name = aligned_read.qname
            alignments_for_read.append(aligned_read)
        
    except StopIteration:
            # Add alignments for final read
            for alignment in alignments_for_read:
                add_alignment_to_profile(profile, alignment, len(alignments_for_read))
                
def extract_region_profile(profile, strand, start, end):
    '''
    Extract a profile for a specific region. Strand must be 'fwd' or 'rev'.
    '''
    
    region_profile = {}
    
    #we need to adjust the start and end to account for the ribosomal A-site
    #positioning. Weissman found that it was between -8 and -11 (thats where
    #the strong pause at stop codons was) and I found it similarly in three
    #example sequences. 
    
    #if we are on the reverse strand then we need to reverse the order
    if strand == 'rev':
        (start, end) = (end-1, start-1)
        direction = -1
    elif strand == 'fwd':
        direction = 1
    else:
        raise UserError('strand must be either "fwd" or "rev"')
        
    start = start+(A_SITE_OFFSET*direction)
    end = end+(A_SITE_OFFSET*direction)
    
    for region_i, genome_i in enumerate(range(start, end, direction)):
        region_profile[region_i] = profile[strand][genome_i]
    
    return region_profile
    
def create_empty_profile():
    '''
    Creates an empty per-base ribosome profile holding 'fwd' and 'rev' strands.
    '''
    return {'fwd': defaultdict(float), 'rev': defaultdict(float)}

def create_genomewide_profile(bamfile_names):
    '''
    Load in all the bam/sam read files and create a ribosome profile.
    '''
    
    profile = create_empty_profile()
    
    for bamfile_name in bamfile_names:
        bamfile = pysam.Samfile(bamfile_name, "rb")
        add_alignments_by_read(profile, bamfile)
    
    return profile

def save_profile(profile, file_name):
    pickle.dump(profile, open(file_name, 'wb'))

def load_profile(file_name):
    return pickle.load(open(file_name, 'rb'))

def merge_profiles(profile_a, profile_b):
    for strand in ('fwd','rev'):
        for base, score in profile_a[strand].items():
            profile_b[strand][base] += score
    
    return profile_b

def merge_all_profiles(profile_filenames):
    
    profiles = []
    
    for filename in profile_filenames:
        profiles.append(load_profile(filename))
    
    merged_profile = profiles.pop()
    
    for profile in profiles:
        merged_profile = merge_profiles(profile, merged_profile)
    
    return merged_profile
    
def add_translation_info(genome, profile):
    '''
    This function generates a large amount of per-gene and per-codon data 
    based on ribosome profiling, per-codon occupancy, NTE/CTE data, etc. 
    
    It creates three objects, (genome, codons, genes, bad_cds) containing
    this data. The new genome SeqRecord object has additional data appended
    to each SeqFeature regarding occupancy and mRNA transcription level. 
    
    Additionally, return a codon table with CTE and NTE usage for E. coli.
    '''
    
    codons = defaultdict(lambda: {
        'count' : 0,
        'count_total' : 0,
        'count_nterm' : 0,
        'rnaseq_count' : 0,
        'occup' : 0,
        'logmean_occup' : 0,
        'mean_occup' : 0,
        'cte' : 0 ,
        'rpkm' : 0,
        'dpkm' : 0,
        'rpkm_nterm' : 0,
        'dpkm_nterm' : 0,
        'occup_nterm' : 0,
        'logmean_occup_nterm' : 0,
        'mean_occup_nterm' : 0})
        
    genes = defaultdict(lambda: {
        'logmean_occup' : 0,
        'mean_occup' : 0,
        'logmean_occup_nterm' : 0,
        'mean_occup_nterm' : 0,
        'rpkm' : 0,
        'dpkm' : 0})
        
    hexamers = defaultdict(lambda: {
        'count' : 0,
        'count_total' : 0,
        'count_nterm' : 0,
        'rnaseq_count' : 0,
        'occup' : 0,
        'logmean_occup' : 0,
        'mean_occup' : 0,
        'occup_nterm' : 0,
        'logmean_occup_nterm' : 0,
        'mean_occup_nterm' : 0,
        'rpkm' : 0,
        'dpkm' : 0,
        'rpkm_nterm' : 0,
        'dpkm_nterm' : 0})
        
    
    cds = defaultdict(dict)
    
    bad_cds = {}

    for feature in genome.features:
                
        #only take CDS features:
        if feature.type != 'CDS': continue
        
        #skip if no gene id
        try: gene_id = feature.qualifiers['gene'][0]
        except(KeyError): continue
        
        #Li/Weissman don't use these genes due to frameshift/homology
        if gene_id in ('tufA','tufB','prfX','dnaX'):
            continue
            
        #get the profile for this region
        strand = 'fwd' if feature.strand == 1 else 'rev'
        
        #skip if this gene is too short - less than 30 aa
        if len(feature) < 90: 
            bad_cds[feature.qualifiers['gene'][0]] = 'Too short (%d)' % len(feature)
            continue
            
        #get the DNA sequence for this feature
        feature_location_slice = slice(
            feature.location.start.position,
            feature.location.end.position)
        feature_seq = genome.seq[feature_location_slice]
        if strand == 'rev': feature_seq = feature_seq.reverse_complement()
        
        #skip if this CDS is not a multiple of 3
        if not (len(feature_seq) % 3 == 0):
            bad_cds[feature.qualifiers['gene'][0]] = 'Not codon divisible'
            continue

        #skip if this CDS does not start with 'NTG'
        if not (feature_seq[1:3].tostring() == 'TG'): 
            bad_cds[feature.qualifiers['gene'][0]] = 'First codon is ' + \
                    feature_seq[0:3].tostring()
            continue

        #skip if this CDS does not end with a stop codon (TGA, TAT, TAG)
        final_codon = feature_seq[len(feature)-3:len(feature)].tostring()
        if not (final_codon in ['TGA','TAA','TAG']):
            bad_cds[feature.qualifiers['gene'][0]] = 'Last codon is ' + \
                    final_codon
            continue
                
        #get the ribosomal occupancy
        region_profile = extract_region_profile(
            profile, 
            strand, 
            feature.location.start.position, 
            feature.location.end.position)
                
        #get ribosome occupancy for this gene, but we want to ignore the 
        #first 10 and last 10 codons as Li does
        occupancies = region_profile.values()[30:len(region_profile)-30]
        occupancies_nterm = region_profile.values()[0:30]
        
        #replace 0 occupancy with float min to avoid log of 0 errors
        replace_zero = lambda occ: sys.float_info[3] if occ == 0 else occ
        occupancies = map(replace_zero, occupancies)
        occupancies_nterm = map(replace_zero, occupancies_nterm)
        
        #mean occup for this gene (and for first 10 aa)
        mean_occup = float(sum(occupancies)) / float(len(occupancies))
        mean_occup_nterm = float(sum(occupancies_nterm)) / float(30)
        
        #logmean occup for this gene (and for first 10 aa)
        logmean_occup = math.exp(
            sum(map(math.log, occupancies)) / float(len(occupancies)))
        logmean_occup_nterm = math.exp(
            sum(map(math.log, occupancies_nterm)) / float(30))
        
        #create new rbs occupancy qualifiers for this region
        feature.qualifiers['mean_rbs_occupancy'] = mean_occup
        feature.qualifiers['logmean_rbs_occupancy'] = logmean_occup
        genes[gene_id]['mean_rbs_occupancy'] = mean_occup
        genes[gene_id]['logmean_rbs_occupancy'] = logmean_occup
        genes[gene_id]['mean_rbs_occupancy_nterm'] = mean_occup_nterm
        genes[gene_id]['logmean_rbs_occupancy_nterm'] = logmean_occup_nterm
        if 'rpkm' in feature.qualifiers:        
            genes[gene_id]['rpkm'] = feature.qualifiers['rpkm']
            genes[gene_id]['dpkm'] = feature.qualifiers['dpkm']
        
        #update the scores for each codon in this CDS
        feature_codons = enumerate(
            zip(range(0, len(feature)-3, 3), 
            range(3, len(feature), 3)))
        
        #for each codon in this feature, update occupancy and mRNA stats
        for codon_i, (codon_start, codon_end) in feature_codons:
        
            codon = feature_seq[codon_start:codon_end].tostring() 
            codon_occup = sum(occupancies[codon_start:codon_end+1]) / 4                       
            codons[codon]['count_total'] += 1
            codons[codon]['mean_occup'] += mean_occup
            codons[codon]['logmean_occup'] += logmean_occup
            
            #skip this if this feature does not have rna-seq data
            if 'rpkm' in feature.qualifiers:
                codons[codon]['rnaseq_count'] += 1
            
            #ignore the first and last 10 codons for occupancy measures
            if codon_start > 30 and (len(feature) - codon_end > 30):
                codons[codon]['count'] += 1
                codons[codon]['occup'] += codon_occup / mean_occup
                codons[codon]['mean_occup'] += mean_occup
                codons[codon]['logmean_occup'] += logmean_occup
                if 'rpkm' in feature.qualifiers:
                    codons[codon]['rpkm'] += feature.qualifiers['rpkm'] 
                    codons[codon]['dpkm'] += feature.qualifiers['dpkm']
                
            #but record the first 10 codons separately
            if codon_start < 30:
                codons[codon]['count_nterm'] += 1
                codons[codon]['occup_nterm'] += codon_occup / mean_occup
                codons[codon]['mean_occup_nterm'] += mean_occup
                codons[codon]['logmean_occup_nterm'] += logmean_occup
                if 'rpkm' in feature.qualifiers:
                    codons[codon]['rpkm_nterm'] += feature.qualifiers['rpkm'] 
                    codons[codon]['dpkm_nterm'] += feature.qualifiers['dpkm']
        
        #do the same thing for all sliding hexamers
        feature_hexamers = enumerate(
            zip(range(0, len(feature)), 
            range(6, len(feature)+6)))
        
        #for each codon in this feature, update occupancy and mRNA stats
        for hexamer_i, (hex_start, hex_end) in feature_hexamers:
        
            hexamer = feature_seq[hex_start:hex_end].tostring()   
            hex_occup = sum(occupancies[hex_start:hex_end+1]) / 7                     
            hexamers[hexamer]['count_total'] += 1
            hexamers[hexamer]['mean_occup'] += mean_occup
            hexamers[hexamer]['logmean_occup'] += logmean_occup
            
            #some features do not have rnaseq data
            if 'rpkm' in feature.qualifiers:
                hexamers[hexamer]['rnaseq_count'] += 1
            
            #ignore the first and last 10 codons for occupancy measures
            if hex_start > 25 and (len(feature) - hex_end > 30):
                hexamers[hexamer]['count'] += 1
                hexamers[hexamer]['occup'] += hex_occup / mean_occup
                hexamers[hexamer]['mean_occup'] += mean_occup
                hexamers[hexamer]['logmean_occup'] += logmean_occup
                if 'rpkm' in feature.qualifiers:
                    hexamers[hexamer]['rpkm'] += feature.qualifiers['rpkm'] 
                    hexamers[hexamer]['dpkm'] += feature.qualifiers['dpkm']
                
            #but record the first 10 codons separately
            if hex_start < 25:
                hexamers[hexamer]['count_nterm'] += 1
                hexamers[hexamer]['occup_nterm'] += codon_occup / mean_occup
                hexamers[hexamer]['mean_occup_nterm'] += mean_occup
                hexamers[hexamer]['logmean_occup_nterm'] += logmean_occup
                if 'rpkm' in feature.qualifiers:                
                    hexamers[hexamer]['rpkm_nterm'] += feature.qualifiers['rpkm'] 
                    hexamers[hexamer]['dpkm_nterm'] += feature.qualifiers['dpkm']
            
    #Next, let's compute the codon usage per amino acid
    codon_lookup = CodonUsageMemex(build_codon_usage_dict())
    aa_usage_counts = defaultdict(int)
    aa_usage_counts_nterm = defaultdict(int)
    
    #calc aa usages 
    for codon in codons.keys():
        aa = codon_lookup.codon_to_aa_dict[codon]
        aa_usage_counts[aa] += codons[codon]['count']
        aa_usage_counts_nterm[aa] += codons[codon]['count_nterm']
    
    # divide by aa usages to get frequencies, 
    # normalize occupancy and rna seq to count
    for codon, stats in codons.items():
        aa = codon_lookup.codon_to_aa_dict[codon]
        aa_total = aa_usage_counts[aa]
        aa_total_nterm = aa_usage_counts_nterm[aa]
        stats['freq'] = float(stats['count']) / float(aa_total)
        if aa_total_nterm > 0:
            stats['freq_nterm'] = (float(stats['count_nterm'])
                    / float(aa_total_nterm))
        else:
            stats['freq_nterm'] = 0
        
        stats['norm_occup'] = stats['occup'] / stats['count']
        stats['norm_occup_nterm'] = stats['occup_nterm'] / stats['count']
        
        stats['rpkm'] = stats['rpkm'] / stats['rnaseq_count']
        stats['dpkm'] = stats['dpkm'] / stats['rnaseq_count']
        stats['rpkm_nterm'] = stats['rpkm_nterm'] / stats['rnaseq_count']
        stats['dpkm_nterm'] = stats['dpkm_nterm'] / stats['rnaseq_count']
        
        stats['aa'] = aa
    
    #Finally, we need to rescale these into what Frydman calls the cu_i, 
    #which is the translational level rescaled to 1 per AA
    
    #get adjusted vlues relative to max
    get_max = lambda a_dict,a_key: max([v[a_key] for v in a_dict.values()] + 
            [sys.float_info[3]])
    max_m_occup = get_max(codons,'mean_occup')
    max_m_occup_nt = get_max(codons,'mean_occup_nterm')
    max_lm_occup = get_max(codons,'logmean_occup')
    max_lm_occup_nt = get_max(codons,'logmean_occup_nterm')
    max_drs_occup = get_max(codons,'dpkm')
    max_drs_occup_nt = get_max(codons,'dpkm_nterm')
    max_rs_occup = get_max(codons,'rpkm')
    max_rs_occup_nt = get_max(codons,'rpkm_nterm')

    for stats in codons.values():
        stats['mean_cu_i'] = stats['mean_occup'] / max_m_occup
        stats['logmean_cu_i'] = stats['logmean_occup'] / max_lm_occup
        stats['dpkm_cu_i'] = stats['dpkm'] / max_drs_occup
        stats['rpkm_cu_i'] = stats['rpkm'] / max_rs_occup
        
        stats['mean_cu_i_nt'] = stats['mean_occup_nterm'] / max_m_occup_nt
        stats['logmean_cu_i_nt'] = stats['logmean_occup_nterm'] / max_lm_occup_nt
        stats['dpkm_cu_i_nt'] = stats['dpkm_nterm'] / max_drs_occup_nt
        stats['rpkm_cu_i_nt'] = stats['rpkm_nterm'] / max_rs_occup_nt
        
    #get cTE scores and create nTE' scores
    cte_file = open(CTE_SCORE_FILE)
    for line in cte_file:
        (codon, cte) = line.split()
        codons[codon]['cte'] = float(cte)
    
    for codon, stats in codons.items():
        stats['mean_nte'] = stats['cte'] / stats['mean_cu_i']
        stats['logmean_nte'] = stats['cte'] / stats['logmean_cu_i']
        stats['drs_nte'] = stats['cte'] / stats['dpkm_cu_i']
        stats['rs_nte'] = stats['cte'] / stats['rpkm_cu_i']
        
        stats['mean_nte_nt'] = stats['cte'] / (stats['mean_cu_i_nt'] + sys.float_info[3])
        stats['logmean_nte_nt'] = stats['cte'] / (stats['logmean_cu_i_nt'] + sys.float_info[3])
        stats['drs_nte_nt'] = stats['cte'] / (stats['dpkm_cu_i_nt'] + sys.float_info[3])
        stats['rs_nte_nt'] = stats['cte'] / (stats['rpkm_cu_i_nt'] + sys.float_info[3])
    
    #get adjusted nTE' out of max to find nTE
    max_m_nte = get_max(codons,'mean_nte')
    max_lm_nte = get_max(codons,'logmean_nte')
    max_rs_nte = get_max(codons,'rs_nte')
    max_drs_nte = get_max(codons,'drs_nte')

    max_m_nte_nt = get_max(codons,'mean_nte_nt')
    max_lm_nte_nt = get_max(codons,'logmean_nte_nt')
    max_rs_nte_nt = get_max(codons,'rs_nte_nt')
    max_drs_nte_nt = get_max(codons,'drs_nte_nt')
    
    
    for stats in codons.values():
        stats['mean_nte'] = stats['mean_nte'] / max_m_nte
        stats['logmean_nte'] = stats['logmean_nte'] / max_lm_nte
        stats['rs_nte'] = stats['rs_nte'] / max_rs_nte
        stats['drs_nte'] = stats['drs_nte'] / max_drs_nte
        
        stats['mean_nte_nt'] = stats['mean_nte_nt'] / max_m_nte_nt
        stats['logmean_nte_nt'] = stats['logmean_nte_nt'] / max_lm_nte_nt
        stats['rs_nte_nt'] = stats['rs_nte_nt'] / max_rs_nte_nt
        stats['drs_nte_nt'] = stats['rs_nte_nt'] / max_rs_nte_nt
    
    return (genome, codons, genes, hexamers, bad_cds)
    
def get_frame_occupancies(genome, profile, border = 15):
    
    frame_counts = defaultdict(int)
    
    
    for feature in genome.features:
        #only take CDS features:
        if feature.type != 'CDS': continue

        #skip if no gene id
        try: gene_id = feature.qualifiers['gene']
        except(KeyError): continue

        #Li/Weissman don't use these genes due to frameshift/homology
        if gene_id in ('tufA','tufB','prfX','dnaX'): continue
        
        codon_data = get_feature_occupancy(genome, profile, feature, False, 0)
        
        blank_local_codons = [0,0,0]
        
        #calculate the per_frame occupancies
        for line in codon_data:                    
            blank_local_codons[line[3]] += line[4]
            if line[3] == 2:
                for i in (0,1,2):
                    frame_counts[i] += blank_local_codons[i] / (float(sum(blank_local_codons)+1))
        
        
    return frame_counts

def get_feature_occupancy(genome, profile, feature, 
        to_table=False, border=15):
    '''
    This function creates a nice per base table of occupancy,
    and lines it up with the nucleotide, the amino acid frame, and can give
    borders before and after. It was useful in lining up the A-site offset
    and finding interesting features like the SecM pause.
    '''

    strand = 'fwd' if feature.strand == 1 else 'rev'

    region_profile = extract_region_profile(
        profile, 
        strand, 
        feature.location.start.position-border, 
        feature.location.end.position+border)

    feature_location_slice = slice(
        feature.location.start.position-border,
        feature.location.end.position+border)
    feature_seq = genome.seq[feature_location_slice]
    
    if strand == 'rev': feature_seq = feature_seq.reverse_complement()

    nt_range = range(-border,len(feature_seq.tostring())+border)
    aa_range = [math.floor(float(1 + (i/3))) for i in nt_range]

    codon_data = list(itertools.izip(
        feature_seq.tostring(), 
        nt_range, 
        aa_range, 
        itertools.cycle([0,1,2]), 
        region_profile.values()))

    #print a nice tab delimited table
    if to_table:
        print '\t'.join(('Base','NT','AA','Frame','Occupancy'))
        for line in codon_data:
            if line[1] == '0':
                print '<START>'
            print '\t'.join([str(i) for i in line])
            if line[1] == len(feature)-1:
                print '<END>'
    
    return codon_data


def test_secM_occupancy(genome, profile):
    for feature in genome.features:
        if 'gene' in feature.qualifiers and feature.qualifiers['gene'][0] == 'secM':
            break
    return feature
            

if __name__ == '__main__':

    usage_string =  '''
    
    Scripts for loading ribosome profile data and annotating codons and
    gene features using this data.
    
    Usage:
    
    Creates profile from bam file(s):
        weissman_profile create file1.bam .. fileN.bam profile_output
    
    Merges profiles into one:
        weissman_profile create profile1 ... profileN profile_output
    
    Print per-codon/pre-gene/per-hex gene translation table:
        This creates three files, one with per-codon data, another with
        per-CDS rpkm, log-mean and mean occupancy, and a third with
        the codon-type data but per hexanucleotide, in all frames.
    
        weissman_profile codons genome.gbk profile codon_output gene_output
    '''
        
    if len(sys.argv) == 1:    
        print usage_string
        
    elif sys.argv[1] == 'create':
        bam_files = sys.argv[2:-1]
        profile_output = sys.argv[-1]
        save_profile(create_genomewide_profile(bam_files, profile_output),
            profile_output)
    
    elif sys.argv[1] == 'merge':
        profile_files = sys.argv[2:-1]
        profile_output = sys.argv[-1]
        save_profile(merge_all_profiles(profile_files), profile_output)
    
    elif sys.argv[1] == 'codons':
        
        if len(sys.argv) != 7: 
            print 'Invalid usage.'
            print usage_string
        
        (genome_file, profile_file, codon_out_file, gene_out_file, 
            hex_out_file) = sys.argv[2:]
        
        #load genome
        genome = SeqIO.read(genome_file, 'genbank')
        
        #load pickled per base ribosome profile
        profile = load_profile(profile_file)
        
        #add mRNA data
        genome = get_per_gene_mrna_data(genome)
        
        #calculate per codon, genome, and gene based data
        translation_info = add_translation_info(genome, profile)
        (genome, codons, genes, hexamers, bad_cds) = translation_info

        #print genes that had problems to stderr
        for gene_err in bad_cds.items():
            print >> sys.stderr, '\t'.join(gene_err)
        
        #print codon data
        with open(codon_out_file, 'w') as codon_fh:
            #print codon data
            print >> codon_fh, '\t'.join(['codon'] + codons['ATG'].keys())
            #header
            #codon rows
            for codon, stats in codons.items():
                print >> codon_fh, '\t'.join([codon] 
                        + map(str, stats.values()))
        
        #print gene data
        with  open(gene_out_file, 'w') as gene_fh:
        #header
            print >> gene_fh, '\t'.join(['gene'] + codons.values()[0].keys())
            #gene rows
            for gene, stats in genes.items():
                print >> gene_fh, '\t'.join([gene] 
                        + map(str, stats.values()))
        
        #print hex data
        with  open(hex_out_file, 'w') as hex_fh:
        #header
            print >> hex_fh, '\t'.join(['hexamer'] + hexamers.values()[0].keys())
            #gene rows
            for hexamer, stats in hexamers.items():
                print >> hex_fh, '\t'.join([hexamer] 
                        + map(str, stats.values()))
    
    else:
        print usage_string
        
        




