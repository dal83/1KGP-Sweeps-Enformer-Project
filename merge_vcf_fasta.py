"""
* Filename: merge_vcf_fasta.py
* Part of: 1KGP_sweeps_enformer_project
* Contains: Script to create FASTA files based on VFC data.
* Debbie Lilly, July 2023
"""

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import vcf
import os

# Read in FASTA file sequence (change to your FASTA file)
ref_genome = "../data/hg38_fasta_files/ref_genome_high_cov_LCT_2Mb.fa"
# Establishing FASTA file to write output to (change to your FASTA output file)
file_out = "fasta_files/complete_genome_records_hc2.fa"
# Clearing contents of the FASTA output file to ensure it's empty
with open(file_out, 'w') as f_out:
    SeqIO.write('', f_out, 'fasta')


# Open VCF file
vcfFullPath = '../data/complete_vcf_high_coverage/high_cov_2Mb_LCT_region.vcf'  # Replace with your VCF file
vcfRecords = vcf.Reader(open(vcfFullPath, 'r'))
allSampleNames = vcfRecords.samples

# Inputs: Allele (either "1" or "2"), vcf file specific to given allele of given sample.
# Returns: None.
# Does: Uses consensus to create unique FASTA file based on given vcf file.
#       Appends this FASTA entry to FASTA output file file_out.
#       ID for each FASTA entry: chr#:startPosition-endPosition_sampleName_allele.
#           e.g.: chr2:35600000-37600000_HG00096_1
# Note: I temporarily store each FASTA file in a separate FASTA file before
#       appending it to file_out. There's probably a better way to do this.
def consensus_proc(allele, sampVcf):
    temp_fa_path = "fasta_files/temp_fasta2.fa" # Replace with your temporary FASTA file location.
    os.system("touch " + temp_fa_path)
    os.system("bcftools consensus -f " + ref_genome + " --haplotype " + allele + " " +
              sampVcf + ".gz > " + temp_fa_path)
    temp_genome = SeqIO.read(temp_fa_path, "fasta")
    fasta_fields = SeqRecord(
        temp_genome.seq, id=(temp_genome.id + "_" + samp + "_" + allele), description="")
    with open(file_out, 'a') as f_out:
        SeqIO.write(fasta_fields, f_out, 'fasta')


# Parsing through each sample name in the VCF file. Creates a unique VCF file for each sample.
for samp in allSampleNames:
    # Replace with your directory for storing vcf files.
    vcfSampPath = '../results/vcf_samps_high_cov/' + samp + '.vcf'
   
   # Filtering out non-SNP variants, multiallelic variants, and any entries for which no samples contain alt allele.
   # New VCF file will only contain variant info for current sample (samp).
    os.system("bcftools view -Ou --samples " + samp + " --types snps " +
              vcfFullPath + " | bcftools view -Ou -i \'GT=\"alt\"\' | bcftools norm --multiallelics - > " + vcfSampPath)
    os.system("bgzip --force " + vcfSampPath)
    os.system("tabix -p vcf " + vcfSampPath + ".gz")

    # Creating unique FASTA files for each allele of given samp.
    consensus_proc("1", vcfSampPath)
    print(samp + "allele 1 complete")
    consensus_proc("2", vcfSampPath)
    print(samp + "allele 2 complete")
