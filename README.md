# 1KGP Sweeps Enformer Project Steps

→ Path to project directory: `/gpfs/gibbs/project/reilly/dal83/1KGP_sweeps_enformer_project`

## Setting up environment

1. `**ssh` into your YCRC McCleary cluster account**
2. **Set up a conda environment for your project**
    - [https://docs.ycrc.yale.edu/clusters-at-yale/guides/conda/](https://docs.ycrc.yale.edu/clusters-at-yale/guides/conda/)
    - Recommend setting up environment with `python`, `biopython`, `numpy`, `samtools`, `bcftools`, `tabix`, `r`, `rstudio-desktop`, `seaborn` , `pandas`, `matplotlib`, `jupyter`, `cudatoolkit=11.8.0`

## Retrieving FASTA sequences

1. **Download indexed VCF data (run a slurm job using `wget` and `sbatch`)**
    - I used the high coverage data from 1KGP: [http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/)
2. `**tabix` the region of interest: `tabix -h [.vcf.gz file] [chromosome and sequence range] >> [location for .vcf file]`**
    - eg. `tabix -h [http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz) chr2:134851076-136851076 >> ../complete_vcf_high_coverage/high_cov_2Mb_LCT_region.vcf`
3. **Download reference fasta file from UCSC (should be in format `chr#.fa.gz`)**
    1. Download from: [https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/)
    2. `gunzip` the fasta file
    3. Extract region of interest using `samtools faidx`
        1. `samtools faidx [fasta file] [chr#: region]`
            1. eg. `samtools faidx ../chr2.fa chr2:135608645-137608646 > ref_genome_LCT_2Mb.fa`
        2. for help: [https://manpages.ubuntu.com/manpages/focal/man1/samtools-faidx.1.html](https://manpages.ubuntu.com/manpages/focal/man1/samtools-faidx.1.html)
4. **Create FASTA files for each sample given in the VCF file.**
    - Run Python script `[merge_vcf_fasta.py](https://www.notion.so/merge_vcf_fasta-py-651eb755d7994507956818d51ef3af53?pvs=21)`
        
        [merge_vcf_fasta.py](https://www.notion.so/merge_vcf_fasta-py-651eb755d7994507956818d51ef3af53?pvs=21)
        
    - Useful bcftools documentation: [https://samtools.github.io/bcftools/bcftools.html](https://samtools.github.io/bcftools/bcftools.html)
    - Useful Biopython/SeqIO documentation: [http://biopython.org/DIST/docs/tutorial/Tutorial.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf)

## Prepare FASTA files for Enformer

1. **Run Python script `[sequence_segments.py](https://www.notion.so/sequence_segments-py-803693cc54564691952a838effa8be2b?pvs=21)` to split FASTA entries into overlapping regions of 393,216 bp.**
    
    [sequence_segments.py](https://www.notion.so/sequence_segments-py-803693cc54564691952a838effa8be2b?pvs=21)
    
    → *FASTA file entries need to be split up into smaller entries because Enformer takes in sequences 393,216 bp long.*
    
    → I stored the overlapping fasta regions in `/data/hg38_fasta_files/segmented_seq_chr2.fa`
    

## Running Enformer

1. **Download Enformer colab notebook.**
    - [https://colab.research.google.com/github/deepmind/deepmind_research/blob/master/enformer/enformer-usage.ipynb#scrollTo=Q48earqRyFa6](https://colab.research.google.com/github/deepmind/deepmind_research/blob/master/enformer/enformer-usage.ipynb#scrollTo=Q48earqRyFa6)
2. **Request a GPU (A100 recommended) on the cluster.**
    - `salloc -C "a100" --gpus=1 --time=6:00:00 --partition gpu`
    - More info: [https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/resource-requests/#request-gpus](https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/resource-requests/#request-gpus)
3. **Create a TensorFlow conda environment.**
    
    ```bash
    # This creates a new conda environment from scratch.
    # Just activate your conda environment if you already established one.
    module load miniconda
    conda create --name tf-condacuda python numpy pandas matplotlib jupyter cudatoolkit=11.8.0
    conda activate tf-condacuda
    
    pip install nvidia-cudnn-cu11==8.6.0.163
    
    # Store system paths to cuda libraries for gpu communication
    mkdir -p $CONDA_PREFIX/etc/conda/activate.d
    echo 'CUDNN_PATH=$(dirname $(python -c "import nvidia.cudnn;print(nvidia.cudnn.__file__)"))' >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib/:$CUDNN_PATH/lib' >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
    
    # Install tensorflow
    pip install tensorflow==2.12.*
    ```
    
    - More info: [https://docs.ycrc.yale.edu/clusters-at-yale/guides/gpus-cuda/](https://docs.ycrc.yale.edu/clusters-at-yale/guides/gpus-cuda/)
4. `**pip install` or `conda install` remaining packages (listed in the Enformer Python script).**
    - All my dependencies if you have trouble:
        - [Enformer Dependencies](https://www.notion.so/Enformer-Dependencies-711842e61c2b4d2c86ecaddfa78c5080?pvs=21)
    - *Fixing potential cuDNN error:* `conda install cudnn=8.8.0`
    - *Fixing potential “Start cannot spawn child process” error:* `conda install -c nvidia cuda-nvcc`
5. **Modify Enformer script**
    1. Set `fasta_file` variable to your fasta.
    2. I added a function to read in the fasta file entries and make predictions. You can delete most functions, classes, and variables. See my modified script:
        
        [enformer_usage.py](https://www.notion.so/enformer_usage-py-32c43ee98448435db6deeab4def7fe8b?pvs=21)
        
        → *Make sure sequences passed into Enformer are all uppercase (using `.upper()` method on a string).*
        
    3. In `enformer_usage.py`, I only looked at one 393,216bp region at a time, but this script can easily be adjusted to provide predictions for all 393,216bp regions in your fasta files.
6. **If you experience the error “ValueError: Trying to load a model of incompatible/unknown type”: `export TFHUB_CACHE_DIR=$HOME/palmer_scratch/tfhub_modules`**
    1. If this works, then put it into your `.bashrc` file (`nano $HOME/.bashrc`, then add this line to the end) so that it gets set automatically.
7. **At this point, you should have Enformer predictions for all of the fasta sequences stored in separate compressed npz files!**

## Data Analysis in R

1. **I recommend using RStudio Desktop through McCleary OnDemand**
    - [https://ood-mccleary.ycrc.yale.edu/pun/sys/dashboard/batch_connect/sys/ycrc_rstudio_conda_r/session_contexts/new](https://ood-mccleary.ycrc.yale.edu/pun/sys/dashboard/batch_connect/sys/ycrc_rstudio_conda_r/session_contexts/new)
2. **Install necessary R packages: `gplots`, `reticulate`, `readr`, `tidyverse`,`ggridges` by either:**
    - `conda install r-[package]` on previous conda environment
    - create a new conda environment: `conda create -c conda-forge --name r-env r-base r-essentials r-gplots r-reticulate r-readr r-tidyverse r-ggridges`
3. **Run script to put predictions into one table and normalize values.**
    
    [normalize_preds.R](https://www.notion.so/normalize_preds-R-6634fb63138645e5b01f19ac345e613d?pvs=21)
    
4. **Create heatmap.**
    - I created a separate R file to do this, but you could save a few steps by just doing this in the normalize_preds.R file.
        - At this point, the row names and column names in the matrix might be empty, so make sure to fill those in with the prediction tracks and human samples.
            
            [heatmap.R](https://www.notion.so/heatmap-R-636d8cae07d7462ebec548b50b149945?pvs=21)
            
5. **I changed the x-axis labels based on sample populations. Use `heatmap.2()` to specify label colors.**
    
    [heatmap.R [changing axis colors]](https://www.notion.so/heatmap-R-changing-axis-colors-569a5100add94c73815470435b7444eb?pvs=21)
    

Here are some examples of heatmaps I created with subsets of the normalized prediction data:

![*Enformer_Subset_200_700_Alt_Ref.pdf.* Red = alt allele, black = ref allele.]
(<img width="897" alt="Screen_Shot_2023-07-28_at_3 58 09_PM" src="https://github.com/dal83/1KGP-Sweeps-Enformer-Project/assets/126291855/ad169182-6a45-4ffb-bacf-8dea2bd2a8cb">)

*Enformer_Subset_200_700_Alt_Ref.pdf.* Red = alt allele, black = ref allele.

![*Enformer_Subset_200_700_Pops.pdf.* Yellow = Europe, Green = South Asia, Blue = Americas]
(<img width="932" alt="Screen_Shot_2023-07-28_at_3 58 21_PM" src="https://github.com/dal83/1KGP-Sweeps-Enformer-Project/assets/126291855/dc487f3c-d468-4516-8027-9d10ec531f1b">)
*Enformer_Subset_200_700_Pops.pdf.* Yellow = Europe, Green = South Asia, Blue = Americas

## Looking at Correlations Between Enformer Tracks

1. I went through the above steps to retrieve Enformer predictions and a normalized prediction matrix for the prediction window surrounding the transcription start site (located in my `/data/LCT_transcription_start_site/` directory).
2. I then created a heatmap to view correlations between the DNASE predictions for the window containing the LCT causal SNP and the CAGE predictions for the transcription start site window.
    1. Only looking at reference LCT allele samples to compare haplotype background (filtered out the samples that contain the alt LCT allele here).
    2. I used `cor()` in R to do this. See [correlation_map.R](https://www.notion.so/correlation_map-R-d706744b31594240beac21562fd30e49?pvs=21).
        
        ![*DNASExCAGE_Predictions.png*. Blue = positive values, red = negative values.]
       ![DNASExCAGE_Predictions](https://github.com/dal83/1KGP-Sweeps-Enformer-Project/assets/126291855/3399b87a-b819-4e9b-8c60-c19387d25bdf)

        
        *DNASExCAGE_Predictions.png*. Blue = positive values, red = negative values.
        
        [correlation_map.R](https://www.notion.so/correlation_map-R-d706744b31594240beac21562fd30e49?pvs=21)
