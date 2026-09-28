## pyrho 
```bash
module load BCFtools/1.22-GCC-12.3.0
module load Miniconda3/23.10.0-1
source activate /nesi/nobackup/uoa04053/tram_ont/tools/pyrho_env

outdir=/nesi/nobackup/ga03714/tram/recombination/pyrho/Heliocidaris_tuberculata_mu5
mkdir -p $outdir
cd $outdir

i=$SLURM_ARRAY_TASK_ID
chr=$(sed "${i}q;d" ${outdir}/chrom_list.txt)


chrdir=${outdir}/${chr}
mkdir -p $chrdir
cd $chrdir

## remove Rangitahua and related samples for Ht
# cd /nesi/nobackup/ga03714/tram/wgs_vcf
# tabix Ht_filt_uniqID_noRel_HWDrecomb.vcf.gz
# bcftools query -l Ht_filt_uniqID_noRel_HWDrecomb.vcf.gz > Ht_filt_uniqID_noRel_HWDrecomb.ind

bgvcf=/nesi/nobackup/ga03714/tram/wgs_vcf/Ht_filt_uniqID_noRel_HWDrecomb.vcf.gz
bgvcf_id=/nesi/nobackup/ga03714/tram/wgs_vcf/Ht_filt_uniqID_noRel_HWDrecomb.ind



n_ploids=120
N_size=150

pop=Ht_rep5
pop_id=${pop}.ind
if [ ! -f $pop_id ]; then
    echo "shuf $bgvcf_id | head -n 60 > $pop_id"
fi

mu_rate=8.6e-9
# mu_rate=5e-9
# mu_rate=1e-9
# mu_rate=9.9e-9 ## test for Cr


## Demographic inference
module load Apptainer/1.2.5
### convert vcf to smc
samples=$(cat $pop_id | paste -sd ,)

if [ ! -f "${pop}_${chr}.smc.gz" ]; then
    singularity run /nesi/nobackup/uoa04053/tram_ont/tools/smcpp.sif \
    vcf2smc --cores $SLURM_CPUS_PER_TASK $bgvcf ${pop}_${chr}.smc.gz ${chr} ${pop}:${samples}
fi

# ### estimate pop size
if [ ! -f ${pop}_${chr}_mu${mu_rate}.final.json ]; then
    singularity run /nesi/nobackup/uoa04053/tram_ont/tools/smcpp.sif \
    estimate --cores $SLURM_CPUS_PER_TASK -o ./${pop}_${mu_rate} --base ${pop}_${chr}_mu${mu_rate} ${mu_rate} ${pop}_${chr}.smc.gz
    mv ./${pop}_${mu_rate}/${pop}_${chr}_mu${mu_rate}.final.json .
    rm -rf ./${pop}_${mu_rate}
fi

### plot inference
if [ ! -f ${pop}_${chr}_mu${mu_rate}_smcpp.csv ]; then
    singularity run /nesi/nobackup/uoa04053/tram_ont/tools/smcpp.sif \
    plot --cores $SLURM_CPUS_PER_TASK --csv ${pop}_${chr}_mu${mu_rate}_smcpp.pdf ${pop}_${chr}_mu${mu_rate}.final.json
fi

#################################
#################################
## Recombination estimate
if [ ! -f ${pop}_${chr}_mu${mu_rate}_lookuptable.hdf ]; then
    pyrho make_table --samplesize $n_ploids --approx  --moran_pop_size $N_size \
    --numthreads $SLURM_CPUS_PER_TASK --mu ${mu_rate} --outfile ${pop}_${chr}_mu${mu_rate}_lookuptable.hdf --smcpp_file ${pop}_${chr}_mu${mu_rate}_smcpp.csv
fi

if [ ! -f ${pop}_${chr}_mu${mu_rate}_hyperparam.tsv ]; then
    pyrho hyperparam --samplesize $n_ploids \
    --numthreads $SLURM_CPUS_PER_TASK \
    --ploidy 2 \
    --tablefile ${pop}_${chr}_mu${mu_rate}_lookuptable.hdf \
    --mu ${mu_rate} \
    --blockpenalty 25,50,75,100 \
    --windowsize 25,50,75,100 \
    --num_sims 10 \
    --smcpp_file ${pop}_${chr}_mu${mu_rate}_smcpp.csv \
    --outfile ${pop}_${chr}_mu${mu_rate}_hyperparam.tsv
fi

opt_block=$(awk -F'\t' 'NR==2 || $12 < min { min=$12; val=$1 } END { print val }' ${pop}_${chr}_mu${mu_rate}_hyperparam.tsv | sed "s/\\.0//")
opt_window=$(awk -F'\t' 'NR==2 || $12 < min { min=$12; val=$2 } END { print val }' ${pop}_${chr}_mu${mu_rate}_hyperparam.tsv | sed "s/\\.0//")

if [ ! -f ${pop}_${chr}.vcf.gz ]; then
    bcftools view -i 'MAF[0]>0.05' -r ${chr} -S $pop_id -Oz -o ${pop}_${chr}.vcf.gz  $bgvcf
fi

if [ ! -f ${pop}_${chr}_mu${mu_rate}.rmap ]; then
    pyrho optimize --tablefile ${pop}_${chr}_mu${mu_rate}_lookuptable.hdf \
    --numthreads $SLURM_CPUS_PER_TASK \
    --ploidy 2 \
    --vcffile ${pop}_${chr}.vcf.gz \
    --outfile ${pop}_${chr}_mu${mu_rate}.rmap \
    --blockpenalty $opt_block --windowsize $opt_window \
    --logfile .
fi
```

## Calculate mean pyrho
```bash
library(data.table)
library(tidyverse)

# for Heliocidaris tuberculata
## recombination estimation for NZ and OZ
rmap_files <- c(list.files("/nesi/nobackup/ga03714/tram/recombination/pyrho/Heliocidaris_tuberculata_mu5",
                         "Ht.*_mu(1|5).*.rmap", recursive = T, full.names = T),
                list.files("/nesi/nobackup/ga03714/tram/recombination/pyrho/Heliocidaris_tuberculata",
                           "Ht.*_mu8.*.rmap", recursive = T, full.names = T))
chrom <- unique(basename(dirname(rmap_files)))

### average rate across replicates for each chrom
rec_chrom <- do.call(rbind, lapply(chrom, function(c) {
  chrom_rmap_files <- grep(c, rmap_files, value = T)
  rmap_rep <- do.call(rbind, lapply(chrom_rmap_files, function(f) {
    rmap <- fread(f)
    rec_int <- rmap$V2 - rmap$V1
    rec <- sum(rmap$V3 * rec_int) / sum(rec_int) * 10^8
    rec_df <- data.frame(rec = rec,
                         length = sum(rec_int),
                         chrom = c,
                         mu = gsub(".+_mu|.rmap", "", basename(f)))
    return(rec_df)
  }))
  return(rmap_rep)
}))

# rec_chrom$mu <- c("1e-9", "6e-9", "8.6e-9")

# rec_chrom %>%
#   pivot_wider(names_from = "mu", values_from = "rec") %>%
#   View()

rec_chrom_avg <- rec_chrom %>%
  group_by(chrom, length, mu) %>%
  summarise(rec_avg = mean(rec))

rec_chrom_avg %>% 
  pivot_wider(names_from = mu, values_from = rec_avg) %>%
  write.csv(., "/nesi/nobackup/ga03714/tram/recombination/pyrho/pyrho_average_Ht.csv",
            quote = F, row.names = F)

### average rate genome-wide
# sum(rec_chrom_avg$rec_avg*rec_chrom_avg$length)/sum(rec_chrom_avg$length)

rec_chrom_avg %>%
  mutate(s = length*rec_avg) %>%
  group_by(mu) %>%
  dplyr::summarise(gw_rec = sum(s)/sum(length))



### fine-scale rate (average across replicates)

rec_map <- do.call(rbind, lapply(chrom, function(c) {
  chrom_rmap_files <- grep(c, rmap_files, value = T)
  rmap_rep <- do.call(rbind, lapply(chrom_rmap_files, function(f) {
    rmap <- fread(f) %>%
      mutate(chr = c,
             start = V1,
             end = V2,
             r = V3 * 10^8) %>%
      select(chr, start, end, r)
    return(rmap)
  }))
  rmap_rep <- rmap_rep %>%
    group_by(chr, start, end) %>%
    mutate(rate = mean(r))
  return(rmap_rep)
}))

## average rate genome-wide


rec_map_wd <- rec_map %>%
  group_by(chr) %>%
  mutate(bin = cut_interval(start, length = 0.3e6)) %>%
  group_by(chr, bin) %>%
  summarise(avg_rate = mean(rate),
            pos = median(start))

outlier_snps <- read.table("/nesi/nobackup/ga03714/Ht_raw_lcWGS_data/Ht_pop_analysis/LDprune_data/Ht_combined_outliers_1.txt", sep = "_", header = F, col.names = c("chr", "position"))
outlier_snps <- outlier_snps %>%
  group_by(chr) %>%
  mutate(bin = cut_interval(position, length = 0.3e6)) %>%
  group_by(chr, bin) %>%
  summarise(n_out = n(),
            pos = median(position))


ggplot(rec_map_wd, aes(pos/1e6, avg_rate)) +
  facet_wrap(vars(chr), scales = "free") +
  geom_line() +
  geom_point() +
  # geom_point(data = outlier_snps, aes(x = pos/1e6, y = n_out), color ="red") +
  theme_minimal() +
  xlab("Position (Mb)") +
  ylab("Recombination rate (cM/Mb) per 300kb window")
```
