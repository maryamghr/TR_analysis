vcf_dir=/confidential/tGenVar/vntr/output_maryam/tools/run_all_tools/output/hgsvc/TandemTwist/asm/
sample_info_file=/confidential/tGenVar/vntr/output_maryam/tools/run_all_tools/input_folder/samples/hgsvc/HGSVC2_sample_info.tsv
pop_file=output/pop_efficient/population_df.csv
cumm_CN_file=output/pop_efficient/cumulative_CN_allele_classification.csv

log_cum_CN=output/log/pop_efficient/cumulative_CN_allele_classification.log

time python3 -m run_population_analysis compute-population-df --vcf_dir $vcf_dir --sample_info_file $sample_info_file --output_pop_df $pop_file

# time python3 -m run_population_analysis comm-allele-classification --population_df_file $pop_file \
    # --sample_info_file $sample_info_file --col_prefix CN --output_cum_allele $cumm_CN_file