vcf_dir=/confidential/tGenVar/vntr/output_maryam/tools/run_all_tools/output/hgsvc/TandemTwist/asm/
sample_info_file=/confidential/tGenVar/vntr/output_maryam/tools/run_all_tools/input_folder/samples/hgsvc/HGSVC2_sample_info.tsv
pop_file=output/pop_efficient/population_df.csv
outdir_cumm_allele=output/pop_efficient/cumm_allele


log_cum_CN=output/log/pop_efficient/cumulative_CN_allele_classification.log
log_cum_len=output/log/pop_efficient/cumulative_length_allele_classification.log
log_cum_seq=output/log/pop_efficient/cumulative_sequence_allele_classification.log

# time python3 -m run_population_analysis compute-population-df --vcf_dir $vcf_dir --sample_info_file $sample_info_file --output_pop_df $pop_file

(time python3 -m run_population_analysis comm-allele-classification --population_df_file $pop_file \
    --sample_info_file $sample_info_file --col_prefix CN --output_dir $outdir_cumm_allele \
    > $log_cum_CN 2>&1) &

(time python3 -m run_population_analysis comm-allele-classification --population_df_file $pop_file \
    --sample_info_file $sample_info_file --col_prefix len --output_dir $outdir_cumm_allele \
    > $log_cum_len 2>&1) &

(time python3 -m run_population_analysis comm-allele-classification --population_df_file $pop_file \
    --sample_info_file $sample_info_file --col_prefix seq --output_dir $outdir_cumm_allele \
    > $log_cum_seq 2>&1) &