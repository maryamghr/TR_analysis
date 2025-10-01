import argparse
import os
import pandas as pd
import importlib
from pathlib import Path
import population_analysis_efficient
importlib.reload(population_analysis_efficient)
from population_analysis_efficient import get_sample_info_df, compute_population_df, count_alleles, plot_allele_count, compute_cumulative_allele_classification, plot_cumulative_allele_classification


def main():
    parser = argparse.ArgumentParser(
        description="Run functions from population_analysis via subcommands."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    
    # Subcommand: compute_population_df
    parser_pop = subparsers.add_parser("compute-population-df", help="Compute population dataframe")
    parser_pop.add_argument("--vcf_dir", help="Input directory", required=True)
    parser_pop.add_argument("--sample_info_file", help="Sample info file (CSV/TSV)", required=True)
    parser_pop.add_argument("--cores", type=int, default=1, help="Number of cores to use")
    parser_pop.add_argument("--output_pop_df", default="output/population/population_df.csv", help="Output CSV path")
    #parser_pop.add_argument("count_alleles", help="Count alleles in dataframe", default=False)
    #parser_pop.add_argument("output_allele_counts", help="output directory for allele counts", default="output/population/allele_counts")
    #parser_pop.add_argument("output_allele_counts_plot", help="output file for allele counts plot", default="output/population/multi_facet_allele_barplots.pdf")

    # Subcommand: count_alleles
    parser_count = subparsers.add_parser("count-alleles", help="Count alleles in dataframe")
    parser_count.add_argument("--population_df_file", help="CSV with populationdataframe", required=True)
    parser_count.add_argument("--sample-size", type=int, default=70, help="Sample size")
    parser_count.add_argument("output_allele_counts", help="output directory for allele counts")
    parser_count.add_argument("output_allele_counts_plot", help="output file for allele counts plot")

    # Subcommand: commulative-allele-classification
    parser_cum = subparsers.add_parser("comm-allele-classification", help="Compute commulative allele classification")
    parser_cum.add_argument("--population_df_file", help="CSV with populationdataframe", required=True)
    parser_cum.add_argument("--sample_info_file", help="Sample info file (CSV/TSV)", required=True)
    parser_cum.add_argument("--col_prefix", help="Column prefix for population dataframe (CN, len, seq)", default="CN")
    parser_cum.add_argument("--output_cum_allele", help="Output file for commulative allele classification", required=True)
    # parser_cum.add_argument("--output_cum_allele_plot", help="Output file for commulative allele classification plot", required=True)

    args = parser.parse_args()
    

    if args.command == "compute-population-df":
        sample_info = get_sample_info_df(args.sample_info_file)
        pop_df = compute_population_df(args.vcf_dir, output=args.output_pop_df)
        pop_df.to_csv(args.output_pop_df, index=False, sep='\t')
        print(f"Saved population_df to {args.output_pop_df}")

        # if args.count_alleles:
        #     CN_allele_counts = count_alleles(pop_df, column="CN", sample_size=args.sample_size, output=os.path.join(args.output_allele_counts, "CN_allele_counts.csv"))
        #     len_allele_counts = count_alleles(pop_df, column="len_", sample_size=args.sample_size, output=os.path.join(args.output_allele_counts, "len_allele_counts.csv"))
        #     motif_allele_counts = count_alleles(pop_df, column="motifs", sample_size=args.sample_size, output=os.path.join(args.output_allele_counts, "motif_allele_counts.csv"))
        #     plot_allele_count(CN_allele_counts, len_allele_counts, motif_allele_counts, output=args.output_allele_counts_plot)

    elif args.command == "count-alleles":
        pop_df = pd.read_csv(args.population_df_file, sep='\t')
        CN_allele_counts = count_alleles(pop_df, column="CN", sample_size=args.sample_size, output=os.path.join(args.output_allele_counts, "CN_allele_counts.csv"))
        len_allele_counts = count_alleles(pop_df, column="len_", sample_size=args.sample_size, output=os.path.join(args.output_allele_counts, "len_allele_counts.csv"))
        motif_allele_counts = count_alleles(pop_df, column="motifs", sample_size=args.sample_size, output=os.path.join(args.output_allele_counts, "motif_allele_counts.csv"))
        plot_allele_count(CN_allele_counts, len_allele_counts, motif_allele_counts, output=args.output_allele_counts_plot)
    
    elif args.command == "comm-allele-classification":
        pop_df = pd.read_csv(args.population_df_file, sep='\t')
        print(pop_df.head())
        sample_info=get_sample_info_df(args.sample_info_file)
        cumm_df = compute_cumulative_allele_classification(pop_df, sample_info, col_prefix=args.col_prefix, output=args.output_cum_allele)
        # plot_cumulative_allele_classification(cum_df, output=args.output_cum_allele_plot)


if __name__ == "__main__":
    main()

# vcf_dir="/confidential/tGenVar/vntr/output_maryam/tools/run_all_tools/output/hgsvc_test/TandemTwist/asm/"
# sample_info_file="/confidential/tGenVar/vntr/output_maryam/tools/run_all_tools/input_folder/samples/hgsvc/HGSVC2_sample_info.tsv"
# output_pop_df='output/population/population_df_test.csv'
# pop_df_file='output/population/population_df_test.csv'
# output_allele_counts='output/population/allele_counts'
# output_allele_counts_plot='output/population/multi_facet_allele_barplots.pdf'
# out_cum_allele=
# output_cum_allele_plot=
# time python run_population_analysis.py compute_population_df --vcf_dir $vcfdir --sample_info_file $sampleinfo --cores 12 --output $output
# time python run_population_analysis.py comm-allele-classification --population_df_file $popdf --sample_info_file $sampleinfo --output_cum_allele $outputcum --output_cum_allele_plot $outputcumalleleplot
