import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pathlib import Path
from multiprocessing import Pool
from functools import partial
from parse_vcf_files import parse_asm_vcf_tandemtwister




directory = "/confidential/tGenVar/vntr/output_maryam/tools/run_all_tools/output/hgsvc_test/TandemTwist/asm/"
sample_info_file = "/confidential/tGenVar/vntr/output_maryam/tools/run_all_tools/input_folder/samples/hgsvc/HGSVC2_sample_info.tsv"
output_file = 'output/merged_pop_df.csv'

def get_sample_info_df(sample_info_file):
    '''
    reads sample_info file --> columns: sample, sex, pop, superpop
    This function duplicates the df with sample names with two added suffices: _h1 and _h2
    '''

    print('reading sample info')

    sample_info = pd.read_csv(sample_info_file, sep='\t')
    h2 = sample_info.copy()
    h2['sample'] = h2['sample'] + '_h2'

    sample_info['sample'] = sample_info['sample'] + '_h1'

    return pd.concat([sample_info, h2], ignore_index=True)

def sample_to_superpop(sample_info_file):
    '''
    reads sample_info file --> columns: sample, sex, pop, superpop
    returns a dictionary: sample -> superpop
    '''

    print('reading sample info')

    sample_info = pd.read_csv(sample_info_file, sep='\t')
    sample_to_superpop = dict(zip(sample_info['sample'], sample_info['superpop']))

    return sample_to_superpop


# Function to process a single file
def process_vcf(file_path):                        #, sample_info
    sample_name = os.path.basename(file_path)[:-7] # Remove '.vcf.gz'
    print(f"reading vcf file for {sample_name}")
    df = parse_asm_vcf_tandemtwister(pysam.VariantFile(file_path))
    df["sample"] = sample_name

    return df
    # return df.merge(sample_info, on='sample', how='inner')



def compute_population_df(directory, output):
    '''
    input: a dictionary including these vcf files:
    {sample}_h1.vcf.gz and {sample}_h2_vcf.gz for all samples in the cohort

    sample info:
    A file containing the sample information in the population
    
    output:
        cols: chrom,start,end,motifs,motif_ids,CN_ref,sample1_h1, sample1_h2, sample2_h1, sample2_h2, ... 
        sample_h possible values: CN,GT,len_,seq=motif_concat)
        one separate output file per measured value (CN, len_, seq)
    '''

    vcf_files = list(Path(directory).glob("*.vcf.gz"))

    df_pop = pd.DataFrame()

    for f in vcf_files:

        if not f.is_file():
            raise FileNotFoundError(f"VCF file not found: {f}")
        

        sample = os.path.basename(f)[:-7]  # Remove '.vcf.gz'
        print(f"Processing {sample}")

        df = process_vcf(f)
        df = df[["chrom", "start", "end", "motifs", "CN_ref", "CN", "len_", "motif_concat"]]

        # rename the last three columns to len and seq; add sample suffix
        df.rename(columns={"CN": f"CN_{sample}"}, inplace=True)
        df.rename(columns={"len_": f"len_{sample}"}, inplace=True)
        df.rename(columns={"motif_concat": f"seq_{sample}"}, inplace=True)

        if df_pop.empty:
            df_pop = df
        else:
            df_pop = pd.merge(df_pop, df, on=["chrom", "start", "end", "motifs", "CN_ref"], how="outer", suffixes=('', ''))

    df_pop.to_csv(output, index=False, sep='\t')

    return df_pop


def count_alleles(df, column="CN", sample_size=70, output="output/population/allele_counts.csv"):
    '''
    df: population dataframe
    column: column to count alleles on (possible values: CN, len_, motifs)
    sample_size: sample size

    output:
        allele_count_summary: summary of the number of alleles in the population
    '''

    pass

    return allele_count_summary


def plot_allele_count(cn_allele_counts, len_allele_counts, motif_allele_counts, n=5, output="output/population/multi_facet_allele_barplots.pdf"):
    '''
    makes a bar plot of the number of regions with 1,2, ..., n, >n allele counts
    '''

    def group_allele_counts(count_dict, n=4):
        series = pd.Series(count_dict)
        grouped = series[series.index <= n]
        gtn_sum = series[series.index > n].sum()
        plot_series = pd.concat([grouped, pd.Series({f'>{n}': gtn_sum})])
        plot_series.index = plot_series.index.astype(str)
        return plot_series

    # Prepare data
    vectors = {
        "CN alleles": group_allele_counts(cn_allele_counts, n),
        "length alleles": group_allele_counts(len_allele_counts, n),
        "motif alleles": group_allele_counts(motif_allele_counts, n),
    }

    # Setup plot
    sns.set(style="whitegrid", font_scale=1.2)
    fig, axes = plt.subplots(1, 3, figsize=(7, 5), sharey=True)
    fig.supxlabel("Number of Alleles", fontsize=12)

    colors = sns.color_palette("Set2")

    # Plot each vector
    for ax, (title, data), color in zip(axes, vectors.items(), colors):
        sns.barplot(x=data.index, y=data.values, ax=ax, color=color)
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("")
        ax.set_ylabel("Number of TR regions" if ax == axes[0] else "")  # only first has ylabel
        ax.tick_params(axis='x', rotation=0)

    # Overall adjustments
    plt.suptitle("Allele Count Distribution Across TR Regions", fontsize=12)
    # plt.xlabel("Number of alleles")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output, format='pdf', bbox_inches='tight')
    # plt.show()

    #return plt

def subset_columns(df, columns=["chrom", "start", "end", "motifs", "CN_ref"], prefix_col="CN_"):
    """
    Subsets the dataframe to only include specified columns and those starting with a given prefix.

    Parameters:
        df (pd.DataFrame): Input dataframe
        columns (list): List of specific columns to retain
        prefix_col (str): Prefix to identify additional columns to retain. Possible values: CN_, len_, seq_

    Returns:
        pd.DataFrame: Subsetted dataframe
    """

    cols_to_keep = [col for col in df.columns if col in columns or col.startswith(prefix_col)]

    return df[cols_to_keep]


def compute_cumulative_allele_classification(df, sample_info, col_prefix="CN", sample_order=None, superpop_order=['AFR', 'AMR', 'EUR', 'EAS', 'SAS'], 
                                             output='output/population/cumulative_allele.csv'): # value_col="CN"
    """
    Creates a stacked bar plot showing cumulative allele diversity (singleton, biallelic, multiallelic) across samples.

    Parameters:
        df: population df, cols: chrom,start,end,motifs,CN_ref,CN_{sample},len_{sample},seq_{sample} for all {sample} in samples
        prefix_col (str): Column to count alleles from ('CN_', 'len_', or 'seq_')
    
    output:
        summary_df (pd.DataFrame): DataFrame with columns ['sample', 'singleton', 'biallelic', 'multiallelic', 'superpop']
    """

    # TODO: only for CN so far (for now), we need to extend it to len_ and motifs
    # TODO: This part is very slow, especillay reshaping part... try to improve... think about multi-sample vcf file

    #all_samples = df['sample'].drop_duplicates().tolist()

    # get all columns starting with 'CN_'
    all_samples = [col.replace('CN_', '') for col in df.columns if col.startswith('CN_') and col != 'CN_ref']
    print('all samples:', all_samples)

    subsample_info = sample_info[sample_info['sample'].isin(all_samples)]

    sample_superpop = subsample_info[['sample', 'superpop']]

    if superpop_order:
        sample_superpop['superpop'] = pd.Categorical(sample_superpop['superpop'], categories=superpop_order, ordered=True)
    sample_superpop = sample_superpop.sort_values(by=['superpop', 'sample'])
    if sample_order is None:
        sample_order = sample_superpop['sample'].tolist()

    
    # create two dataframes for computing allele sets and counting the number of singleton, biallelic, and multiallelic
    df_alleles_set = df[["chrom", "start", "end"]]
    df_num_alleles = df[["chrom", "start", "end"]]

    # computing and counting the set of alleles seen up to each added sample
    prev_sample = None
    for sample in sample_order:
        print(f'processing sample {sample}')

        if df_alleles_set.shape[1] == 3:
            # no sample is added yet
            df_alleles_set[[f"alleles_upto_{sample}"]] = df[f'{col_prefix}_{sample}'].apply(lambda x: set([x]) if pd.notna(x) else set())
            df_num_alleles[f"num_alleles_upto_{sample}"] = df_alleles_set[f"alleles_upto_{sample}"].apply(len)
        else:
            df_alleles_set[f"alleles_upto_{sample}"] = df_alleles_set[f"alleles_upto_{prev_sample}"].combine(df[f'{col_prefix}_{sample}'], 
                                                                            lambda s1, s2: s1.union(set([s2])) if pd.notna(s2) else s1)
            df_num_alleles[f"num_alleles_upto_{sample}"] = df_alleles_set[f"alleles_upto_{sample}"].apply(len)

        prev_sample = sample


    # Select sample columns
    sample_cols = df_num_alleles.columns[3:]  # columns from sample1 onwards
    df_samples = df_num_alleles[sample_cols]

    # Vectorized counts
    singleton = (df_samples == 1).sum()
    biallelic = (df_samples == 2).sum()
    multiallelic = (df_samples > 2).sum()

    # Assemble summary dataframe
    summary_df = pd.DataFrame({
        "sample": sample_cols,
        "singleton": singleton.values,
        "biallelic": biallelic.values,
        "multiallelic": multiallelic.values
    })

    summary_df = pd.merge(summary_df, subsample_info, on='sample')

    summary_df.to_csv(output, index=False, sep='\t')

    return summary_df


def plot_cumulative_allele_classification(df, sample_order=None, superpop_order=None, palette=None, figsize=(15, 6), output='output/population/commulative_alleles.pdf'):
    """
    Creates a stacked bar plot showing cumulative allele diversity (sample, singleton, biallelic, multiallelic) across samples.

    Parameters:
        df (pd.DataFrame): DataFrame with columns ['sample', 'singleton', 'biallelic', 'multiallelic', 'superpop']
        value_col (str): Column to count alleles from ('CN' or 'len_')
        sample_order (list): List of sample names in desired plotting order (optional)
        superpop_order (list): Order of superpopulations for sorting (optional)
        palette (dict): Dictionary mapping class ('Singleton', etc) to color
        figsize (tuple): Figure size
        save_path (str): If provided, saves the figure as PDF
    """
    
    # Custom superpop color map
    superpop_colors = {
        'AFR': 'gold',
        'AMR': 'red',
        'EAS': 'green',
        'EUR': 'blue',
        'SAS': 'purple'
    }

    # Allele type bar colors (clearly distinct)
    allele_colors = {
        'singleton': 'gold', # '#e41a1c',    # red
        'biallelic': 'skyblue', # '#377eb8',    # blue
        'multiallelic': 'tomato' # '#4daf4a'  # green
    }

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Stacked bar chart
    bottom = [0] * len(df)
    for col in ['singleton', 'biallelic', 'multiallelic']:
        ax.bar(df['sample'], df[col], bottom=bottom, label=col, color=allele_colors[col], alpha=1)
        bottom = [i + j for i, j in zip(bottom, df[col])]

    # Add superpop color bands below bars
    y_offset = -max(bottom) * 0.1
    bar_height = max(bottom) * 0.05

    # Draw one color band per superpop across its samples
    unique_superpops = df['superpop'].unique()
    for superpop in unique_superpops:
        indices = df.index[df['superpop'] == superpop].tolist()
        start = indices[0] - 0.4
        end = indices[-1] + 0.4
        width = end - start
        ax.add_patch(plt.Rectangle((start, y_offset), width, bar_height,
                                color=superpop_colors.get(superpop, 'gray'),
                                transform=ax.transData, clip_on=False))

    # Superpop legend
    handles_superpop = [
        plt.Rectangle((0, 0), 1, 1, color=superpop_colors[sp])
        for sp in df['superpop'].unique()
    ]
    labels_superpop = list(df['superpop'].unique())
    legend_superpop = ax.legend(handles_superpop, labels_superpop, title='Superpopulation',
                                loc='upper left', bbox_to_anchor=(1.02, 1))

    # Allele type legend (moved outside)
    handles_alleles = [
        plt.Rectangle((0, 0), 1, 1, color=allele_colors[atype])
        for atype in allele_colors
    ]
    labels_alleles = list(allele_colors.keys())
    legend_alleles = ax.legend(handles_alleles, labels_alleles, title='Allele Type',
                            loc='upper left', bbox_to_anchor=(1.02, 0.6))

    ax.add_artist(legend_superpop)

    # Plot styling
    ax.set_ylabel("Count")
    ax.set_xlabel("Sample")
    ax.set_title("Allele Types per Sample (Stacked)")
    plt.xticks(rotation=45)

    # Space for annotations
    plt.subplots_adjust(bottom=0.2, right=0.75)

    ax.set_xticklabels(ax.get_xticklabels(), fontsize=8)
    plt.savefig(output, format='pdf')
    
    return plt



def compute_pca(df, output=None):
    '''
    df = df_all, concat of all sample dfs, generated by compute_population_df
    '''

    # # Step 1: Filtering
    # # filtering ----
    df = df[df['chrom'].isin([f'chr{i}' for i in range(1, 23)])]

    # # Step 2: Create a unique region ID to pivot
    df['region'] = df['chrom'].astype(str) + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
    # 126 s

    # # Step 3: Pivot to sample x region matrix (samples as rows)
    pivot_df = df.pivot_table(index='sample', columns='region', values='CN', aggfunc='mean')
    # 67 s

    # # Step 4: Fill missing CNs (regions not observed for a sample)
    pivot_df = pivot_df.fillna(0)

    std_threshold = 0
    region_std = pivot_df.std(axis=0)
    high_var_regions = region_std[region_std > std_threshold].index
    pivot_df = pivot_df[high_var_regions]
    # 1.9 s


    # # 547949 regions with more than 0 std


    # # Step 5: Normalize (z-score scaling)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(pivot_df)

    # # Step 6: PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    ## 8.4 s

    # # Step 7: Prepare for plotting
    pca_df = pd.DataFrame(X_pca, columns=['PC1', 'PC2'])
    pca_df['sample'] = pivot_df.index

    # # Map sample → superpop
    superpop_map = sample_info[['sample', 'superpop']].drop_duplicates().set_index('sample')
    pca_df['superpop'] = pca_df['sample'].map(superpop_map['superpop'])

    # Step 8: Plot
    # remove the outlier sample
    pca_df_filt = pca_df[pca_df['sample']!='HG00732_h1']
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=pca_df_filt, x='PC1', y='PC2', hue='superpop', palette='tab10')
    plt.title('PCA of CN Profiles by Sample')
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.legend(title='Superpopulation')
    plt.tight_layout()
    # plt.show()
    
    return pca_df, plt




