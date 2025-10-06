import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import LabelEncoder, StandardScaler
from pathlib import Path
from multiprocessing import Pool
from functools import partial
from parse_vcf_files import parse_asm_vcf_tandemtwister
from matplotlib.ticker import ScalarFormatter



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

    for i, f in enumerate(vcf_files):

        if not f.is_file():
            raise FileNotFoundError(f"VCF file not found: {f}")
        

        sample = os.path.basename(f)[:-7]  # Remove '.vcf.gz'
        print(f"Processing sample {i}: {sample}")

        #df = process_vcf(f)

        df = parse_asm_vcf_tandemtwister(pysam.VariantFile(f))
        df = df[["chrom", "start", "end", "motifs", "CN_ref", "CN", "len_", "motif_concat"]]

        # rename the last three columns to len and seq; add sample suffix
        df.rename(columns={"CN": f"CN_{sample}"}, inplace=True)
        df.rename(columns={"len_": f"len_{sample}"}, inplace=True)
        df.rename(columns={"motif_concat": f"seq_{sample}"}, inplace=True)

        if df_pop.empty:
            df_pop = df
        else:
            df_pop = pd.merge(df_pop, df, on=["chrom", "start", "end", "motifs", "CN_ref"], how="inner", suffixes=('', ''))

    df_pop.to_csv(output, index=False, sep='\t')

    return df_pop



def compute_cumulative_allele_classification(df, sample_info, col_prefix="CN", sample_order=None, 
                                             superpop_order=['AFR', 'AMR', 'EUR', 'EAS', 'SAS'], 
                                             outputdir='output/population/cumm_allele'):
    """
    Computes cumulative allele diversity (singleton, biallelic, multiallelic) across samples.
    It shows how this diversity changes as we add samples one by one.

    Parameters:
        df: population df, cols: chrom,start,end,motifs,CN_ref,CN_{sample},len_{sample},seq_{sample} for all {sample} in samples
        sample_info: sample information dataframe
        col_prefix (str): Column to count alleles from ('CN_', 'len_', or 'seq_')
        sample_order (list): List of sample names in desired order (optional)
        superpop_order (list): Order of superpopulations for sorting (optional)
        output (str): Path to save the cumulative allele classification CSV

    Returns:
        summary_df (pd.DataFrame): DataFrame with columns ['sample', 'singleton', 'biallelic', 'multiallelic', 'superpop']
    """

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
    for i, sample in enumerate(sample_order):
        print(f'processing sample {i}: {sample}')

        if df_alleles_set.shape[1] == 3:
            # no sample is added yet
            df_alleles_set.loc[:, f"{sample}"] = df[f'{col_prefix}_{sample}'].apply(lambda x: set([x]) if pd.notna(x) else set())
        
        else:
            df_alleles_set.loc[:, f"{sample}"] = df_alleles_set[f"{prev_sample}"].combine(df[f'{col_prefix}_{sample}'], 
                                                                            lambda s1, s2: s1.union(set([s2])) if pd.notna(s2) else s1)
        
        df_num_alleles.loc[:, f"{sample}"] = df_alleles_set[f"{sample}"].apply(len)

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

    summary_df.to_csv(os.path.join(outputdir, f'cumm_{col_prefix}_allele_class_summary.csv'), index=False, sep='\t')
    df_alleles_set.to_csv(os.path.join(outputdir, f'cumm_{col_prefix}_allele_set.csv'), index=False, sep='\t')
    df_num_alleles.to_csv(os.path.join(outputdir, f'cumm_{col_prefix}_num_alleles.csv'), index=False, sep='\t')

    return summary_df


def plot_cumulative_allele_classification(df, sample_order=None, superpop_order=None, value_col='CN',
                                          palette=None, figsize=(15, 6), outputdir='output/population/cumm_allele'):
    """
    Creates a stacked bar plot showing cumulative allele diversity (sample, singleton, biallelic, multiallelic) across samples.

    Parameters:
        df (pd.DataFrame): summary DataFrame (returned by compute_cumulative_allele_classification) 
            with columns ['sample', 'singleton', 'biallelic', 'multiallelic', 'superpop']
        value_col (str): Column used to count alleles from ('CN', 'len', or 'seq')
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
        'singleton': 'teal', # '#e41a1c',    # red
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
                                loc='upper left', bbox_to_anchor=(1.1, 1))

    # Allele type legend (moved outside)
    handles_alleles = [
        plt.Rectangle((0, 0), 1, 1, color=allele_colors[atype])
        for atype in allele_colors
    ]
    labels_alleles = list(allele_colors.keys())
    legend_alleles = ax.legend(handles_alleles, labels_alleles, title='Allele Type',
                            loc='upper left', bbox_to_anchor=(1.1, 0.5))

    ax.add_artist(legend_superpop)

    # Plot styling
    ax.set_ylabel("Counts")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0)) # use scientific notation

    # Add secondary y-axis (percentages)
    ax2 = ax.twinx()

    total = sum(df.loc[0, ['singleton', 'biallelic', 'multiallelic']])

    ax2.set_ylim(0, ax.get_ylim()[1] / total * 100)  # scale to percentages
    ax2.set_ylabel("Percentage (%)")

    # Optionally, set ticks to 0-20-40-...-100
    ax2.set_yticks(np.linspace(0, 100, 6))

    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df['sample'], rotation=90, ha="center", fontsize=6)

    ax.set_title(f"TR {value_col} alleles polymorphism in population (per added sample)")

    ax.margins(x=0.01)
    plt.tight_layout()

    # Space for annotations
    plt.subplots_adjust(bottom=0.2, right=0.75)

    # plt.savefig(os.path.join(outputdir, f'cumm_alleles_{value_col}_barplot.pdf'), format='pdf')
    
    return plt

def plot_pie_chart_allele_counts(cumm_df_CN, cumm_df_len, cumm_df_seq):
    """
    Plots pie charts showing the overall distribution of singleton, biallelic, and multiallelic loci
    for CN, length, and sequence alleles.

    Parameters:
        cumm_df_CN (pd.DataFrame): Cumulative allele classification DataFrame for CN alleles, computed by compute_cumulative_allele_classification
        cumm_df_len (pd.DataFrame): Cumulative allele classification DataFrame for length alleles, computed by compute_cumulative_allele_classification
        cumm_df_seq (pd.DataFrame): Cumulative allele classification DataFrame for sequence alleles, computed by compute_cumulative_allele_classification
    Returns:
        plt: Matplotlib pyplot object containing the pie charts
    """
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Define colors to match your allele categories
    colors = ['teal', 'skyblue', 'tomato']
    labels = ['Singleton', 'Biallelic', 'Multiallelic']

    # Pie chart for CN alleles
    axes[0].pie(
        cumm_df_CN[['singleton', 'biallelic', 'multiallelic']].tail(1).values.flatten(),
        # labels=labels,
        autopct='%1.1f%%',
        colors=colors
    )
    axes[0].set_title('CN Alleles')

    # Pie chart for length alleles
    axes[1].pie(
        cumm_df_len[['singleton', 'biallelic', 'multiallelic']].tail(1).values.flatten(),
        # labels=labels,
        autopct='%1.1f%%',
        colors=colors
    )
    axes[1].set_title('Length Alleles')

    # Pie chart for sequence alleles
    axes[2].pie(
        cumm_df_seq[['singleton', 'biallelic', 'multiallelic']].tail(1).values.flatten(),
        # labels=labels,
        autopct='%1.1f%%',
        colors=colors
    )
    axes[2].set_title('Sequence Alleles')

    # Add a shared legend for allele types
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colors]
    fig.legend(handles, labels, title='Allele Type', loc='upper right')

    plt.tight_layout()

    return fig

def compute_pca(df, sample_info, value_col='CN', std_threshold=0, output=None):
    '''
    df = df_pop, merged sample dfs, generated by compute_population_df,
    col: 'CN', 'len', 'seq': the column prefix to use for PCA
    output: if provided, saves the PCA plot as a PDF
    '''

    # # Step 1: Filtering -- subset rows and columns; only autosomes; only the columns with the value_col prefix for computing PCA

    df = df[df['chrom'].isin([f'chr{i}' for i in range(1, 23)])]

    sample_cols = [col for col in df.columns if col.startswith(f'{value_col}_') and col != 'CN_ref']
    
    # Step 2: reshape to matrix: samples x loci
    df_mat = df[sample_cols].T
    df_mat.index = df_mat.index.str.replace(f'{value_col}_', '', regex=False)

    # Step 2a: Encode sequence values if value_col == 'seq'
    if value_col == 'seq':
        # Each column is a locus; label-encode each locus separately
        for i, col in enumerate(df_mat.columns):
            le = LabelEncoder()
            # Fill NaNs with a placeholder string
            df_mat[col] = le.fit_transform(df_mat[col].fillna('NA'))
    else:
        # numeric values (CN, len)
        df_mat.fillna(0, inplace=True)

    # # 547949 regions with more than 0 std
    # 404030 regions with more than 0 std --> much more PCA variance explained

    scaler = StandardScaler(with_mean=True, with_std=True)
    X_scaled = scaler.fit_transform(df_mat)

    # # Step 4: PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    pca_df = pd.DataFrame(X_pca, columns=['PC1', 'PC2'])
    pca_df['sample'] = df_mat.index

    # 12 s

    # # Map sample → superpop
    pca_df = pd.merge(pca_df, sample_info[['sample', 'superpop']], on='sample', how='left')


    # Step 5: Plotting
    # remove the outlier sample
    # pca_df_filt = pca_df[pca_df['sample']!='HG00732_h1']

    plt.figure(figsize=(10, 6))

    # Custom superpop color map
    superpop_colors = {
        'AFR': 'gold',
        'AMR': 'red',
        'EAS': 'green',
        'EUR': 'blue',
        'SAS': 'purple'
    }
    
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='superpop', palette=superpop_colors)
    plt.title(f'PCA of {value_col} Profiles by Sample')
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.legend(title='Superpopulation')
    plt.tight_layout()
    # plt.show()
    
    return pca_df, plt




