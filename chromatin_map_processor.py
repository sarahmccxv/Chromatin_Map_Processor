from hicmatrix import HiCMatrix as hm
import pandas as pd
import math
import scipy.stats as stats
import matplotlib.pyplot as plt

def df_builder(file_path):
    """Returns a pandas dataframe for the bin x bin gene contact data stored in file_path
    >>> df_builder("500kbp_raw_human_trans_aggregate.h5")
    """
    hic_object = hm.hiCMatrix(file_path)
    bins = [(x[0], x[1], x[2]) for x in hic_object.cut_intervals]
    hic_matrix = hic_object.matrix.toarray()
    return pd.DataFrame(hic_matrix, index=bins, columns=bins)

def _gene_data_processor():
    """Returns a pandas dataframe containing the Gene stable ID, chromosome number, gene start site, and gene end site
    >>> _gene_data_processor()
    """
    dtype_mapping = {
        "Chromosome/scaffold name": str,
        "Gene start (bp)": int,
        "Gene end (bp)": int,
        "Gene stable ID": str
    }
    gene_df = pd.read_csv("biomart_GRCh38.p14.txt", sep='\t',
                          usecols=["Chromosome/scaffold name", "Gene start (bp)", "Gene end (bp)", "Gene stable ID"],
                          dtype=dtype_mapping)
    gene_df = gene_df[gene_df["Chromosome/scaffold name"].isin([str(i) for i in range(1,22)])]
    gene_df.columns = ["gene ID", "start", "end", "chromosome"]
    gene_df.set_index("gene ID", inplace=True)
    return gene_df.drop_duplicates()

def _get_bin(bin_x_bin_contacts, chr_number, base):
    """Returns the name of the bin that a base would fall in"""
    bin_size = bin_x_bin_contacts.index[0][2]
    start = math.floor(base / bin_size) * bin_size
    end = start + bin_size
    end_bins = pd.read_csv("HGA_CHR_LENGTHS_GRCh38.p14.csv").set_index("Chromosome").to_dict()['Total length (bp)']
    if end > end_bins[chr_number]:
        end = end_bins[chr_number]
    return ("chr" + chr_number, start, end)
    
def _calculate_contact_number(bin_x_bin_contacts, gene_data):
    """Returns the total number of trans contacts a gene has"""
    chr_number = gene_data[3]
    start_bin = _get_bin(bin_x_bin_contacts, chr_number, gene_data[1])
    end_bin = _get_bin(bin_x_bin_contacts, chr_number, gene_data[2])
    return bin_x_bin_contacts[start_bin: end_bin].sum().sum()

def get_gene_contacts_data(file_path, genes):
    """Returns a pandas dataframe containing processed data displaing the total number of trans contacts that a set of
    give genes have based on the HiC file given."""
    bin_x_bin_contacts = df_builder(file_path)
    genes = set(genes)
    
    gene_df = _gene_data_processor()
    gene_set = set(gene_df.index.tolist())
    found_genes = [gene for gene in genes if gene in gene_set]
    not_found_genes = genes - set(found_genes)

    result_data = []
    for gene_data in gene_df.loc[found_genes].itertuples():
        if gene_data.Index in genes:
            result_data.append({
                "gene": gene_data.Index,
                "number of trans contacts": _calculate_contact_number(bin_x_bin_contacts, gene_data)
            })

    result = pd.DataFrame(result_data)

    if not_found_genes:
        print("Data not available for ", len(not_found_genes), " genes")

    return result

def stat_significance(file_path, exp_genes, control_genes):
    """Returns the T-statistic and the p-value"""
    exp_contacts = get_gene_contacts_data(file_path, exp_genes)
    control_contacts = get_gene_contacts_data(file_path, control_genes)
    stat = stats.ttest_ind(a=exp_contacts["number of trans contacts"], b=control_contacts["number of trans contacts"])
    return (stat[0], stat[1])
    
def display_diagram(file_path, exp_genes, control_genes):
    exp_contacts = get_gene_contacts_data(file_path, exp_genes)
    control_contacts = get_gene_contacts_data(file_path, control_genes)
    p = stat_significance(file_path, exp_genes, control_genes)[1]
    
    data = pd.DataFrame({f'Experimental Group (N={len(exp_contacts)})': exp_contacts["number of trans contacts"],
                         f'Control Group (N={len(control_contacts)})': control_contacts["number of trans contacts"]})
    
    sns.boxplot(data=data)
    plt.ylabel('Total Number of Trans Contacts')
    plt.show()
    print(f"p = {p:.2e}")
