from Bio import SeqIO
from collections import defaultdict
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress, pearsonr
from kneed import KneeLocator
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import argparse

# argument para

par = argparse.ArgumentParser(description= "Codon Bias and Clsutring")
par.add_argument("--input", required= True, help= "Import your fasta file")
par.add_argument("--output", required= True, help= "path to output file")
arg = par.parse_args()

# creating output directory if not specified 
os.makedirs(arg.output, exist_ok=True)

# extracting codons
def codon_count(file):
    """
    //Counting the codons in given fasta file
    //giving path is neccasry
    //It returns a dictonary with codon counts
    """
    if not os.path.exists(file):
        raise FileNotFoundError

    codon_counts_dic = defaultdict(int) # dictonray to store codons 

    for record in SeqIO.parse(file, "fasta"): # parasing the fasta file by use of seqio library from biopython
        sequence = str(record.seq).upper() # convering to uppper case.
        if "N" in sequence:
            print(f"Skipping sequence {record.id} due to ambiguous bases.")
            continue
        if len(sequence) % 3 == 0: # checking if lenght is divisble ny 3 as only if this is true we can go ahead
            for i in range(0, len(sequence), 3): # taking 3 sequence at a tile
                codon = sequence[i:i+3] # extaring that by using slipicing
                codon_counts_dic[codon] += 1 # storing in our dictonary

    return codon_counts_dic

# file input
codon_counts = codon_count(arg.input) # applying function to our input file

#file = "C:\\Users\\vishw\\OneDrive\\Desktop\\New folder\\CBP\\GCF_000005845.2_ASM584v2_cds_from_genomic.fna"
# codon_counts = codon_count(file)
# print(codon_counts)



def frequency_of_codon(codon_counts):
    """ This function takes codon list and sums it calcualtes frequency of the codon
    """


    codon_count_num = codon_counts.values() # extracting only values from our list
    Total_codon_count =sum(codon_count_num) # suming the total values of codons

    codon_frequency = {} # empty dictonary to store codon frquency 

#  looping through our codn dictonary to calculate frequencey

    for codon, count in codon_counts.items():
        codon_frequency[codon] = count / Total_codon_count

    
    return codon_frequency

frequency = frequency_of_codon(codon_counts)
# mapping codons to amino acids 

codon_table = {
    "Phe": ["TTT", "TTC"],
    "Leu": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "Ile": ["ATT", "ATC", "ATA"],
    "Met": ["ATG"],  # Start codon
    "Val": ["GTT", "GTC", "GTA", "GTG"],
    "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "Pro": ["CCT", "CCC", "CCA", "CCG"],
    "Thr": ["ACT", "ACC", "ACA", "ACG"],
    "Ala": ["GCT", "GCC", "GCA", "GCG"],
    "Tyr": ["TAT", "TAC"],
    "Stop": ["TAA", "TAG", "TGA"],  # Stop codons
    "His": ["CAT", "CAC"],
    "Gln": ["CAA", "CAG"],
    "Asn": ["AAT", "AAC"],
    "Lys": ["AAA", "AAG"],
    "Asp": ["GAT", "GAC"],
    "Glu": ["GAA", "GAG"],
    "Cys": ["TGT", "TGC"],
    "Trp": ["TGG"],
    "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "Gly": ["GGT", "GGC", "GGA", "GGG"],
}
# Calculate total codon counts for each amino acd
aa_totals = {}
for amino_acid in codon_table: # using for loop to iterate over our codon table
    total_count = 0
    for codon in codon_table[amino_acid]:
        total_count += codon_counts.get(codon, 0)  # Use .get() to handle missing codons
    aa_totals[amino_acid] = total_count

# Display the results
# print("Total codon counts for each amino acid")
# for amino_acid, count in aa_totals.items():
#     print(f"{amino_acid} -  {count}")

rscu_values = {} # creating rscu dictoanry
for aa, codons in codon_table.items():
    # Total frequency for the amino acid
    total_aa_count = aa_totals[aa] 
    num_codons = len(codons) # Number of codons for the amino acid
    expected_freq = total_aa_count / num_codons     # Expected frequency for each codon
    
# Step 3: Calculate RSCU for each codon
    for codon in codons:
        observed_freq = codon_counts.get(codon, 0)  
        
        rscu = observed_freq / expected_freq
        # print(rscu)
        rscu_values[codon] = {
        "Amino Acid": aa,
        "Observed Count": observed_freq,
        "Expected Count": expected_freq,
        "RSCU": rscu
     }

rscu_df = pd.DataFrame.from_dict(rscu_values, orient="index")
rscu_df.reset_index(inplace=True)
rscu_df.rename(columns={"index": "Codon"}, inplace=True)

# Save to CSV
# rscu_csv_path = os.path.join(arg.output, "codon_values_RSCU.csv")
# rscu_df.to_csv("codon_values_RSCU.csv", index=False)
# print("\nRSCU values saved to 'codon_rscu_values.csv'")

# Sorting the data first by Amino Acid, then Codon
rscu_df_sorted = rscu_df.sort_values(by=["Amino Acid", "Codon"])

# Generate a bar chart for RSCU values
plt.figure(figsize=(16, 8))

# Sort by amino acids and codons for better visualization
rscu_df_sorted = rscu_df.sort_values(by=["Amino Acid", "Codon"])

#  the bar plot
sns.barplot(
    x="Codon",
    y="RSCU",
    hue="Amino Acid",
    data=rscu_df_sorted,
    dodge=False,  # Group by amino acid without separating bars
    palette="Set3"
)

# chart details
plt.title("RSCU Values by Codon and Amino Acid", fontsize=16)
plt.xlabel("Codon", fontsize=12)
plt.ylabel("RSCU Value", fontsize=12)
plt.xticks(rotation=45, fontsize=10)
plt.legend(title="Amino Acid", loc="upper right", fontsize=10)

# Saving plot to output directory
bar_chart_path = os.path.join(arg.output, "rscu_bar_chart.png") 
plt.tight_layout()
plt.savefig(bar_chart_path)
plt.close()

print(f"Bar chart saved to: {bar_chart_path}")

def GC_content_from_sequence(sequence):
    """
    Calculate GC content at the first (GC1), second (GC2), and third (GC3) codon positions,
    as well as GC12 (average of GC1 and GC2), directly from a raw DNA sequence.

    It gives : 
        gc1_per, gc2_per, gc3_per, gc12_per: GC content percentages.
    """
    sequence = sequence.upper()  # Ensure uppercase
    if "N" in sequence or len(sequence) % 3 != 0:
        return None  # Skip invalid sequences

    # Initialize counters
    gc1 = gc2 = gc3 = 0
    total_codons = len(sequence) // 3  # Total codons

    # Loop through the sequence codon by codon
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        if codon[0] in "GC":
            gc1 += 1
        if codon[1] in "GC":
            gc2 += 1
        if codon[2] in "GC":
            gc3 += 1

    # Calculate percentages
    gc1_per = (gc1 / total_codons) * 100
    gc2_per = (gc2 / total_codons) * 100
    gc3_per = (gc3 / total_codons) * 100
    gc12_per = (gc1_per + gc2_per) / 2

    return gc1_per, gc2_per, gc3_per, gc12_per
gc3_values = []
gc12_values = []

# Iterate over each gene in the FASTA file
for record in SeqIO.parse(arg.input, "fasta"):
    result = GC_content_from_sequence(str(record.seq))  # Pass sequence as a string
    if result:  # If the result is valid
        _, _, gc3, gc12 = result
        gc3_values.append(gc3)
        gc12_values.append(gc12)

gc_output_data = os.path.join(arg.output, "GC3_GC12_Values.txt")
with open(gc_output_data, "w") as file:
    file.write("GC3 Values:\n" + "\n".join(map(str, gc3_values)) + "\n\nGC12 Values:\n" + "\n".join(map(str, gc12_values))) #  writting the output in not txt file


# Initialize lists to store GC3 and GC12 values
gc3_values = []
gc12_values = []

# Flags to identify where GC12 values start
gc_section = None  # Will hold either GC3 or GC12 when processing

# Open the file and read line by line
with open(gc_output_data, 'r') as file:
    for line in file:
        line = line.strip()  # Remove any extra whitespace
        if line:  # Check if the line is not empty
            # If the line starts with GC3 or GC12, set the appropriate flag
            if line.startswith("GC3 Values:"):
                gc_section = "GC3"
                continue  # Skip this header line
            elif line.startswith("GC12 Values:"):
                gc_section = "GC12"
                continue  # Skip this header line

            # If the line is numeric (GC value)
            try:
                value = float(line)
                # Append the value to the correct list based on the current section
                if gc_section == "GC3":
                    gc3_values.append(value)
                elif gc_section == "GC12":
                    gc12_values.append(value)
            except ValueError:
                # Handle any lines that aren't valid numbers
                continue

# Ensure GC3 and GC12 lists have the same length
min_length = min(len(gc3_values), len(gc12_values))
gc3_values = gc3_values[:min_length]
gc12_values = gc12_values[:min_length]

# Display the number of values extracted for validation
print(f"Number of GC3 values - {len(gc3_values)}")
print(f"Number of GC12 values - {len(gc12_values)}")

rscu_data = rscu_df.drop(columns=["Codon", "Amino Acid"])  # Keep only numerical RSCU values

inertia = []
k_range = range(1, 10)  # Check K from 1 to 9
for k in k_range:
    kmeans = KMeans(n_clusters=k, random_state=42)
    kmeans.fit(rscu_data)
    inertia.append(kmeans.inertia_)

# Plot Elbow Curve
plt.figure(figsize=(8, 5))
plt.plot(k_range, inertia, marker='o')
plt.xlabel('Number of Clusters (K)')
plt.ylabel('Inertia')
plt.title('Elbow Method for Optimal K')
plt.grid(alpha=0.3)
elbow_curve = os.path.join(arg.output, "Elbow_curve.png")
plt.savefig(elbow_curve)
plt.close()


# Automatically select optimal K based on the elbow point
knee_locator = KneeLocator(range(1, 10), inertia, curve="convex", direction="decreasing")
optimal_k = knee_locator.knee
print(f"Optimal number of clusters (K) detected: {optimal_k}")


# Perform linear regression
slope, intercept, r_value, p_value, std_err = linregress(gc3_values, gc12_values)

# Create the regression line
regression_line = [slope * x + intercept for x in gc3_values]

correlation, corr_p_value = pearsonr(gc3_values, gc12_values)

# Plotting the data
plt.figure(figsize=(10, 6))

# Scatter plot for the actual data points
plt.scatter(gc3_values, gc12_values, alpha=0.5, s=10, label='Data Points')

# Regression line
plt.plot(gc3_values, regression_line, color='red', label=f'Regression Line (y = {slope:.2f}x + {intercept:.2f})')

# Adding labels, title, and grid
plt.title("Neutrality Graph with Regression Line", fontsize=14)
plt.xlabel("GC3 Content (%)", fontsize=12)
plt.ylabel("GC12 Content (%)", fontsize=12)
plt.grid(alpha=0.3)

# Adding statistical results on the graph
stats_text = f"R² = {r_value**2:.2f}\nCorrelation = {correlation:.2f}\np-value = {corr_p_value:.2e}"
plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes,
         fontsize=10, verticalalignment='top', bbox=dict(boxstyle="round", alpha=0.1))

plt.legend()

# saving 
neutrality_graph_path = os.path.join(arg.output, "neutrality_graph.png")
plt.savefig(neutrality_graph_path)
plt.close()
print("Analysis complete. Outputs saved to:", arg.output)


# Prepare RSCU data for clustering
rscu_data = rscu_df.drop(columns=["Codon", "Amino Acid"])  # Keep only numerical RSCU values

pca = PCA(n_components=2)  # Reduce data to 2 components
pca_result = pca.fit_transform(rscu_data)  # Transform RSCU data
rscu_df["PC1"] = pca_result[:, 0]  #  PC1 to the dataframe
rscu_df["PC2"] = pca_result[:, 1]  # PC2 to the dataframe


# apply K-means clustering with optimal K
kmeans = KMeans(n_clusters=optimal_k, random_state=42)
rscu_df["Cluster"] = kmeans.fit_predict(rscu_data)
# saving the cluster map
plt.figure(figsize=(10, 7))
sns.scatterplot(x="PC1", y="PC2", hue="Cluster", data=rscu_df, palette="viridis", s=100, alpha=0.7)
plt.title("Clustering of Genes Based on RSCU Values")
plt.xlabel("PCA 1")
plt.ylabel("PCA 2")
plt.legend(title="Cluster")
plt.grid(alpha=0.3)
clustering_path = os.path.join(arg.output, "Clustering_graph.png")
plt.savefig(clustering_path)
plt.close()
print("Analysis complete. Outputs saved to:", arg.output)

print("Completed all steps")