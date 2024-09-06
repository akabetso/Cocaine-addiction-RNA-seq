import pandas as pd
import io

#import the raw data
raw_metadata = pd.read_csv("../data/ref_genome/SraRunTable.txt", sep=',')
print(raw_metadata.head())
print(raw_metadata.columns)

#Extract the columns as needed and add replicates based on the biosamples
columns_to_keep = ["Run", "BioSample", "Sample Name", "sex", "tissue"]
shrunken_metadata = raw_metadata[columns_to_keep]
shrunken_metadata['replicates'] = shrunken_metadata.groupby('BioSample').cumcount() + 1
print(shrunken_metadata.head(10))

#Define the drug and the control
shrunken_metadata['Sample Name'] = shrunken_metadata['Sample Name'].str.replace(r'.*Saline.*', 'Saline', case=False, regex=True)
shrunken_metadata['Sample Name'] = shrunken_metadata['Sample Name'].str.replace(r'.*cocaine.*', 'Cocaine', case=False, regex=True)
print(shrunken_metadata)

#extract the runs runs\experiment to work with from the rest of the database
runs_to_keep = [
    "SRR15502495", "SRR15502506", "SRR15502562", "SRR15502573",
    "SRR15502584", "SRR15502595", "SRR15502832", "SRR15502843",
    "SRR15502854", "SRR15502865", "SRR15502965", "SRR15502976",
    "SRR15502988", "SRR15502999", "SRR15503143", "SRR15503154"
]
# Filter the DataFrame to keep only rows where 'Run' is in runs_to_keep
metadata = shrunken_metadata[shrunken_metadata['Run'].isin(runs_to_keep)]
metadata = metadata.drop(metadata.columns[1], axis=1)
metadata.to_csv("../results/featureCounts/metadata.csv", index=False)
metadata.rename(columns={'Sample Name': 'sample_name'}, inplace=True)
metadata

#load the feature counts data
counts_df = pd.read_csv("../results/featureCounts/cleaned_fcounts.txt", sep='\t')
counts_df.head()

#rename sample names to correspond with metadata
strip_extension = [name.replace('_aligned.bam', '') for name in counts_df.columns]
counts_df.columns = strip_extension
counts_df = counts_df.iloc[:, :-1]
counts_df.to_csv("../results/featureCounts/counts_matrix.csv", index=False)
counts_df
