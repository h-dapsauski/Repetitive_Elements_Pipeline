import os
import pandas as pd 

#set the folder where all quast results are
#change it based on the quast output directory
base_folder = 'Quast_Output'

#metrics to extract from each report.tsv
metrics = {
    'N50':'N50',
    'Genome fraction (%)':'Genome fraction (%)',
    '# misassemblies':'misassemblies',
    '# mismatches per 100 kbp': '# mismatches per 100 kbp'
}

#empty list to collect the result
all_samples =[]

#each folder in quast output directory
for folder in os.listdir(base_folder):
    folder_path = os.path.join(base_folder, folder)
    report_file= os.path.join(folder_path, "report.tsv")
    #skipping if not exist
    if not os.path.exists(report_file):
        continue

    #extract sample ID and depth
    parts = folder.split('_')
    if len(parts) < 4:
        continue
    sample_id = '_'.join(parts[:4])
    coverage = parts[4] + 'x'

    #start a new row
    row = {"Sample ID": sample_id, "Coverage": coverage}

    #open and read the report.tsv file line by line
    with open(report_file, 'r') as file:
        for line in file:
            if '\t' not in line:
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            metric_name = parts[0].strip()
            spades_value = parts[1].strip()
            unicycler_value = parts[2].strip()

            #choose if it is the metric that we want
            if metric_name in metrics:
                row[metrics[metric_name] + "SPAdes"] = float(spades_value)
                row[metrics[metric_name] + "Unicycler"] = float(unicycler_value)

    all_samples.append(row)

#save the result
df = pd.DataFrame(all_samples)
df.to_csv("quast_summary.csv", index=False)
print("Done")


