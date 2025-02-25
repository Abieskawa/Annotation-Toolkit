import pandas as pd
import re
import argparse

def extract_query_number(query):
    match = re.search(r'(\d+)', query)
    return int(match.group(1)) if match else float('inf')

def simple_interproscan_tsv(interproscan_df, prefix):
    output_file = f"{prefix}_IPR-GO_ann.tsv"

    df_filtered = interproscan_df[["Query", "Source","Start", "End", "IPR_ID","IPR_Description", "GO_Terms"]].copy()

    df_filtered.loc[:, "GO_Terms"] = df_filtered["GO_Terms"].replace("-", "No GO terms")

    df_filtered["Position"] = df_filtered.apply(lambda row: "-" if row["IPR_ID"] == "-" else f'{row["Start"]}-{row["End"]}', axis=1)

    def format_go_terms(go_terms):
        if pd.isna(go_terms) or go_terms == "No GO terms":
            return "No GO terms"
        go_ids = re.findall(r"GO:\d+", go_terms)  
        unique_go_ids = sorted(set(go_ids))  
        return ",".join(unique_go_ids) if unique_go_ids else "No GO terms"

    df_filtered.loc[:, "GO_Terms"] = df_filtered["GO_Terms"].apply(format_go_terms)

    df_filtered =  df_filtered[~(( df_filtered["IPR_ID"] == "-") & ( df_filtered["GO_Terms"] == "No GO terms"))]

    df_group1 = df_filtered.sort_values("Start").groupby(["Query", "Source", "IPR_ID", "IPR_Description"], as_index=False).agg({
        "GO_Terms": lambda x: ",".join(sorted(set(filter(lambda go: go != 'No GO terms', ",".join(x).split(","))))) or 'No GO terms',
        "Position": lambda x: ",".join(x) if "-" not in x.values else "-"
    })

    df_group1.loc[:, "Source_Position"] = df_group1.apply(lambda row: f'{row["Source"]} ({row["Position"]})', axis=1)

    df_group2 = df_group1.groupby(["Query", "IPR_ID", "IPR_Description"], as_index=False).agg({
        "GO_Terms": lambda x: ",".join(sorted(set(filter(lambda go: go != 'No GO terms', ",".join(x).split(","))))) or 'No GO terms',
        "Source_Position": lambda x: "|".join(x) if "-" not in x.values else "-"
    })
    
    df_group2["Query_sort"] = df_group2["Query"].apply(extract_query_number)

    output_df = df_group2.loc[:, ["Query", "IPR_ID", "IPR_Description", "Source_Position", "GO_Terms", "Query_sort"]].sort_values(["Query_sort", "IPR_ID"]).drop(columns=["Query_sort"])


    output_df.to_csv(output_file, sep="\t", index=False)
    

def main():
    parser = argparse.ArgumentParser(description="Extract IPR and GO annotations from InterProScan TSV file.")
    parser.add_argument("-in", "--input_file", required=True, help="Path to the input InterProScan TSV file.")
    parser.add_argument("-p", "--prefix", required=True, help="Prefix name for the output TSV file.")
    args = parser.parse_args()

    input_file = args.input_file
    prefix = args.prefix

    columns = [
        "Query", "MD5", "SeqLength", "Source", "Signature", "Description", "Start", "End", "E-value", "Status",
        "Date", "IPR_ID", "IPR_Description", "GO_Terms", "Misc"
    ]
    interproscan_df = pd.read_table(input_file, names=columns)

    simple_interproscan_tsv(interproscan_df, prefix)

if __name__ == "__main__":
    main()