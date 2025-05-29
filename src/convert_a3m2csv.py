import argparse
import subprocess
import json
import pandas as pd
from pathlib import Path
from tqdm import tqdm

def fasta_to_list(fasta_string):
    keys = []
    values = []
    current_key = None
    current_seq = []

    # Handle empty input
    if not fasta_string.strip():
        return keys, values

    # Process each line
    for line in fasta_string.strip().split('\n'):
        line = line.strip()

        # Skip empty lines
        if not line:
            continue

        # Header line
        if line.startswith('>'):
            # Save the previous sequence if it exists
            if current_key:
                keys.append(current_key)
                values.append(''.join(current_seq))

            # Start new sequence
            current_key = line[1:]  # Remove the '>' character
            current_seq = []

        # Sequence line
        else:
            current_seq.append(line)

    # Save the last sequence
    if current_key:
        keys.append(current_key)
        values.append(''.join(current_seq))

    return keys, values

def convert_to_csv(inputs, msa_dir, prefix = "sequence"):
    # alphafold3の入力ファイルからboltz-1のMSA入力のためのCSVファイルの作成

    # Process each sequence in the inputs
    for i, sequence in enumerate(inputs["sequences"]):
        if "protein" not in sequence:
            continue

        msa_list = []

        # Process paired MSA
        if "pairedMsa" in sequence["protein"]:
            keys, values = fasta_to_list(sequence["protein"]["pairedMsa"])
            # paired配列の処理は、DUMMY配列以外、indexを揃える
            for j in range(len(values)):
                if keys[j] != "DUMMY":
                    msa_list.append({
                        "key": j,
                        "sequence": values[j]
                    })

        # Process unpaired MSA
        if "unpairedMsa" in sequence["protein"]:
            keys, values = fasta_to_list(sequence["protein"]["unpairedMsa"])
            # unpaired配列の処理は、indexを-1とする
            for j in range(len(values)):
                if keys[j] != "DUMMY":
                    msa_list.append({
                        "key": -1,
                        "sequence": values[j]
                    })

        if msa_list:
            print(f"Converting MSA to CSV: {prefix}_{i}.csv")
            df = pd.DataFrame(msa_list)
            df.to_csv(msa_dir/f"{prefix}_{i}.csv", index=False)

def convert_msa_to_json(input_msa_file, json_file):
    # msatojsonを用いてMSA情報を変換
    ret = subprocess.run(["msatojson", "-i", input_msa_file, "-o", json_file])
    if ret.returncode != 0:
        raise Exception(f"Error: {ret.stderr}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_msa_file", type=str)
    parser.add_argument("working_dir", type=str)
    args = parser.parse_args()

    input_msa_file = Path(args.input_msa_file)
    if not input_msa_file.exists():
        raise FileNotFoundError(f"Input MSA file {input_msa_file} does not exist")
    working_dir = Path(args.working_dir)
    working_dir.mkdir(parents=True, exist_ok=True)
    json_file = working_dir/(input_msa_file.stem + ".json")

    print(f"Converting MSA file to JSON: {input_msa_file}")
    convert_msa_to_json(input_msa_file, json_file)
    with open(json_file, "r") as f:
        inputs = json.load(f)
    print(f"Converting MSA to CSV: {json_file}")
    convert_to_csv(inputs, working_dir, prefix=input_msa_file.stem)
    print("MSA conversion completed successfully")
