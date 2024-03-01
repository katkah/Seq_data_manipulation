import random
import subprocess
import sys
import argparse



def read_fastq(file_path):
    """
    Read a FASTQ file and return a list of names
    """
    sequences= []
    with open(file_path, 'r') as file:
        while True:
            seq_id = file.readline().strip()
            if not seq_id:
                break
            _ = file.readline() #skip seq
            _ = file.readline()  # skip the '+' separator line
            _ = file.readline() #skip qual_scores
            sequences.append(seq_id)
    return sequences

def split_fastq(sequences, number, file_name):
    """
    Randomly split the FASTQ file into two separate files.
    """

    random.shuffle(sequences)
    split_index = len(sequences) // number
    sequences1 = sequences[:split_index]
    sequences2 = sequences[split_index:]
    
    output_file1 = file_name + ".names1.txt"
    output_file2 = file_name + ".names2.txt"

    write_names(sequences1, output_file1)
    write_names(sequences2, output_file2)

def write_names(sequences, output_file):
    """
    Write sequence names
    """
    with open(output_file, 'w') as file:
        for i, seq in enumerate(sequences):
            if i == len(sequences) - 1:  # If it's the last element
                file.write(seq)
            else:
                file.write(seq + '\n')


def call_seqtk(file_name, number):
    for i in range(1,number + 1):
        output_file = file_name + ".names" + str(i) + ".txt"  
        parts = file_name.split('.')
        base_name = parts[0]
        final = base_name + "_part" + str(i) + ".fastq" 
        cmd = f"seqtk subseq {file_name} {output_file} > {final}"
        print(f"running {cmd}")   
        #Execute the subprocess command and capture the output
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()
    
        seqtk_result = ""
        # Check if the command was successful
        if process.returncode == 0:
        # Store the output in the results dictionary
            seqtk_result = output.decode('utf-8')
        else:
        # Print an error message if the command failed
            print(f"Error processing {output_file}: {error.decode('utf-8')}")




def main():

    # Usage: python3 script.py -i input_fastq_file
    # Create ArgumentParser object
    parser = argparse.ArgumentParser(description="Tool processing fastq file to divide it randomly in two files")
    
    # Add arguments
    parser.add_argument('-i', '--input_fastq_file', required=True, help='Fastq file to be randomly divided.')

    # Parse command-line arguments
    args = parser.parse_args()

    # Access arguments
    file_name = args.input_fastq_file
    
    random.seed(4)
    number = 2 #Into how many parts we want to divide the original file
    sequences = read_fastq(file_name)
    split_fastq(sequences, number, file_name)
    call_seqtk(file_name, number)
    print("FASTQ file successfully split into two files:", output_file1, "and", output_file2)


if __name__ == "__main__":
    main()

