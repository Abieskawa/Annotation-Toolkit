#!/usr/bin/env python3
import sys

#####<<<Define functions>>>#####

def Count(data):
    """
    Count the sequence length by removing newline characters.
    data: the chunk of a FASTA entry (header + sequence lines).
    """
    # Find the first newline (end of header)
    v = data.find('\n')
    # Sequence lines start after the header
    x = data[v+1:]
    # Count how many newline characters are in the sequence lines
    y = x.count('\n')
    # Return length of the sequence lines after removing newlines
    return len(x) - y

def AverageSequenceLength(sequencelengthlist):
    """
    Return the average of values in a list of sequence lengths.
    """
    return float(sum(sequencelengthlist)) / float(len(sequencelengthlist))

def Min(sequencelengthlist):
    """
    Return the minimum value in a list of sequence lengths.
    """
    sequencelengthlistsorted = sorted(sequencelengthlist)
    return sequencelengthlistsorted[0] if sequencelengthlistsorted else 0

def Max(sequencelengthlist):
    """
    Return the maximum value in a list of sequence lengths.
    """
    sequencelengthlistsorted = sorted(sequencelengthlist)
    return sequencelengthlistsorted[-1] if sequencelengthlistsorted else 0

def NX0(seqlenlist, cutoff):
    """
    Compute the Nx statistic:
    - Nx length: length of the scaffold/contig at the position 
      where the cumulative length reaches 'cutoff' * sum(all lengths).
    - Lx: index (count) of how many sequences needed to reach that cutoff.
    """
    # Sort descending
    seqlist_sorted = sorted(seqlenlist, reverse=True)
    count = 0
    total_len_sum = 0
    Nx_length = 0
    Lx_count = 0
    
    total_sum = sum(seqlist_sorted)
    for i, length in enumerate(seqlist_sorted):
        total_len_sum += length
        if float(total_len_sum) / float(total_sum) >= cutoff:
            Nx_length = length
            Lx_count = i + 1
            break
    
    return [Nx_length, Lx_count]

def Percent(cutoff, seqlenlist, datalist):
    """
    Calculate:
      - # of sequences >= cutoff
      - sum of bases in those sequences
      - % of total length (>= cutoff vs all)
      - # of N's and % of N's
      - total bases without N for those sequences
    """
    results = []
    
    # Sort lengths ascending
    seqlenlistsorted = sorted(seqlenlist)
    # Sort datalist by length of the chunk (like original code)
    datalistsorted = sorted(datalist, key=len)
    
    count = 0
    temp_index = -1
    
    # Find first index where length >= cutoff
    for i, length in enumerate(seqlenlistsorted):
        if length >= cutoff:
            temp_index = i
            break
    
    if temp_index == -1:
        # Means no sequences >= cutoff
        Numofscaffold = 0
        Totalbase = 0
        percentage = 0
        TotalN = 0
        Npercent = 0
    else:
        Numofscaffold = len(seqlenlistsorted[temp_index:])
        Totalbase = sum(seqlenlistsorted[temp_index:])
        percentage = float(Totalbase) / float(sum(seqlenlistsorted)) * 100
        
        # Join the sequences (already sorted by length) that are >= cutoff
        datalist_join = ''.join(datalistsorted[temp_index:])
        TotalN = datalist_join.count('N')
        Npercent = float(TotalN) / float(Totalbase) * 100
    
    TotalbaseNoN = Totalbase - TotalN
    
    results.append(Numofscaffold)
    results.append(Totalbase)
    results.append(percentage)
    results.append(TotalN)
    results.append(Npercent)
    results.append(TotalbaseNoN)
    return results

#####<<<Main>>>#####

def main():
    # usage: python3 assembly_quality_check_v5.py [Input FASTA file] [data type (optional)]
    
    if len(sys.argv) < 2:
        print("Usage: python3 assembly_quality_check_v5.py [FASTA file] [data type (optional)]")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    
    if len(sys.argv) < 3:
        datatype = 'Scaffold'
    else:
        datatype = sys.argv[2]
    
    # Read the file using readline in a loop
    rawdata = ""
    with open(input_fasta, 'r') as f:
        while True:
            line = f.readline()
            if not line:  # End of file
                break
            rawdata += line
    
    # Count number of sequences (i.e. '>' lines)
    seqnu = rawdata.count('>')
    print(f'Number of {datatype}: {seqnu}')
    SeqNum = f'Number of {datatype}: {seqnu}'
    
    # Split into chunks by '>'
    datalist = rawdata.split('>')
    
    # Collect lengths
    seq_lengths = []
    datanum = 1
    while datanum <= seqnu:
        b = datalist[datanum]
        c = Count(b)
        seq_lengths.append(c)
        datanum += 1
    
    # Prepare output file
    Newname = input_fasta + '.quality.txt'
    with open(Newname, 'w') as output:
        # Write number of sequences
        output.write(SeqNum + '\n')
        
        # Average
        ave_len = AverageSequenceLength(seq_lengths)
        ASL = f"Average {datatype} Length : {ave_len}"
        print(ASL)
        output.write(ASL + '\n')
        
        # Min
        mini = Min(seq_lengths)
        Mini = f"Minimum {datatype} Length : {mini}"
        print(Mini)
        output.write(Mini + '\n')
        
        # Max
        maxi = Max(seq_lengths)
        Maxi = f"Maximum {datatype} Length : {maxi}"
        print(Maxi)
        output.write(Maxi + '\n')
        
        # N50, L50
        N50_val, L50_val = NX0(seq_lengths, 0.5)
        N50length = f"N50 (bp) : {N50_val}\nL50 : {L50_val}"
        print(N50length)
        output.write(N50length + '\n')
        
        # N75, L75
        N75_val, L75_val = NX0(seq_lengths, 0.75)
        N75length = f"N75 (bp) : {N75_val}\nL75 : {L75_val}"
        print(N75length)
        output.write(N75length + '\n')
        
        # N90, L90
        N90_val, L90_val = NX0(seq_lengths, 0.9)
        N90length = f"N90 (bp) : {N90_val}\nL90 : {L90_val}"
        print(N90length)
        output.write(N90length + '\n')
        
        # N95, L95
        N95_val, L95_val = NX0(seq_lengths, 0.95)
        N95length = f"N95 (bp) : {N95_val}\nL95 : {L95_val}"
        print(N95length)
        output.write(N95length + '\n')
        
        # N99, L99
        N99_val, L99_val = NX0(seq_lengths, 0.99)
        N99length = f"N99 (bp) : {N99_val}\nL99 : {L99_val}"
        print(N99length)
        output.write(N99length + '\n')
        
        # >= 300bp
        per300ans = Percent(300, seq_lengths, datalist)
        per300 = (
            f"Number of {datatype} >= 300bp : {per300ans[0]}\n"
            f"Total bases in {datatype}s >= 300bp : {per300ans[1]}\n"
            f"The percentage of {datatype} length >= 300bp : {per300ans[2]}\n"
            f"Number of N in >= 300bp : {per300ans[3]}\n"
            f"N-percent in >= 300bp : {per300ans[4]}\n"
            f"Total bases without N in {datatype}s >= 300bp : {per300ans[5]}"
        )
        print(per300)
        output.write(per300 + '\n')
        
        # >= 100kb
        per100kbans = Percent(100000, seq_lengths, datalist)
        per100kb = (
            f"Number of {datatype} >= 100kb : {per100kbans[0]}\n"
            f"Total bases in {datatype}s >= 100kb : {per100kbans[1]}\n"
            f"The percentage of {datatype} length >= 100kb : {per100kbans[2]}\n"
            f"Number of N in >= 100kb : {per100kbans[3]}\n"
            f"N-percent in >= 100kb : {per100kbans[4]}\n"
            f"Total bases without N in {datatype}s >= 100kb : {per100kbans[5]}"
        )
        print(per100kb)
        output.write(per100kb + '\n')
        
        # >= 1Mb
        per1Mans = Percent(1000000, seq_lengths, datalist)
        per1M = (
            f"Number of {datatype} >= 1Mb : {per1Mans[0]}\n"
            f"Total bases in {datatype}s >= 1Mb : {per1Mans[1]}\n"
            f"The percentage of {datatype} length >= 1Mb : {per1Mans[2]}\n"
            f"Number of N in >= 1Mb : {per1Mans[3]}\n"
            f"N-percent in >= 1Mb : {per1Mans[4]}\n"
            f"Total bases without N in {datatype}s >= 1Mb : {per1Mans[5]}"
        )
        print(per1M)
        output.write(per1M + '\n')


if __name__ == "__main__":
    main()
