def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement.get(base, 'N') for base in reversed(seq)])

def read_contigs(contig_file):
    contigs = {}
    current_name = None
    with open(contig_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_name = line[1:]
                contigs[current_name] = []
            else:
                contigs[current_name].append(line)
    # Merge multi-line sequences
    for name in contigs:
        contigs[name] = ''.join(contigs[name])
    return contigs

def process_order(order_file, contigs, output_file):
    scaffold_num = 0
    with open(order_file, 'r') as f, open(output_file, 'w') as out:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Split scaffold part and contigs
            parts = line.split(': ')
            if len(parts) != 2:
                continue  # Skip invalid lines
            scaffold_part, contigs_part = parts
            contig_list = contigs_part.split()
            
            # Split into groups of two contiguous contigs
            groups = []
            for i in range(0, len(contig_list), 2):
                if i + 1 >= len(contig_list):
                    raise ValueError(f"Invalid contig pair in line: {line}")
                groups.append((contig_list[i], contig_list[i+1]))
            
            # Process each group to get the sequence
            scaffold_seqs = []
            for name1, name2 in groups:
                # Extract base contig name and suffixes
                base1, suffix1 = name1.rsplit('_', 1)
                base2, suffix2 = name2.rsplit('_', 1)
                if base1 != base2:
                    raise ValueError(f"Contig pair {name1} and {name2} are not from the same original contig.")
                original_name = base1
                if original_name not in contigs:
                    raise KeyError(f"Contig {original_name} not found.")
                seq = contigs[original_name]
                
                # Determine orientation
                if suffix1 == '1' and suffix2 == '2':
                    scaffold_seqs.append(seq)
                elif suffix1 == '2' and suffix2 == '1':
                    scaffold_seqs.append(reverse_complement(seq))
                else:
                    raise ValueError(f"Invalid suffixes: {suffix1} and {suffix2} in {name1} {name2}")
            
            # Join sequences with 500 N's between each contig group
            final_seq = ('N' * 10).join(scaffold_seqs)
            
            # Write to output
            scaffold_num += 1
            out.write(f">scaffold{scaffold_num}\n{final_seq}\n")

def main():
    import sys
    if len(sys.argv) != 4:
        print("Usage: python script.py <contig_file> <order_file> <output_file>")
        sys.exit(1)
    contig_file, order_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]
    contigs = read_contigs(contig_file)
    process_order(order_file, contigs, output_file)

if __name__ == "__main__":
    main()