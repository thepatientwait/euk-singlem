import sys
import argparse

def parse_fasta(filepath):
    sequences = {}
    current_id = None
    current_seq = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq).upper()
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = "".join(current_seq).upper()
    return sequences

def main():
    parser = argparse.ArgumentParser(description="Find top closest sequences to a target in a concatenated alignment.")
    parser.add_argument("fasta", help="Path to concatenated alignment FASTA file")
    parser.add_argument("target", help="Target genome ID to test")
    parser.add_argument("-n", "--top", type=int, default=10, help="Number of top matches to print (default: 10)")
    
    args = parser.parse_args()
    
    seqs = parse_fasta(args.fasta)
    if args.target not in seqs:
        print(f"Error: Target '{args.target}' not found in FASTA.")
        print(f"Sample IDs available: {list(seqs.keys())[:5]}...")
        sys.exit(1)
        
    target_seq = seqs[args.target]
    results = []
    
    for other_id, other_seq in seqs.items():
        if other_id == args.target:
            continue
        
        matches = 0
        comparable = 0
        for a, b in zip(target_seq, other_seq):
            if a != '-' and b != '-':
                comparable += 1
                if a == b:
                    matches += 1
                    
        identity = (matches / comparable) if comparable > 0 else 0
        results.append((other_id, identity))
        
    # Sort by highest identity
    results.sort(key=lambda x: x[1], reverse=True)
    
    print(f"\nTop {min(args.top, len(results))} closest matches to '{args.target}':")
    print("-" * 35)
    for gid, ident in results[:args.top]:
        print(f"{ident*100:6.2f}%  {gid}")

if __name__ == "__main__":
    main()
