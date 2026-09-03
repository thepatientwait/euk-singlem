# init
# ----
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import subprocess
import tempfile
import os

# logging init
# ------------
from log_config import configure_logger

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)

# snakefile parameters
# --------------------
FAA_FILE = snakemake.input['faa']
GENE_ID = snakemake.wildcards['gene_id']

VALIDATION_OUTPUT = snakemake.output['validation']

MIN_SEQUENCES = snakemake.params.get('min_sequences', 5)
MIN_ALIGNMENT_LENGTH = snakemake.params.get('min_alignment_length', 100)
MAX_GAP_FRACTION = snakemake.params.get('max_gap_fraction', 0.5)
MIN_ALIGNMENT_QUALITY = snakemake.params.get('min_alignment_quality', 0.3)


# functions
# ---------
def count_sequences(faa_file):
    """Count number of sequences in fasta file"""
    count = 0
    with open(faa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def calculate_alignment_stats(alignment):
    """Calculate statistics for a multiple sequence alignment"""
    if len(alignment) == 0:
        return {
            'length': 0,
            'num_sequences': 0,
            'mean_gap_fraction': 1.0,
            'conservation_score': 0.0
        }
    
    aln_length = alignment.get_alignment_length()
    num_seqs = len(alignment)
    
    # Calculate gap fractions per sequence
    gap_fractions = []
    for record in alignment:
        gaps = str(record.seq).count('-')
        gap_fractions.append(gaps / aln_length)
    
    mean_gap_fraction = sum(gap_fractions) / len(gap_fractions)
    
    # Calculate conservation score (fraction of positions with >50% identical residues)
    conserved_positions = 0
    for i in range(aln_length):
        column = alignment[:, i]
        # Count most common non-gap character
        char_counts = {}
        for char in column:
            if char != '-':
                char_counts[char] = char_counts.get(char, 0) + 1
        
        if char_counts:
            max_count = max(char_counts.values())
            non_gap_count = sum(char_counts.values())
            if non_gap_count > 0 and max_count / non_gap_count > 0.5:
                conserved_positions += 1
    
    conservation_score = conserved_positions / aln_length if aln_length > 0 else 0
    
    return {
        'length': aln_length,
        'num_sequences': num_seqs,
        'mean_gap_fraction': mean_gap_fraction,
        'conservation_score': conservation_score
    }


def create_hmm_alignment(faa_file):
    """Create HMM alignment using hmmbuild and hmmalign"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # First create initial alignment with muscle or mafft
        initial_aln = os.path.join(tmpdir, 'initial.aln')
        try:
            # Try mafft first
            subprocess.run(
                ['mafft', '--auto', '--quiet', faa_file],
                stdout=open(initial_aln, 'w'),
                stderr=subprocess.DEVNULL,
                check=True
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("mafft failed or not available, trying muscle")
            try:
                subprocess.run(
                    ['muscle', '-in', faa_file, '-out', initial_aln, '-quiet'],
                    check=True,
                    stderr=subprocess.DEVNULL
                )
            except (subprocess.CalledProcessError, FileNotFoundError):
                logger.error("Both mafft and muscle failed, cannot validate alignment")
                return None
        
        # Build HMM from initial alignment
        hmm_file = os.path.join(tmpdir, 'gene.hmm')
        try:
            subprocess.run(
                ['hmmbuild', '--amino', hmm_file, initial_aln],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
        except subprocess.CalledProcessError:
            logger.error("hmmbuild failed")
            return None
        
        # Align sequences to HMM
        stockholm_aln = os.path.join(tmpdir, 'hmm.sto')
        try:
            subprocess.run(
                ['hmmalign', '--amino', '-o', stockholm_aln, hmm_file, faa_file],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
        except subprocess.CalledProcessError:
            logger.error("hmmalign failed")
            return None
        
        # Read stockholm alignment
        try:
            alignment = AlignIO.read(stockholm_aln, 'stockholm')
            return alignment
        except Exception as e:
            logger.error(f"Failed to read alignment: {e}")
            return None


def validate_gene(faa_file, gene_id):
    """Validate if gene sequences are suitable for SingleM package"""
    validation_result = {
        'gene_id': gene_id,
        'num_sequences': 0,
        'passes': False,
        'failures': [],
        'warnings': [],
        'alignment_length': 0,
        'mean_gap_fraction': 0.0,
        'conservation_score': 0.0
    }
    
    # Check 1: Minimum number of sequences
    num_seqs = count_sequences(faa_file)
    validation_result['num_sequences'] = num_seqs
    
    if num_seqs < MIN_SEQUENCES:
        validation_result['failures'].append(f"Too few sequences ({num_seqs} < {MIN_SEQUENCES})")
        return validation_result
    
    # Check 2: Create HMM-based alignment and validate
    logger.info(f"Creating HMM alignment for {gene_id}")
    alignment = create_hmm_alignment(faa_file)
    
    if alignment is None:
        validation_result['failures'].append("Failed to create HMM alignment")
        return validation_result
    
    stats = calculate_alignment_stats(alignment)
    validation_result['alignment_length'] = stats['length']
    validation_result['mean_gap_fraction'] = stats['mean_gap_fraction']
    validation_result['conservation_score'] = stats['conservation_score']
    
    # Check 3: Minimum alignment length
    if stats['length'] < MIN_ALIGNMENT_LENGTH:
        validation_result['failures'].append(
            f"Alignment too short ({stats['length']} < {MIN_ALIGNMENT_LENGTH})"
        )
    
    # Check 4: Maximum gap fraction
    if stats['mean_gap_fraction'] > MAX_GAP_FRACTION:
        validation_result['failures'].append(
            f"Too many gaps (mean {stats['mean_gap_fraction']:.2%} > {MAX_GAP_FRACTION:.0%})"
        )
    
    # Check 5: Minimum conservation
    if stats['conservation_score'] < MIN_ALIGNMENT_QUALITY:
        validation_result['failures'].append(
            f"Poor conservation ({stats['conservation_score']:.2%} < {MIN_ALIGNMENT_QUALITY:.0%})"
        )
    
    # Warnings for borderline cases
    if MIN_ALIGNMENT_QUALITY <= stats['conservation_score'] < MIN_ALIGNMENT_QUALITY + 0.1:
        validation_result['warnings'].append("Low conservation score")
    
    if MAX_GAP_FRACTION - 0.1 < stats['mean_gap_fraction'] <= MAX_GAP_FRACTION:
        validation_result['warnings'].append("High gap fraction")
    
    # Gene passes if no failures
    validation_result['passes'] = len(validation_result['failures']) == 0
    
    return validation_result


# main
# ----
def main():
    logger.info(f"Validating gene sequences for {GENE_ID}")
    
    result = validate_gene(FAA_FILE, GENE_ID)
    
    # Write validation result
    with open(VALIDATION_OUTPUT, 'w') as f:
        f.write(f"gene_id\t{result['gene_id']}\n")
        f.write(f"num_sequences\t{result['num_sequences']}\n")
        f.write(f"alignment_length\t{result['alignment_length']}\n")
        f.write(f"mean_gap_fraction\t{result['mean_gap_fraction']:.4f}\n")
        f.write(f"conservation_score\t{result['conservation_score']:.4f}\n")
        f.write(f"passes\t{result['passes']}\n")
        
        if result['failures']:
            f.write(f"failures\t{'; '.join(result['failures'])}\n")
        
        if result['warnings']:
            f.write(f"warnings\t{'; '.join(result['warnings'])}\n")
    
    if result['passes']:
        logger.info(f"✓ {GENE_ID} passed validation")
    else:
        logger.warning(f"✗ {GENE_ID} failed validation: {', '.join(result['failures'])}")


if __name__ == '__main__':
    main()
