#!/usr/bin/env python3

import polars as pl
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
import re
import sys
from joblib import Parallel, delayed

# Import the logging configuration
sys.path.append(str(Path(__file__).parent))
from log_config import configure_logger

# Configure logger
logger = configure_logger(__file__)



def parse_busco_file(busco_file: Path, use_multi_copy: bool) -> Tuple[str, str, pl.DataFrame]:
    """
    Parse a BUSCO results file and extract group, genome_acc, and complete genes.
    
    Args:
        busco_file: Path to BUSCO full_table.tsv file
        
    Returns:
        Tuple of (group, genome_acc, complete_genes_dataframe)
    """
    # Extract group and genome accession from path
    # Path format: outputs/busco/{group}/{genome_acc}/run_eukaryota_odb12/full_table.tsv
    path_parts = busco_file.parts
    group = path_parts[-4]
    genome_acc = path_parts[-3]
    
    logger.info(f"Parsing BUSCO file for {group}/{genome_acc}")
    
    # Read BUSCO full table, skipping comment lines that start with #
    busco_df = pl.read_csv(busco_file, separator='\t', skip_rows=2).rename({'# Busco id': 'Busco id'})
    
    # Debug: log the column names to understand the structure
    logger.debug(f"BUSCO file columns: {busco_df.columns}")
    
    # Validate that required columns exist
    if 'Status' not in busco_df.columns:
        raise ValueError(f"Missing 'Status' column in {busco_file}. Found columns: {busco_df.columns}")
    
    if 'Busco id' not in busco_df.columns:
        raise ValueError(f"Missing 'Busco id' column in {busco_file}. Found columns: {busco_df.columns}")
    
    if 'Score' not in busco_df.columns:
        raise ValueError(f"Missing 'Score' column in {busco_file}. Found columns: {busco_df.columns}")

    # Keep every complete or duplicated hit, with copy 1 being the best-scoring hit.
    selected_genes = (
        busco_df
        .filter(pl.col('Status').is_in(['Complete'] + (['Duplicated'] if use_multi_copy else [])))
        .with_columns(pl.col('Score').cast(pl.Float64, strict=False))
        .sort(['Busco id', 'Score'], descending=[False, True], nulls_last=True)
        .unique(subset=['Busco id'], keep='first')
    )
    
    logger.info(f"Found {len(selected_genes)} complete/duplicated genes for {group}/{genome_acc}")
    return group, genome_acc, selected_genes


def extract_gene_sequences(busco_file: Path, selected_genes: pl.DataFrame,
                          group: str, genome_acc: str) -> Dict[str, Dict]:
    """
    Extract gene sequences from BUSCO output directory.
    
    Args:
        busco_file: Path to BUSCO full_table.tsv file
        selected_genes: DataFrame of complete and duplicated BUSCO genes,
            ordered by BUSCO ID and descending score
        group: Organism group name
        genome_acc: Genome accession identifier
        
    Returns:
        Dictionary mapping busco_id to sequence information
    """
    seq_dirs = [
        busco_file.parent / "busco_sequences" / "single_copy_busco_sequences",
        busco_file.parent / "busco_sequences" / "multi_copy_busco_sequences",
    ]
    
    existing_seq_dirs = [seq_dir for seq_dir in seq_dirs if seq_dir.exists()]
    if not existing_seq_dirs:
        raise FileNotFoundError(f"Sequences directories not found: {seq_dirs}")
    
    gene_data = {}
    sequences_found = 0
    
    for busco_id in selected_genes['Busco id'].unique().to_list():
        seq_file = next(
            (seq_dir / f"{busco_id}.faa" for seq_dir in existing_seq_dirs
             if (seq_dir / f"{busco_id}.faa").exists()),
            None,
        )
        if seq_file is None:
            logger.debug(f"Sequence file not found for {busco_id}")
            continue

        records = list(SeqIO.parse(seq_file, 'fasta'))
        busco_rows = selected_genes.filter(pl.col('Busco id') == busco_id)
        used_record_ids = set()

        for copy_number, row in enumerate(busco_rows.iter_rows(named=True), start=1):
            record = _match_busco_record(records, row, used_record_ids)
            if record is None:
                logger.warning(f"No FASTA record matched {busco_id} copy {copy_number} in {seq_file}")
                continue

            used_record_ids.add(record.id)
            sequence = f">{record.description}\n{record.seq}"
            sequence_id = f"{genome_acc}|{busco_id}|{copy_number}"
            gene_data[sequence_id] = {
                'sequence': sequence,
                'sequence_id': sequence_id,
                'group': group,
                'genome_acc': genome_acc,
                'busco_id': busco_id,
                'copy_number': copy_number,
                'score': row['Score'],
                'length': len(record.seq),
            }
            sequences_found += 1
    
    logger.info(f"Extracted {sequences_found} gene sequences for {group}/{genome_acc}")
    return gene_data


def _match_busco_record(records, row: Dict, used_record_ids: set):
    """Match a BUSCO table row to its FASTA record using sequence ID syntax."""
    sequence_name = str(row.get('Sequence', ''))
    start = row.get('Gene Start')
    end = row.get('Gene End')
    strand = str(row.get('Strand', ''))
    try:
        coordinate_key = (sequence_name, min(int(start), int(end)), max(int(start), int(end)), strand)
    except (TypeError, ValueError):
        coordinate_key = None

    for record in records:
        if record.id in used_record_ids:
            continue
        match = re.search(r'\|([^|:]+):(\d+)-(\d+)\|([+-])(?:\s|$)', record.description)
        if match is None:
            continue
        record_key = (match.group(1), min(int(match.group(2)), int(match.group(3))),
                      max(int(match.group(2)), int(match.group(3))), match.group(4))
        if coordinate_key == record_key:
            return record
    return None


def calculate_sequence_length(sequence: str) -> int:
    """Calculate the actual sequence length, handling FASTA format."""
    if '>' in sequence:
        # FASTA format: extract just the sequence part
        lines = sequence.strip().split('\n')
        return len(''.join(lines[1:]))
    else:
        # Plain sequence
        return len(sequence.replace('\n', ''))



def process_single_busco_file(busco_file_str: str, use_multi_copy: bool) -> Tuple[Dict[str, Dict], Optional[str]]:
    """
    Process a single BUSCO file and return results.
    
    Args:
        busco_file_str: String path to BUSCO full_table.tsv file
        
    Returns:
        Tuple of (gene_data_dict, error_message_or_None)
    """
    busco_file = Path(busco_file_str)
    
    # Check if file exists first
    if not busco_file.exists():
        error_msg = f"BUSCO file does not exist: {busco_file}"
        logger.error(error_msg)
        return {}, error_msg
    
    try:
        # Parse BUSCO results
        group, genome_acc, selected_genes = parse_busco_file(busco_file, use_multi_copy)
        
        # Extract gene sequences
        gene_data = extract_gene_sequences(busco_file, selected_genes, group, genome_acc)
        
        # Transform to the format expected by the aggregator
        # busco_id -> {sequence_id: gene_info}
        result = defaultdict(dict)
        for sequence_id, info in gene_data.items():
            result[info['busco_id']][sequence_id] = info
        
        return dict(result), None
        
    except Exception as e:
        error_msg = f"Critical error processing {busco_file}: {e}"
        logger.error(error_msg)
        return {}, error_msg


def merge_gene_sequences(results: List[Tuple[Dict[str, Dict], Optional[str]]]) -> Tuple[Dict[str, Dict], List[str]]:
    """
    Merge results from parallel processing into a single gene_sequences dictionary.
    
    Args:
        results: List of tuples from process_single_busco_file
        
    Returns:
        Tuple of (merged_gene_sequences, critical_errors)
    """
    gene_sequences = defaultdict(dict)
    critical_errors = []
    
    for gene_data, error_msg in results:
        if error_msg:
            critical_errors.append(error_msg)
        else:
            # Merge gene data
            for busco_id, taxon_data in gene_data.items():
                gene_sequences[busco_id].update(taxon_data)
    
    return dict(gene_sequences), critical_errors



def process_all_busco_results(busco_files: List[str], use_multi_copy: bool) -> Dict[str, Dict]:
    """
    Process all BUSCO result files and extract gene sequences using multiprocessing.
    
    Args:
        busco_files: List of paths to BUSCO full_table.tsv files
        
    Returns:
        Dictionary mapping busco_id to taxon sequences
    """
    logger.info(f"Processing {len(busco_files)} BUSCO result files")
    
    # Get number of threads from snakemake
    n_jobs = getattr(snakemake, 'threads', 1)
    logger.info(f"Using {n_jobs} parallel jobs for processing")
    
    # Process files in parallel
    logger.info("Starting parallel processing of BUSCO files...")
    results = Parallel(n_jobs=n_jobs, prefer="threads", verbose=1)(
        delayed(process_single_busco_file)(busco_file, use_multi_copy) 
        for busco_file in busco_files
    )
    
    # Merge results and collect errors
    gene_sequences, critical_errors = merge_gene_sequences(results)
    
    processed_count = len(busco_files) - len(critical_errors)
    
    # Fail if we have too many critical errors
    error_rate = len(critical_errors) / len(busco_files) if busco_files else 0
    if error_rate > 0.1:  # Fail if more than 10% of files have critical errors
        raise RuntimeError(f"Too many critical errors ({len(critical_errors)}/{len(busco_files)} = {error_rate:.1%}). "
                          f"First few errors: {critical_errors[:3]}")
    
    if critical_errors:
        logger.warning(f"Had {len(critical_errors)} critical errors, but continuing since error rate is acceptable ({error_rate:.1%})")
    
    logger.info(f"Processing complete: {processed_count} genomes processed, "
               f"{len(critical_errors)} errors, {len(gene_sequences)} unique genes found")
    
    return gene_sequences


def deduplicate_sequences(genome_sequences: Dict) -> Tuple[Dict, int]:
    """
    Deduplicate sequences based on sequence ID and content.
    
    Args:
        genome_sequences: Dictionary mapping genome_acc to sequence info
        
    Returns:
        Tuple of (deduplicated_sequences, duplicates_removed_count)
    """
    seen_sequences = set()
    deduplicated = {}
    duplicates_removed = 0
    
    for sequence_id, info in genome_sequences.items():
        sequence = info['sequence']
        
        # Extract sequence ID and content
        if '>' in sequence:
            lines = sequence.strip().split('\n')
            seq_content = ''.join(lines[1:])
        else:
            seq_content = sequence
        
        seq_id = info.get('sequence_id', f"{info['genome_acc']}|{info['busco_id']}|1")
        unique_key = (seq_id, seq_content)
        
        if unique_key not in seen_sequences:
            seen_sequences.add(unique_key)
            # Store with the prefixed sequence ID as the key for consistency
            deduplicated[seq_id] = {
                **info,
                'original_genome_acc': info['genome_acc'],
                'sequence_id': seq_id
            }
        else:
            duplicates_removed += 1
            
    return deduplicated, duplicates_removed


USABLE_TAX_STATUSES = ('complete', 'partial')


def load_usable_taxonomy(taxonomy_df: pl.DataFrame) -> Dict[str, str]:
    """Return {genome acc: taxonomy string} for genomes with a usable taxonomy.

    "Usable" is defined once, by ``tax_status`` in the mapping written by
    fetch_ncbi_taxonomy.py; genomes absent from the mapping have no taxonomy.
    """
    required = {'acc', 'parsed_taxon_string', 'tax_status'}
    missing = required - set(taxonomy_df.columns)
    if missing:
        raise ValueError(
            f"Taxonomy mapping is missing columns {sorted(missing)}; "
            "regenerate it with the fetch_ncbi_taxonomy rule"
        )
    usable = taxonomy_df.filter(pl.col('tax_status').is_in(USABLE_TAX_STATUSES))
    return dict(zip(usable['acc'].cast(pl.Utf8).to_list(), usable['parsed_taxon_string'].to_list()))


def resolve_sequence_taxonomy(seq_id: str, info: Dict, usable_taxonomy: Dict[str, str]) -> Optional[str]:
    """Return the taxonomy string for a sequence's genome, or None if it has no usable taxonomy."""
    return usable_taxonomy.get(str(info.get('original_genome_acc', seq_id)))


def filter_resolved_sequences(deduplicated_sequences: Dict, usable_taxonomy: Dict[str, str],
                              busco_id: str, allow_missing_taxonomy: bool = False) -> Dict:
    """
    Keep only sequences whose genome has a usable taxonomy entry.

    FASTA and taxonomy outputs must contain the same sequence set, otherwise
    downstream tools like graftM create fail with 'Unable to find sequence in
    the taxonomy definition'.

    If allow_missing_taxonomy is True (e.g. tree building, which does not need
    per-sequence taxonomy), all sequences are kept. The per-gene taxonomy TSV
    will then only list sequences that do have a usable taxonomy.
    """
    if allow_missing_taxonomy:
        return dict(deduplicated_sequences)

    kept = {}
    skipped = 0

    for seq_id, info in deduplicated_sequences.items():
        if resolve_sequence_taxonomy(seq_id, info, usable_taxonomy) is None:
            skipped += 1
            continue
        kept[seq_id] = info

    if skipped > 0:
        logger.info(
            f"{busco_id}: excluding {skipped} sequences without taxonomy "
            f"from FASTA/taxonomy outputs ({len(kept)} remaining)"
        )

    return kept


def write_gene_fasta(busco_id: str, deduplicated_sequences: Dict, output_dir: Path) -> Path:
    """
    Write FASTA file for a single gene with pre-deduplicated sequences.
    
    Args:
        busco_id: BUSCO gene identifier
        deduplicated_sequences: Pre-deduplicated sequences keyed by sequence_id
        output_dir: Output directory path
        
    Returns:
        Path to written FASTA file
    """
    fasta_file = output_dir / f"{busco_id}.faa"
    
    with open(fasta_file, 'w') as f:
        for seq_id, info in deduplicated_sequences.items():
            sequence = info['sequence']
            
            # Extract just the sequence content, ignoring original header
            if '>' in sequence:
                lines = sequence.strip().split('\n')
                seq_content = ''.join(lines[1:])
            else:
                seq_content = sequence
            
            # Always write with the prefixed seq_id
            f.write(f">{seq_id}\n{seq_content}\n")
    
    return fasta_file


def write_taxonomy_file(busco_id: str, deduplicated_sequences: Dict, 
                       usable_taxonomy: Dict[str, str], output_dir: Path) -> Path:
    """
    Write taxonomy mapping file for pre-deduplicated sequences.
    
    Args:
        busco_id: BUSCO gene identifier
        deduplicated_sequences: Pre-deduplicated sequences keyed by sequence_id
        usable_taxonomy: {genome acc: taxonomy string} for genomes with a usable taxonomy
        output_dir: Output directory path
        
    Returns:
        Path to written taxonomy file
    """
    tax_file = output_dir / f"{busco_id}_taxonomy.tsv"
    written = 0
    skipped = 0

    with open(tax_file, 'w') as f:
        f.write("sequence_id\ttaxonomy\n")

        for seq_id, info in deduplicated_sequences.items():
            taxonomy = resolve_sequence_taxonomy(seq_id, info, usable_taxonomy)
            if taxonomy is None:
                skipped += 1
                logger.debug(f"Skipping {seq_id}: no usable taxonomy")
                continue

            # Write using the deduplicated sequence ID
            f.write(f"{seq_id}\t{taxonomy}\n")
            written += 1

    if skipped > 0:
        logger.info(f"{busco_id}: wrote {written} taxonomy rows and skipped {skipped} unresolved/shallow rows")
    
    return tax_file


def reorganize_for_singlem_genome_input(
    gene_sequences: Dict[str, Dict],
    usable_taxonomy: Dict[str, str],
    output_dir: Path,
    min_sequences_per_gene: int,
    allow_missing_taxonomy: bool = False,
) -> Tuple[List[str], Path, Path]:
    """
    Reorganize BUSCO data for the singlem genome input workflow:
    1. Create genome-centric FASTA files (instead of gene-centric)
    2. Create list of genome protein files
    3. Create HMM names file for BUSCO genes
    
    Args:
        gene_sequences: Dictionary mapping busco_id to taxon sequences
        usable_taxonomy: {genome acc: taxonomy string} for genomes with a usable taxonomy
        output_dir: Output directory path
        min_sequences_per_gene: Minimum number of sequences required per genome
        allow_missing_taxonomy: Keep genomes without a usable taxonomy in the
            genome FASTAs and genome list (they are still omitted from
            genome_taxonomy.tsv)
        
    Returns:
        Tuple of (processed_genomes, genome_file_list_path, hmms_and_names_path)
    """
    logger.info("Reorganizing BUSCO sequences for singlem genome input format")
    
    # Create subdirectories
    genomes_dir = output_dir / "genomes"
    genomes_dir.mkdir(parents=True, exist_ok=True)
    
    # Reorganize by genome accession
    genome_sequences = defaultdict(dict)  # genome_id -> {busco_id: sequence_info}
    
    for busco_id, taxon_seqs in gene_sequences.items():
        for sequence_id, seq_info in taxon_seqs.items():
            genome_acc = seq_info['genome_acc']
            genome_sequences[genome_acc][sequence_id] = seq_info
    
    logger.info(f"Found {len(genome_sequences)} genomes with BUSCO sequences")

    # Genomes without a usable taxonomy are excluded from the genome FASTAs, genome list,
    # and genome taxonomy file alike, unless allow_missing_taxonomy is set.
    resolved_accs = set(usable_taxonomy)

    # Write genome-centric FASTA files
    processed_genomes = []
    genome_file_paths = []
    skipped_no_taxonomy = 0

    for genome_id, busco_seqs in genome_sequences.items():
        if not allow_missing_taxonomy and genome_id not in resolved_accs:
            skipped_no_taxonomy += 1
            continue
        if len(busco_seqs) >= min_sequences_per_gene:
            genome_file = genomes_dir / f"{genome_id}.faa"
            
            with open(genome_file, 'w') as f:
                for sequence_id, seq_info in busco_seqs.items():
                    sequence = seq_info['sequence']
                    
                    # Extract sequence content
                    if '>' in sequence:
                        lines = sequence.strip().split('\n')
                        seq_content = ''.join(lines[1:])
                    else:
                        seq_content = sequence
                    
                    # Create sequence ID: genome|busco_gene|copy_number
                    seq_id = seq_info['sequence_id']
                    f.write(f">{seq_id}\n{seq_content}\n")
            
            processed_genomes.append(genome_id)
            genome_file_paths.append(str(genome_file.absolute()))
            logger.debug(f"Wrote {len(busco_seqs)} sequences for genome {genome_id}")
        else:
            logger.debug(f"Skipping genome {genome_id} (only {len(busco_seqs)} sequences)")

    if skipped_no_taxonomy > 0:
        logger.info(
            f"Excluded {skipped_no_taxonomy} genomes without taxonomy from genome-input outputs"
        )

    # Write genome file list
    genome_list_file = output_dir / "genome_protein_files.txt"
    with open(genome_list_file, 'w') as f:
        for file_path in genome_file_paths:
            f.write(f"{file_path}\n")
    
    # Create HMMs and names file for BUSCO genes
    hmms_and_names_file = output_dir / "hmms_and_names_busco.tsv"
    busco_hmm_dir = Path("busco_downloads/lineages/eukaryota_odb12/hmms")
    
    with open(hmms_and_names_file, 'w') as f:
        f.write("gene\tname\tfilepath\tdescription\n")
        
        for busco_id in gene_sequences.keys():
            hmm_file = busco_hmm_dir / f"{busco_id}.hmm"
            if hmm_file.exists():
                f.write(f"{busco_id}\t{busco_id}\t{hmm_file.absolute()}\t{busco_id}\n")
    
    # Create taxonomy file for genomes (only genomes with a usable taxonomy are listed)
    genome_taxonomy_file = output_dir / "genome_taxonomy.tsv"
    kept_genomes = 0
    with open(genome_taxonomy_file, 'w') as f:
        f.write("genome_id\ttaxonomy\n")
        for genome_id in processed_genomes:
            taxonomy = usable_taxonomy.get(str(genome_id))
            if taxonomy is None:
                continue
            f.write(f"{genome_id}\t{taxonomy}\n")
            kept_genomes += 1

    skipped_genomes = len(processed_genomes) - kept_genomes
    if skipped_genomes > 0:
        logger.info(f"genome_taxonomy.tsv: wrote {kept_genomes} genomes and skipped {skipped_genomes} without usable taxonomy")

    logger.info(f"SingleM genome-input reorganization complete: {len(processed_genomes)} genomes processed")
    return processed_genomes, genome_list_file, hmms_and_names_file


def write_gene_files(gene_sequences: Dict[str, Dict], usable_taxonomy: Dict[str, str], 
                    output_dir: Path, min_sequences: int = 5,
                    allow_missing_taxonomy: bool = False) -> List[str]:
    """
    Write FASTA and taxonomy files for genes with sufficient sequences.
    
    Args:
        gene_sequences: Dictionary mapping busco_id to taxon sequences
        usable_taxonomy: {genome acc: taxonomy string} for genomes with a usable taxonomy
        output_dir: Output directory path
        min_sequences: Minimum number of sequences required per gene
        allow_missing_taxonomy: Keep sequences whose genome has no usable taxonomy
        
    Returns:
        List of processed gene IDs
    """
    processed_genes = []
    
    logger.info(f"Writing gene files (minimum {min_sequences} sequences per gene)")
    
    for busco_id, taxon_seqs in gene_sequences.items():
        if len(taxon_seqs) >= min_sequences:
            try:
                # Deduplicate sequences once
                deduplicated_seqs, duplicates_removed = deduplicate_sequences(taxon_seqs)
                
                if duplicates_removed > 0:
                    logger.info(f"Removed {duplicates_removed} duplicate sequences from {busco_id} "
                               f"({len(deduplicated_seqs)} unique sequences remaining)")

                # Keep only sequences with a usable taxonomy so the FASTA and
                # taxonomy TSV always describe the same sequence set.
                deduplicated_seqs = filter_resolved_sequences(
                    deduplicated_seqs, usable_taxonomy, busco_id, allow_missing_taxonomy
                )
                if not deduplicated_seqs:
                    logger.warning(f"Skipping {busco_id}: no sequences with resolved taxonomy")
                    continue
                
                # Write FASTA file with deduplicated sequences
                fasta_file = write_gene_fasta(busco_id, deduplicated_seqs, output_dir)
                
                # Write taxonomy file with same deduplicated sequences
                tax_file = write_taxonomy_file(busco_id, deduplicated_seqs, usable_taxonomy, output_dir)
                
                processed_genes.append(busco_id)
                logger.debug(f"Wrote files for {busco_id} ({len(deduplicated_seqs)} unique sequences)")
                
            except Exception as e:
                logger.error(f"Error writing files for {busco_id}: {e}")
        else:
            logger.debug(f"Skipping {busco_id} (only {len(taxon_seqs)} sequences)")
    
    logger.info(f"Successfully wrote files for {len(processed_genes)} genes")
    return processed_genes





def main():
    """Main function to process BUSCO results and create SingleM input files."""
    logger.info("Starting BUSCO sequence aggregation")
    
    # Get inputs from snakemake
    manifest_file = Path(snakemake.input.busco_manifest)
    output_dir = Path(snakemake.output.sequences_dir)
    taxonomy_file = snakemake.input.taxonomy
    use_multi_copy = snakemake.params.get('use_multi_copy', False)
    allow_missing_taxonomy = bool(snakemake.params.get('allow_missing_taxonomy', False))
    if allow_missing_taxonomy:
        logger.info("allow_missing_taxonomy enabled: keeping sequences/genomes without taxonomy")

    if not manifest_file.exists():
        raise FileNotFoundError(f"BUSCO manifest file not found: {manifest_file}")

    manifest_df = pl.read_csv(manifest_file, separator='\t')
    if 'busco_result' not in manifest_df.columns:
        raise ValueError(f"Manifest missing 'busco_result' column: {manifest_file}")

    busco_results = sorted(manifest_df['busco_result'].to_list())
    
    # Validate inputs exist
    if not busco_results:
        raise ValueError(f"No BUSCO result files listed in manifest: {manifest_file}")
    
    if not Path(taxonomy_file).exists():
        raise FileNotFoundError(f"Taxonomy file not found: {taxonomy_file}")
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load taxonomy mapping
    logger.info(f"Loading taxonomy mapping from {taxonomy_file}")
    taxonomy_df = pl.read_csv(taxonomy_file, separator='\t')
    usable_taxonomy = load_usable_taxonomy(taxonomy_df)
    logger.info(f"Loaded taxonomy mapping with {len(taxonomy_df)} entries "
                f"({len(usable_taxonomy)} with usable taxonomy)")
    
    # Process all BUSCO results - let errors bubble up
    gene_sequences = process_all_busco_results(busco_results, use_multi_copy)
    
    # Validate we got meaningful results
    if not gene_sequences:
        raise RuntimeError("No gene sequences were successfully processed from any BUSCO files")
    
    # Write traditional gene-centric files (for backward compatibility)
    processed_genes = write_gene_files(
        gene_sequences, usable_taxonomy, output_dir,
        allow_missing_taxonomy=allow_missing_taxonomy,
    )
    
    # Reorganize for singlem genome-input workflow format
    min_sequences_per_gene = snakemake.params.get('min_sequences_per_gene', 5)
    genome_input_dir = Path(snakemake.output.genome_input_genomes).parent
    processed_genomes, genome_list_file, hmms_file = reorganize_for_singlem_genome_input(
        gene_sequences, usable_taxonomy, genome_input_dir, min_sequences_per_gene,
        allow_missing_taxonomy=allow_missing_taxonomy,
    )
    
    logger.info(f"Successfully completed aggregation: "
               f"{len(processed_genes)} genes and {len(processed_genomes)} genomes processed")
    logger.info(f"SingleM genome-input files created in: {genome_input_dir}")
    logger.info(f"Genome file list: {genome_list_file}")
    logger.info(f"HMM names file: {hmms_file}")

if __name__ == "__main__":
    main()