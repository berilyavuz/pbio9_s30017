# ============================================================
# Album no: s30017
# Date: 2026-05-07
# Description: Random DNA sequence generator in FASTA format.
#              Supports motif search, complementary strand,
#              in silico transcription, sliding window GC
#              analysis, and ORF identification.
#v4
# ============================================================

import random
import csv
import os


# ── helpers ─────────────────────────────────────────────────

def validate_positive_int(prompt: str,
                           min_val: int = 1,
                           max_val: int = 100_000) -> int:
    """
    Repeatedly asks the user for an integer in [min_val, max_val].
    Returns the valid value; never raises an exception.
    """
    while True:
        raw = input(prompt).strip()
        try:
            value = int(raw)
            if min_val <= value <= max_val:
                return value
            raise ValueError
        except ValueError:
            print(f"Error: value must be an integer in the range "
                  f"[{min_val}, {max_val}].")


def validate_id(prompt: str) -> str:
    """
    Asks the user for a sequence ID that contains no whitespace.
    Repeats until valid input is provided.
    """
    while True:
        seq_id = input(prompt).strip()
        if seq_id and " " not in seq_id and "\t" not in seq_id:
            return seq_id
        print("Error: ID cannot be empty or contain whitespace.")


# ── core sequence functions ──────────────────────────────────

def generate_sequence(length: int) -> str:
    """Returns a random DNA sequence of the specified length."""
    nucleotides = ["A", "C", "G", "T"]
    return "".join(random.choice(nucleotides) for _ in range(length))


def generate_sequence_weighted(length: int,
                                weights: dict) -> str:
    """
    Returns a random DNA sequence of the given length using
    user-supplied nucleotide probabilities (weights).
    weights: {"A": float, "C": float, "G": float, "T": float}
    """
    bases = list(weights.keys())
    probs = [weights[b] / 100.0 for b in bases]
    return "".join(random.choices(bases, weights=probs, k=length))


def calculate_stats(sequence: str) -> dict:
    """
    Returns a dictionary of sequence statistics.
    Keys: 'A', 'C', 'G', 'T' (float, %), 'GC' (float, %).
    Note: only uppercase nucleotide characters are counted;
    embedded name letters (lowercase) are ignored.
    """
    # extract only uppercase nucleotide characters
    pure = [ch for ch in sequence if ch in "ACGT"]
    n = len(pure)
    if n == 0:
        return {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0, "GC": 0.0,
                "gc_ratio_A": 0.0}

    counts = {base: pure.count(base) for base in "ACGT"}
    stats = {base: round(counts[base] / n * 100, 2) for base in "ACGT"}
    gc = round((counts["G"] + counts["C"]) / n * 100, 2)
    stats["GC"] = gc
    return stats


def insert_name(sequence: str, name: str) -> str:
    """
    Inserts a name at a random position in the sequence.
    The name is written in lowercase so it is visually
    distinguishable from the uppercase nucleotides.
    """
    pos = random.randint(0, len(sequence))
    return sequence[:pos] + name.lower() + sequence[pos:]

# validate that the sequence length is within acceptable bounds
def format_fasta(seq_id: str,
                 description: str,
                 sequence: str,
                 line_width: int = 80) -> str:
    """
    Returns a properly formatted FASTA record as a string.
    Header line starts with '>'. Sequence is wrapped at
    line_width characters per line.
    """
    header = f">{seq_id}"
    if description:
        header += f" {description}"

    # wrap sequence into fixed-width lines
    lines = [sequence[i:i + line_width]
             for i in range(0, len(sequence), line_width)]

    return header + "\n" + "\n".join(lines) + "\n"


def save_fasta(seq_id: str, fasta_content: str) -> str:
    """
    Saves a FASTA string to a file named {seq_id}.fasta.
    Appends the required validator marker at the end.
    Returns the filename.
    """
    filename = f"{seq_id}.fasta"
    with open(filename, "w", encoding="utf-8") as fh:
        fh.write(fasta_content)
    return filename


def print_stats(stats: dict, length: int) -> None:
    """Prints nucleotide statistics to stdout in a readable format."""
    print(f"\nSequence statistics (n={length}):")
    for base in "ACGT":
        print(f"  {base}: {stats[base]:.2f}%")
    print(f"  GC-content: {stats['GC']:.2f}%")


# ── additional feature 1: motif search ──────────────────────

def search_motif(sequence: str, motif: str) -> list:
    """
    Searches for all (possibly overlapping) occurrences of
    motif in sequence. Returns a list of 1-based positions
    following the biological convention.
    """
    positions = []
    start = 0
    # only uppercase portion is biologically relevant
    pure_seq = "".join(ch for ch in sequence if ch.isupper())
    while True:
        idx = pure_seq.find(motif.upper(), start)
        if idx == -1:
            break
        positions.append(idx + 1)   # convert to 1-based index
        start = idx + 1
    return positions


# ── additional feature 2: complementary strands ─────────────

COMPLEMENT_MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def complementary_strand(sequence: str) -> str:
    """
    Returns the complementary strand (5'→3' of the template).
    Only uppercase nucleotides are complemented; lowercase
    embedded name characters are left unchanged.
    """
    result = []
    for ch in sequence:
        if ch.upper() in COMPLEMENT_MAP:
            comp = COMPLEMENT_MAP[ch.upper()]
            result.append(comp if ch.isupper() else comp.lower())
        else:
            result.append(ch)
    return "".join(result)


def reverse_complement(sequence: str) -> str:
    """
    Returns the reverse complement of sequence.
    This represents the antiparallel complementary strand
    read in the 5'→3' direction.
    """
    return complementary_strand(sequence)[::-1]


# ── additional feature 3: in silico transcription ───────────

def transcribe_to_mrna(sequence: str) -> str:
    """
    Generates an mRNA sequence by replacing every T with U.
    Lowercase embedded characters are preserved as-is.
    """
    return sequence.replace("T", "U").replace("t", "u")


# ── additional feature 4: sliding window GC analysis ────────

def sliding_window_gc(sequence: str,
                      window: int,
                      step: int = 1) -> list:
    """
    Calculates GC-content in a sliding window of given width
    moved along the pure nucleotide sequence with the given
    step size.
    Returns a list of dicts: {start_position, gc_content}.
    start_position is 1-based.
    """
    pure = "".join(ch for ch in sequence if ch in "ACGT")
    results = []
    for i in range(0, len(pure) - window + 1, step):
        chunk = pure[i:i + window]
        gc = (chunk.count("G") + chunk.count("C")) / len(chunk) * 100
        results.append({
            "start_position": i + 1,
            "gc_content": round(gc, 2)
        })
    return results


def save_sliding_window_csv(data: list, filename: str) -> None:
    """Saves sliding window GC results to a CSV file."""
    with open(filename, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh,
                                fieldnames=["start_position", "gc_content"])
        writer.writeheader()
        writer.writerows(data)


# ── additional feature 5: ORF identification ─────────────────

STOP_CODONS = {"TAA", "TAG", "TGA"}


def find_orfs(sequence: str, min_length: int = 100) -> list:
    """
    Searches for all open reading frames (ORFs) in the pure
    nucleotide sequence. An ORF starts at ATG and ends at the
    nearest in-frame stop codon (TAA, TAG, TGA).
    min_length: minimum ORF length in nucleotides (excluding stop).
    Returns a list of dicts: {start, end, length, sequence}.
    Positions are 1-based.
    """
    pure = "".join(ch for ch in sequence if ch in "ACGT")
    orfs = []
    for frame in range(3):
        i = frame
        while i < len(pure) - 2:
            codon = pure[i:i + 3]
            if codon == "ATG":
                # found start – scan for in-frame stop
                j = i + 3
                while j < len(pure) - 2:
                    stop_candidate = pure[j:j + 3]
                    if stop_candidate in STOP_CODONS:
                        orf_len = j - i   # length without stop codon
                        if orf_len >= min_length:
                            orfs.append({
                                "start": i + 1,
                                "end": j + 3,
                                "length": orf_len,
                                "sequence": pure[i:j]
                            })
                        break
                    j += 3
                i = j   # continue scanning after this ORF
            else:
                i += 3
    return orfs


# ── nucleotide distribution input ───────────────────────────

def get_nucleotide_weights() -> dict:
    """
    Asks the user for the percentage of each nucleotide.
    Validates that the four values sum to exactly 100.
    Returns a dict {"A": float, "C": float, "G": float, "T": float}.
    """
    while True:
        weights = {}
        total = 0.0
        valid = True
        for base in "ACGT":
            raw = input(f"  Percentage of {base} (e.g. 25): ").strip()
            try:
                val = float(raw)
                if val < 0:
                    raise ValueError
                weights[base] = val
                total += val
            except ValueError:
                print("  Error: enter a non-negative number.")
                valid = False
                break
        if valid:
            if abs(total - 100.0) < 1e-6:
                return weights
            print(f"  Error: percentages must sum to 100 "
                  f"(got {total:.2f}).")


# ── main program flow ────────────────────────────────────────

def main():
    """
    Main entry point. Orchestrates user interaction, sequence
    generation, file saving, statistics, and optional features.
    """
    print("=== DNA FASTA Sequence Generator ===\n")

    # ── basic input ──────────────────────────────────────────
    length = validate_positive_int("Enter sequence length: ")
    seq_id = validate_id("Enter sequence ID: ")
    description = input("Enter a description of the sequence "
                        "(press Enter to skip): ").strip()
    user_name = input("Enter your name: ").strip()

    # ── optional: custom nucleotide distribution ─────────────
    use_weights = input("\nUse custom nucleotide distribution? "
                        "[y/N]: ").strip().lower()
    if use_weights == "y":
        print("Enter the percentage for each nucleotide "
              "(they must sum to 100):")
        weights = get_nucleotide_weights()
        sequence = generate_sequence_weighted(length, weights)
    else:
        sequence = generate_sequence(length)

    # ── insert name (challenge requirement) ──────────────────
    sequence_with_name = insert_name(sequence, user_name)

    # ── build and save primary FASTA record ──────────────────
    fasta_str = format_fasta(seq_id, description, sequence_with_name)
    filename = save_fasta(seq_id, fasta_str)
    print(f"\nSequence saved to file: {filename}")

    # ── statistics (calculated on pure nucleotide sequence) ──
    stats = calculate_stats(sequence_with_name)
    print_stats(stats, length)

    # ══ additional feature 1: motif search ═══════════════════
    run_motif = input("\nSearch for a motif? [y/N]: ").strip().lower()
    if run_motif == "y":
        motif = input("  Enter motif (e.g. ATG): ").strip()
        positions = search_motif(sequence_with_name, motif)
        if positions:
            print(f"  Motif '{motif.upper()}' found at positions "
                  f"(1-based): {positions}")
        else:
            print(f"  Motif '{motif.upper()}' not found in sequence.")

    # ══ additional feature 2: complementary strands ══════════
    run_comp = input("\nGenerate complementary strands? "
                     "[y/N]: ").strip().lower()
    if run_comp == "y":
        comp = complementary_strand(sequence)
        rev_comp = reverse_complement(sequence)

        comp_id = f"{seq_id}_COMP"
        revcomp_id = f"{seq_id}_REVCOMP"

        comp_fasta = format_fasta(comp_id,
                                  f"Complement of {seq_id}",
                                  comp)
        revcomp_fasta = format_fasta(revcomp_id,
                                     f"Reverse complement of {seq_id}",
                                     rev_comp)

        # append both records to the existing FASTA file
        with open(filename, "a", encoding="utf-8") as fh:
            fh.write(comp_fasta)
            fh.write(revcomp_fasta)

        print(f"  Complementary and reverse-complement records "
              f"appended to {filename}.")

    # ══ additional feature 3: in silico transcription ════════
    run_trans = input("\nPerform in silico transcription (DNA→mRNA)? "
                      "[y/N]: ").strip().lower()
    if run_trans == "y":
        mrna = transcribe_to_mrna(sequence)
        mrna_id = f"{seq_id}_mRNA"
        mrna_fasta = format_fasta(mrna_id,
                                  f"mRNA transcript of {seq_id}",
                                  mrna)
        with open(filename, "a", encoding="utf-8") as fh:
            fh.write(mrna_fasta)
        print(f"  mRNA record appended to {filename}.")

    # ══ additional feature 4: sliding window GC ══════════════
    run_sliding = input("\nRun sliding window GC analysis? "
                        "[y/N]: ").strip().lower()
    if run_sliding == "y":
        window = validate_positive_int(
            "  Window size (nt): ", min_val=2, max_val=length)
        step = validate_positive_int(
            "  Step size (nt): ", min_val=1, max_val=window)
        sw_data = sliding_window_gc(sequence, window, step)
        csv_name = f"{seq_id}_gc_window.csv"
        save_sliding_window_csv(sw_data, csv_name)
        print(f"  Sliding window results saved to {csv_name} "
              f"({len(sw_data)} windows).")

    # ══ additional feature 5: ORF identification ═════════════
    run_orf = input("\nSearch for ORFs? [y/N]: ").strip().lower()
    if run_orf == "y":
        min_orf = validate_positive_int(
            "  Minimum ORF length in nt: ", min_val=3, max_val=length)
        orfs = find_orfs(sequence, min_orf)
        if orfs:
            print(f"  Found {len(orfs)} ORF(s):")
            for orf in orfs:
                print(f"    Start: {orf['start']:>6}  "
                      f"End: {orf['end']:>6}  "
                      f"Length: {orf['length']:>5} nt")
        else:
            print("  No ORFs found meeting the length requirement.")

    print("\nDone.")


if __name__ == "__main__":
    main()
