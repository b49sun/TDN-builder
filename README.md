TDN Builder
===========

A Python pipeline for building 3D atomic-level PDB models of tetrahedral DNA nanostructures (TDN).

Dependencies
------------
- Python >= 3.8
- NumPy
- Biopython

Usage
-----
    python TDN_builder.py <input.txt>           # Standard mode
    python TDN_builder.py -autofill <input.txt> # Autofill mode

Input — Standard Mode
---------------------
Line 1: len_hinge (integer 1-4)
Lines 2-5: four DNA sequences A, B, C, D

Example:
    len_hinge = 3
    TCAACTGCCTGGTGATAGGGGGTTTACGACACTACGTGGGAAGGGGGTTTTACTATGGCGGCTCTTCTTTTTTTT
    ACATTCCTAAGTCTGAATTTTTTTTATTACAGCTTGCTACACGGGGGTTTCCCCCTATCACCAGGCAGTTGAAAA
    AAAAATTCAGACTTAGGAATGTTTTCCCCCTTCCCACGTAGTGTCGTAAACCCCCGTGTAGCAAGCTGTAATAAA
    AAAAAGAAGAGCCGCCATAGTAAAAACATTCCTAAGTCTGAATTTTTAAAAAAAATTCAGACTTAGGAATGTAAA

Input — Autofill Mode
---------------------
Line 1: len_hinge = N
Line 2: full A chain
Line 3: B[0] + B[1]
Line 4: C[0]

Example:
    len_hinge = 3
    TCAACTGCCTGGTGATAGGGGGTTTACGACACTACGTGGGAAGGGGGTTTTACTATGGCGGCTCTTCTTTTTTTT
    ACATTCCTAAGTCTGAATTTTTTTTATTACAGCTTGCTACACGGGGGTTT         
    AAAAATTCAGACTTAGGAATGTTTT

Missing segments are auto-completed via reverse complement.

Output
------
- tetrahedron_result.pdb  -- final structure
- tetrahedron_tmp.pdb     -- intermediate (auto-deleted on success)

Terminal prints total_loss, total time, and output path.

Notes
-----
- A-chain length must be divisible by 3.
- All sequences must be equal length in standard mode.
- 5'-terminal phosphate atoms (P, OP1, OP2) on the first residue of each
  chain are automatically removed for GROMACS compatibility.
