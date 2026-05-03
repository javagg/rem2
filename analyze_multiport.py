#!/usr/bin/env python3
"""
Analyze multi-port S-parameter CSV output to diagnose Z-matrix diagonalization issue.
"""
import csv
import sys
from pathlib import Path

def analyze_s_params(csv_path):
    """Read and analyze S-parameter CSV file."""
    print(f"\n=== Analyzing {csv_path} ===\n")
    
    if not Path(csv_path).exists():
        print(f"ERROR: File not found: {csv_path}")
        return
    
    # Read the CSV file
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    if not rows:
        print("ERROR: CSV file is empty")
        return
    
    # Get all columns
    cols = list(rows[0].keys())
    print(f"Columns: {cols}\n")
    
    # Try to identify frequency and S-parameter columns
    freq_col = next((c for c in cols if 'freq' in c.lower()), None)
    s_cols = [c for c in cols if c.startswith('S') and '[' in c]
    
    if not freq_col:
        print("WARNING: Could not find frequency column")
        freq_col = cols[0]
    
    print(f"Frequency column: {freq_col}")
    print(f"S-parameter columns: {s_cols}\n")
    
    if not s_cols:
        print("No S-parameter data found")
        return
    
    # Extract matrix dimension (e.g., S[1][1] -> 2x2, S[1][4] -> 4x4)
    max_idx = 0
    for col in s_cols:
        # Parse S[i][j]
        parts = col.replace('S[', '').replace(']', '').split('][')
        if len(parts) == 2:
            idx = max(int(parts[0]), int(parts[1]))
            max_idx = max(max_idx, idx)
    
    n_ports = max_idx + 1
    print(f"Detected {n_ports}-port network\n")
    
    # Print first frequency point's S-matrix
    if rows:
        freq_val = rows[0].get(freq_col, "N/A")
        print(f"Frequency: {freq_val}")
        print("\nS-matrix at first frequency point:")
        
        # Build matrix display
        matrix_data = {}
        for col in s_cols:
            parts = col.replace('S[', '').replace(']', '').split('][')
            if len(parts) == 2:
                i, j = int(parts[0]), int(parts[1])
                val = rows[0].get(col, "0")
                matrix_data[(i,j)] = val
        
        for i in range(n_ports):
            row_str = ""
            for j in range(n_ports):
                val = matrix_data.get((i,j), "0")
                row_str += f"{val:>15} "
            print(f"  S[{i}] = [{row_str}]")
    
    # Check if matrix is diagonal
    print("\n=== Checking Matrix Structure ===")
    is_diagonal = True
    off_diag_count = 0
    for col in s_cols:
        parts = col.replace('S[', '').replace(']', '').split('][')
        if len(parts) == 2:
            i, j = int(parts[0]), int(parts[1])
            if i != j:  # Off-diagonal
                val_str = rows[0].get(col, "0")
                try:
                    val = complex(val_str.replace('i', 'j') if 'i' in val_str else val_str)
                    if abs(val) > 1e-6:
                        is_diagonal = False
                        off_diag_count += 1
                        print(f"  S[{i}][{j}] = {val} (non-zero!)")
                except:
                    pass
    
    if is_diagonal:
        print(f"  ⚠️  Matrix is DIAGONAL: all {n_ports * (n_ports - 1)} off-diagonal elements are ~zero")
    else:
        print(f"  ✓ Matrix has {off_diag_count} non-zero off-diagonal elements")
    
    # Analyze diagonal values
    print("\n=== Diagonal Values ===")
    for i in range(n_ports):
        col = f"S[{i}][{i}]"
        if col in cols:
            val_str = rows[0][col]
            print(f"  S[{i}][{i}] = {val_str}")
    
    # Check Z-matrix if available
    z_cols = [c for c in cols if c.startswith('Z') and '[' in c]
    if z_cols:
        print("\n=== Z-matrix Data Found ===")
        for i in range(n_ports):
            row_str = ""
            for j in range(n_ports):
                col = f"Z[{i}][{j}]"
                if col in cols:
                    row_str += f"{rows[0][col]:>15} "
            if row_str:
                print(f"  Z[{i}] = [{row_str}]")

if __name__ == "__main__":
    # Look for S-parameter CSV files
    output_dir = Path("output")
    if not output_dir.exists():
        print(f"Output directory not found: {output_dir}")
        sys.exit(1)
    
    s_param_files = list(output_dir.glob("**/s_params.csv"))
    
    if not s_param_files:
        print("No s_params.csv files found in output/")
        print("Searching for any CSV files...")
        csv_files = list(output_dir.glob("**/*.csv"))
        for f in csv_files:
            print(f"  {f}")
        sys.exit(1)
    
    for csv_file in s_param_files:
        analyze_s_params(str(csv_file))
