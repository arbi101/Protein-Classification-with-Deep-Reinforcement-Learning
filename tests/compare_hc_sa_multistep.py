"""Benchmark classical HP-folding algorithms at multiple iteration budgets.

The script loads up to 20 FASTA records, repeatedly runs HC, SA and MC, and
writes both machine-readable CSV data and a human-readable summary. It is an
offline experiment script; Django imports the algorithms, not this main loop.
"""

# time measures wall-clock execution for each independent run.
import time
# statistics.mean summarizes stochastic energy and runtime distributions.
import statistics
# re extracts a compact accession from pipe-delimited FASTA headers.
import re

# Import the shared HP conversion and expose each implementation under a short,
# consistent name for the benchmark loop below.
from test_hc import sequence_to_hp, generate_2d_structure as run_hc
from test_sa import generate_2d_structure_sa as run_sa
from test_mc import generate_2d_structure_mc as run_mc

def fetch_realistic_proteins(file_path="20_proteins.fasta"):
    """Parse up to 20 FASTA records into ``(name, sequence)`` tuples."""

    # The default path assumes this script is launched from the tests folder.
    with open(file_path, 'r') as f:
        fasta_data = f.read()
    
    # These variables implement a small state machine for multi-record FASTA.
    proteins = []
    current_name = ""
    current_seq = ""
    for line in fasta_data.splitlines():
        if line.startswith(">"):
            # A new header closes and stores the previous FASTA record.
            if current_name and current_seq:
                proteins.append((current_name, current_seq))
            # Extract simple name (e.g., QCR10_HUMAN)
            match = re.search(r'\|([^|]+)\|', line)
            current_name = match.group(1) if match else line[1:10]
            current_seq = ""
        else:
            # FASTA sequences may span multiple physical lines.
            current_seq += line.strip()
    # Save the final record because no later header will close it.
    if current_name and current_seq:
        proteins.append((current_name, current_seq))
    
    # The benchmark dataset is intentionally capped for predictable runtime.
    return proteins[:20]

def main():
    """Run the complete convergence experiment and write its reports."""

    print(f"Loading 20 realistic protein sequences (length 50-100 aa) from 20_proteins.fasta...")
    proteins = fetch_realistic_proteins()
    
    # We will test convergence across different iteration levels
    # Every algorithm receives the same proposal budget at each level.
    iterations_list = [1000, 10000, 50000, 100000]
    n_runs = 25  # Increased from 10 to 25 for better statistical significance
    
    print(f"\nStarting tests across iterations: {iterations_list} ({n_runs} runs per algorithm/level)")
    
    # Save a detailed CSV and a readable text summary
    # Opening with "w" replaces previous results, ensuring one coherent run.
    csv_file = open("convergence_data.csv", "w")
    csv_file.write("protein,length,iterations,algo,run_idx,energy,time_s\n")
    
    txt_file = open("convergence_summary.txt", "w")
    txt_file.write(f"=== Convergence Comparison: HC vs SA ===\n")
    txt_file.write(f"Testing {len(proteins)} proteins over iterations {iterations_list}\n\n")

    # Each protein is reduced to the same HP representation before comparison.
    for name, seq in proteins:
        hp_str = sequence_to_hp(seq)
        h_count = hp_str.count('H')
        
        header = f"\nEvaluating {name} - Length: {len(seq)} AA - H count: {h_count}"
        print(header)
        txt_file.write(header + "\n")
        
        # Increasing budgets reveal how solution quality and runtime converge.
        for iters in iterations_list:
            print(f"  [Iterations: {iters}]")
            txt_file.write(f"  [Iterations: {iters}]\n")
            
            # Keep per-run samples so both best-case and mean behavior can be
            # reported after all stochastic repetitions finish.
            hc_energies = []
            hc_times = []
            sa_energies = []
            sa_times = []
            mc_energies = []
            mc_times = []
            
            # Repeated runs reduce the influence of random luck.
            for i in range(n_runs):
                # Run HC
                start = time.time()
                _, hc_en = run_hc(hp_str, iterations=iters)
                hc_t = time.time() - start
                hc_times.append(hc_t)
                hc_energies.append(hc_en)
                csv_file.write(f"{name},{len(seq)},{iters},HC,{i+1},{hc_en},{hc_t:.4f}\n")
                
                # Run SA
                start = time.time()
                # Aggressive cooling schedule: initial_t=30.0, final_t=0.001
                # This gives SA more exploration time to escape local minima
                _, sa_en = run_sa(hp_str, iterations=iters, initial_t=30.0, final_t=0.001)
                sa_t = time.time() - start
                sa_times.append(sa_t)
                sa_energies.append(sa_en)
                csv_file.write(f"{name},{len(seq)},{iters},SA,{i+1},{sa_en},{sa_t:.4f}\n")
                
                # Run Basic MC
                start = time.time()
                _, mc_en = run_mc(hp_str, iterations=iters, temperature=5.0)
                mc_t = time.time() - start
                mc_times.append(mc_t)
                mc_energies.append(mc_en)
                csv_file.write(f"{name},{len(seq)},{iters},MC,{i+1},{mc_en},{mc_t:.4f}\n")
                
            # Report both peak quality and average behavior across runs.
            hc_avg_en = statistics.mean(hc_energies)
            hc_min_en = min(hc_energies)
            hc_avg_t = statistics.mean(hc_times)
            
            sa_avg_en = statistics.mean(sa_energies)
            sa_min_en = min(sa_energies)
            sa_avg_t = statistics.mean(sa_times)
            
            mc_avg_en = statistics.mean(mc_energies)
            mc_min_en = min(mc_energies)
            mc_avg_t = statistics.mean(mc_times)
            
            res_hc = f"    HC -> Best Energy: {hc_min_en:<4} | Avg: {hc_avg_en:<6.1f} | Avg Time: {hc_avg_t:.3f} s"
            res_sa = f"    SA -> Best Energy: {sa_min_en:<4} | Avg: {sa_avg_en:<6.1f} | Avg Time: {sa_avg_t:.3f} s"
            res_mc = f"    MC -> Best Energy: {mc_min_en:<4} | Avg: {mc_avg_en:<6.1f} | Avg Time: {mc_avg_t:.3f} s"
            
            print(res_hc)
            print(res_sa)
            print(res_mc)
            txt_file.write(res_hc + "\n")
            txt_file.write(res_sa + "\n")
            txt_file.write(res_mc + "\n")
            
            # Print a flush so we can see progress in the file right away
            # Flush incremental progress in case a long benchmark is interrupted.
            txt_file.flush()
            csv_file.flush()

    # Explicitly close both output streams once all buffered data is written.
    csv_file.close()
    txt_file.close()
    
    print("\nTests completed! Results saved to convergence_data.csv and convergence_summary.txt")

if __name__ == '__main__':
    # Do not launch the benchmark when importing fetch_realistic_proteins().
    main()
