import time
import statistics
import re

from test_remc import sequence_to_hp, generate_2d_structure_remc as run_remc
from compare_hc_sa_multistep import fetch_realistic_proteins

def main():
    print(f"Loading 20 realistic protein sequences (length 50-100 aa) from 20_proteins.fasta...")
    proteins = fetch_realistic_proteins()
    
    # We will test convergence across different iteration levels
    iterations_list = [1000, 10000, 50000, 100000]
    n_runs = 3 # Reduced to 3 because REMC is computation heavy.
    
    print(f"\nStarting tests across iterations: {iterations_list} ({n_runs} runs per level)")
    
    csv_file = open("remc_convergence_data.csv", "w")
    csv_file.write("protein,length,iterations,algo,run_idx,energy,time_s\n")
    
    txt_file = open("remc_convergence_summary.txt", "w")
    txt_file.write(f"=== Convergence Comparison: REMC ===\n")
    txt_file.write(f"Testing {len(proteins)} proteins over iterations {iterations_list}\n\n")

    for name, seq in proteins:
        hp_str = sequence_to_hp(seq)
        h_count = hp_str.count('H')
        
        header = f"\nEvaluating {name} - Length: {len(seq)} AA - H count: {h_count}"
        print(header)
        txt_file.write(header + "\n")
        
        for iters in iterations_list:
            remc_energies = []
            remc_times = []
            
            for i in range(n_runs):
                start = time.time()
                # Run REMC (5 replicas)
                _, remc_en = run_remc(hp_str, iterations=iters, num_replicas=5)
                remc_t = time.time() - start
                remc_times.append(remc_t)
                remc_energies.append(remc_en)
                csv_file.write(f"{name},{len(seq)},{iters},REMC,{i+1},{remc_en},{remc_t:.4f}\n")
                
            remc_avg_en = statistics.mean(remc_energies)
            remc_min_en = min(remc_energies)
            remc_avg_t = statistics.mean(remc_times)
            
            res_remc = f"  [Iterations: {iters}] -> Best Energy: {remc_min_en:<4} | Avg: {remc_avg_en:<6.1f} | Avg Time: {remc_avg_t:.3f} s"
            print(res_remc)
            txt_file.write(res_remc + "\n")
            
            txt_file.flush()
            csv_file.flush()

    csv_file.close()
    txt_file.close()
    
    print("\nTests completed! Results saved to remc_convergence_data.csv and remc_convergence_summary.txt")

if __name__ == '__main__':
    main()
