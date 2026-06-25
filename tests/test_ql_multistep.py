"""Measure tabular Q-Learning convergence at increasing training budgets.

The website's ``ql`` option runs constructive DQN inference instead of this
offline Q-table experiment.
"""

# time measures complete train-and-exploit runtime for each repetition.
import time
# statistics.mean summarizes stochastic outcomes at each budget.
import statistics

# Import the tabular learner and shared benchmark FASTA loader.
from test_ql import sequence_to_hp, generate_2d_structure_ql as run_ql, exploit_q_table
from compare_hc_sa_multistep import fetch_realistic_proteins

def main():
    """Benchmark Q-Learning quality and runtime over multiple total steps."""

    print(f"Loading 20 realistic protein sequences (length 50-100 aa) from 20_proteins.fasta...")
    proteins = fetch_realistic_proteins()
    
    # We will test convergence across different total iterations
    # iterations_list represents total steps = episodes * max_steps_per_episode
    # These are total state transitions, not the number of episodes.
    iterations_list = [1000, 10000, 50000, 100000]
    max_steps_per_episode = 100
    n_runs = 3 # Keep at 3 for speed, comparable to REMC
    
    print(f"\nStarting Q-Learning tests across total iterations: {iterations_list} ({n_runs} runs per level)")
    
    # CSV supports plotting; TXT provides an immediately readable summary.
    csv_file = open("ql_convergence_data.csv", "w")
    csv_file.write("protein,length,iterations,algo,run_idx,energy,time_s\n")
    
    txt_file = open("ql_convergence_summary.txt", "w")
    txt_file.write(f"=== Convergence Comparison: Tabular Q-Learning ===\n")
    txt_file.write(f"Testing {len(proteins)} proteins over iterations {iterations_list}\n\n")

    # Evaluate the identical budget schedule independently for every protein.
    for name, seq in proteins:
        hp_str = sequence_to_hp(seq)
        h_count = hp_str.count('H')
        
        header = f"\nEvaluating {name} - Length: {len(seq)} AA - H count: {h_count}"
        print(header)
        txt_file.write(header + "\n")
        
        for total_iters in iterations_list:
            # Convert the transition budget into complete training episodes.
            episodes = total_iters // max_steps_per_episode
            
            # Store all repetitions before calculating consistency statistics.
            ql_energies = []
            ql_times = []
            
            # Average several stochastic runs at every budget level.
            for i in range(n_runs):
                start = time.time()
                # Run Q-Learning
                _, train_energy, q_table = run_ql(
                    hp_str, 
                    episodes=episodes, 
                    max_steps_per_episode=max_steps_per_episode,
                    alpha=0.1,
                    gamma=0.95
                )
                
                # Evaluate the table separately after exploratory training.
                _, ex_energy = exploit_q_table(hp_str, q_table, steps=1000)
                final_energy = min(train_energy, ex_energy)
                
                ql_t = time.time() - start
                ql_times.append(ql_t)
                ql_energies.append(final_energy)
                csv_file.write(f"{name},{len(seq)},{total_iters},QL,{i+1},{final_energy},{ql_t:.4f}\n")
                
            # Best energy captures peak quality; means capture expected behavior.
            ql_avg_en = statistics.mean(ql_energies)
            ql_min_en = min(ql_energies)
            ql_avg_t = statistics.mean(ql_times)
            
            res_ql = f"  [Iterations: {total_iters}] -> Best Energy: {ql_min_en:<4} | Avg: {ql_avg_en:<6.1f} | Avg Time: {ql_avg_t:.3f} s"
            print(res_ql)
            txt_file.write(res_ql + "\n")
            
            # Make partial results visible and resilient to interrupted runs.
            txt_file.flush()
            csv_file.flush()

    # Release report files after all proteins and budgets are complete.
    csv_file.close()
    txt_file.close()
    
    print("\nTests completed! Results saved to ql_convergence_data.csv and ql_convergence_summary.txt")

if __name__ == '__main__':
    # Execute only when launched directly, not when imported as a helper.
    main()
