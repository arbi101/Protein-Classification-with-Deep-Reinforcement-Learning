"""Offline grid search for the tabular Q-Learning hyperparameters.

This experiment evaluates learning rate, discount factor and epsilon-decay
speed. It does not tune the constructive DQN used by the Django website.
"""

# time measures complete training plus exploitation for each run.
import time
# itertools.product builds the Cartesian hyperparameter grid.
import itertools

# Reuse the tabular learner and the common FASTA loader.
from test_ql import sequence_to_hp, generate_2d_structure_ql as run_ql, exploit_q_table
from compare_hc_sa_multistep import fetch_realistic_proteins

def main():
    """Evaluate every configured hyperparameter combination repeatedly."""

    print("Loading realistic proteins... picking the first one for Hyperparameter Tuning.")
    proteins = fetch_realistic_proteins()
    # Use one fixed protein so differences primarily reflect hyperparameters.
    name, seq = proteins[0]  # A1L190 (88 AA)
    
    hp_str = sequence_to_hp(seq)
    h_count = hp_str.count('H')
    
    print(f"\nProtein: {name} | Length: {len(seq)} AA | H count: {h_count}")
    
    # Alpha controls update size; gamma controls the value of delayed rewards.
    alphas = [0.1, 0.5, 0.9]
    gammas = [0.5, 0.9, 0.99]
    
    # Hold the computational budget constant across every grid combination.
    total_iters = 50000
    max_steps_per_episode = 100
    episodes = total_iters // max_steps_per_episode
    n_runs = 3 # 3 runs per combination to average out randomness
    
    # Fast decay over 20% of steps, Slow decay over 80% of steps
    eps_fast = int(total_iters * 0.2)
    eps_slow = int(total_iters * 0.8)
    epsilon_decays = [eps_fast, eps_slow]
    
    print(f"Testing {len(alphas) * len(gammas) * len(epsilon_decays)} combinations...")
    print(f"Fixed settings: {total_iters} total steps ({episodes} episodes x {max_steps_per_episode} steps)")
    
    # Detailed CSV rows support later ranking, plotting and statistical checks.
    csv_file = open("ql_hyperparams_data.csv", "w")
    csv_file.write("alpha,gamma,eps_decay,run_idx,energy,time_s\n")
    
    txt_file = open("ql_hyperparams_summary.txt", "w")
    txt_file.write("=== Q-Learning Hyperparameter Tuning ===\n")
    txt_file.write(f"Protein: {name} ({len(seq)} AA)\n")
    txt_file.write(f"Total Steps: {total_iters}\n\n")

    # Cartesian product generates every alpha/gamma/decay combination.
    for alpha, gamma, eps_decay in itertools.product(alphas, gammas, epsilon_decays):
        eps_label = "Fast (20%)" if eps_decay == eps_fast else "Slow (80%)"
        header = f"Testing Alpha: {alpha} | Gamma: {gamma} | Eps Decay: {eps_label}"
        print(header)
        txt_file.write(header + "\n")
        
        # Collect independent outcomes before computing aggregate quality.
        energies = []
        
        # Repeat each setting to reduce sensitivity to random exploration.
        for i in range(n_runs):
            start = time.time()
            _, train_energy, q_table = run_ql(
                hp_str,
                episodes=episodes,
                max_steps_per_episode=max_steps_per_episode,
                alpha=alpha,
                gamma=gamma,
                epsilon_decay=eps_decay
            )
            # Evaluate the learned table after the exploratory training phase.
            _, ex_energy = exploit_q_table(hp_str, q_table, steps=1000)
            # Keep the best fold seen during training or exploitation.
            final_energy = min(train_energy, ex_energy)
            
            elapsed = time.time() - start
            energies.append(final_energy)
            csv_file.write(f"{alpha},{gamma},{eps_label},{i+1},{final_energy},{elapsed:.4f}\n")
            
        # Lower/more-negative energy indicates the stronger fold.
        avg_en = sum(energies) / len(energies)
        min_en = min(energies)
        
        res = f"  -> Best Energy: {min_en} | Avg: {avg_en:.1f}"
        print(res)
        txt_file.write(res + "\n")
        
        # Persist progress because the complete grid search can take a long time.
        txt_file.flush()
        csv_file.flush()
        
    # Closing ensures operating-system buffers are committed to disk.
    csv_file.close()
    txt_file.close()
    print("\nTuning completed! Results saved to ql_hyperparams_data.csv and ql_hyperparams_summary.txt")

if __name__ == '__main__':
    # Keep the expensive grid search from running during module imports.
    main()
