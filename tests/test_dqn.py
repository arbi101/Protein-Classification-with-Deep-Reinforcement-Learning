"""
test_dqn.py
-----------
Training script for the Deep Q-Network (DQN) agent.

This trains the older move-based agent from agent_dqn.py. The current Django
website uses separately trained constructive weights through
agent_dqn_model.py, so this file remains for experiment reproducibility.

This script:
  1. Loads a set of realistic UniProt protein sequences
  2. Trains a SHARED DQN agent across all proteins (transfer learning)
  3. Saves the trained weights to 'dqn_weights.pth' after each protein
  4. Logs the best energy and training time per protein to 'dqn_results.txt'
ACTION SPACE DESIGN (Why these 4 moves?):
  The agent's action space consists of 4 local structural moves, each selected
  from the protein folding literature (Cebián et al., 2012) and proven effective
  for 2D HP lattice protein folding:
  
  1. END_FLIP    : Move a terminal residue to a free adjacent site of its anchor.
                   Efficient for extending the chain; low computational cost.
                   
  2. KINK_JUMP   : Flip an internal corner residue (kink) to the opposite corner.
                   Preserves backbone connectivity; enables local optimizations.
                   
  3. CRANKSHAFT  : Rotate two consecutive internal residues in a U-shaped motif.
                   High probability of acceptance in valid conformations; effective
                   for rearranging hydrophobic cores.
                   
  4. PIVOT       : Rotate a tail of residues ±90° or 180° around a pivot point.
                   Global move with lower acceptance but high impact when valid.
                   Helps escape local minima.
  
  This set balances LOCAL MOVES (end_flip, kink_jump, crankshaft) for fine-tuning
  with GLOBAL MOVES (pivot) for escaping local minima. The randomness in:
    - Move TYPE selection (which move to attempt)
    - Move PARAMETER selection (which residue, which direction/angle)
    - State TRANSITIONS (from experience replay)
  is essential for RL exploration—the agent must sample diverse conformations
  to learn effective folding policies.
Note on the shared agent approach:
  A single DQNAgent is reused across all proteins. The REPLAY BUFFER is
  cleared between proteins (because different proteins produce grids of
  varying effective content, which would cause inconsistency in the buffer).
  However, the NETWORK WEIGHTS are preserved across proteins, allowing
  the agent to transfer learned folding strategies from one sequence to another.
"""

import time
import statistics
import sys
from collections import deque

# Import shared utilities from the Q-Learning module
from test_ql import sequence_to_hp, calculate_energy, apply_move, MOVE_TYPES
# Import the function that loads the 20 benchmark UniProt protein sequences
from compare_hc_sa_multistep import fetch_realistic_proteins
# Import the historical move-based DQN and its 200x200 state encoder. The final
# website model is the separate constructive implementation in agent_dqn_model.
from agent_dqn import DQNAgent, positions_to_grid


def generate_2d_structure_dqn(hp_string, agent=None, episodes=1000,
                               max_steps_per_episode=100, batch_size=64,
                               target_update_freq=20, epsilon_start=1.0,
                               epsilon_end=0.05, epsilon_decay_episodes=None):
    """
    Trains the DQN agent on a single HP sequence and returns the best
    conformation found during training.

    At each episode, the agent starts from a straight-line conformation and
    applies moves selected by the ε-greedy policy. After each step, the
    transition is stored in the replay buffer and a training step is performed.

    Parameters
    ----------
    hp_string               : str       Binary HP sequence to fold
    agent                   : DQNAgent  Pre-created agent (shared across proteins)
                                        If None, a new agent is created
    episodes                : int       Number of training episodes (default: 1000)
    max_steps_per_episode   : int       Max moves per episode (default: 100)
    batch_size              : int       Mini-batch size for each gradient step (default: 64)
    target_update_freq      : int       Episodes between target network syncs (default: 20)
    epsilon_start           : float     Initial exploration rate (default: 1.0 = 100% random)
    epsilon_end             : float     Final exploration rate (default: 0.05 = 5% random)
    epsilon_decay_episodes  : int       Episodes over which ε decays (default: 60% of total)

    Returns
    -------
    best_positions : list of (int, int)  Best conformation found during training
    best_energy    : int                 Minimum energy found during training
    """
    n = len(hp_string)  # number of residues in the sequence

    # Reusing a provided agent preserves neural weights across proteins; passing
    # None starts an independent training experiment.
    if agent is None:
        # lr=5e-4: slightly lower learning rate for more stable training
        agent = DQNAgent(num_actions=len(MOVE_TYPES), lr=5e-4)

    # Set the ε decay horizon: by default, ε decays over the first 60% of episodes
    if epsilon_decay_episodes is None:
        epsilon_decay_episodes = int(episodes * 0.6)

    # Initialize best solution as the straight-line starting conformation
    best_positions = [(i, 0) for i in range(n)]
    best_energy = calculate_energy(best_positions, hp_string)

    losses = []  # track training losses for progress reporting

    # ── Training loop: one episode = one folding attempt from scratch ──────
    # Geometry resets every episode, while replay/network knowledge persists.
    for episode in range(episodes):
        # Reset conformation to straight line at the start of each episode
        current_positions = [(i, 0) for i in range(n)]
        current_energy = calculate_energy(current_positions, hp_string)
        # Encode the initial conformation as a 200×200 grid for the CNN
        current_state = positions_to_grid(current_positions, hp_string)

        # ── Linear ε decay based on episode number ─────────────────────────
        # ε decreases from epsilon_start to epsilon_end over the first
        # epsilon_decay_episodes episodes, then stays at epsilon_end
        # Clamp at epsilon_end so a small amount of exploration always remains.
        epsilon = max(
            epsilon_end,
            epsilon_start - (epsilon_start - epsilon_end) * episode / epsilon_decay_episodes
        )

        # ── Steps within this episode ──────────────────────────────────────
        # One inner iteration collects and learns from one transition.
        for step in range(max_steps_per_episode):
            # Select action using ε-greedy policy
            action_idx = agent.select_action(current_state, epsilon)
            action = MOVE_TYPES[action_idx]  # convert index to move name

            # Apply the selected move to the current conformation
            new_positions = apply_move(current_positions, action, hp_string)

            if new_positions is None:
                # INVALID move (collision / geometrically impossible)
                # Penalize the agent with a negative reward to discourage
                # invalid moves. The state remains unchanged.
                reward = -2.0
                next_state = current_state      # stay in the same state
                next_energy = current_energy
            else:
                # VALID move: compute the new energy
                next_energy = calculate_energy(new_positions, hp_string)
                # Base reward = energy improvement (positive = better fold)
                reward = float(current_energy - next_energy)

                # BONUS reward for discovering a new global minimum.
                # This strongly encourages the agent to explore lower-energy regions.
                if next_energy < best_energy:
                    reward += 5.0  # large bonus for beating the current best

                # Encode the new conformation as a grid
                next_state = positions_to_grid(new_positions, hp_string)

            # ── Store transition in the replay buffer ──────────────────────
            # Each transition (s, a, r, s') is saved for later training
            # Replay storage decouples data collection from mini-batch training.
            agent.memory.push(current_state, action_idx, reward, next_state)

            # ── Train the policy network on a random mini-batch ───────────
            loss = agent.train_step(batch_size)
            if loss > 0:
                losses.append(loss)  # only append if training actually ran

            # ── Transition to next state ───────────────────────────────────
            if new_positions is not None:
                # Move to the new state only if the move was valid
                current_positions = new_positions
                current_energy = next_energy
                current_state = next_state

                # Update global best if new conformation is the best seen so far
                if current_energy < best_energy:
                    best_energy = current_energy
                    best_positions = list(current_positions)

        # ── Target network synchronization ────────────────────────────────
        # Every target_update_freq episodes, copy policy weights to target net.
        # This stabilizes the Q-value targets used in training.
        # Keep the target fixed between periodic policy-to-target copies.
        if episode % target_update_freq == 0:
            agent.update_target_network()

        # ── Progress reporting ────────────────────────────────────────────
        # Print a summary every 10 episodes
        if (episode + 1) % 10 == 0:
            # Average loss over the last 10 training steps
            avg_loss = sum(losses[-10:]) / 10 if losses else 0
            print(f"  Episode {episode+1}/{episodes} | "
                  f"Best Energy: {best_energy} | "
                  f"Eps: {epsilon:.2f} | "
                  f"Avg Loss: {avg_loss:.4f}")

    # The reusable agent already stores learned weights, so return fold data only.
    return best_positions, best_energy


def main():
    """
    Main training and evaluation function.

    Loads all 20 benchmark proteins, trains the shared DQN agent on each,
    saves the weights after every protein, and logs results to dqn_results.txt.
    """
    import torch

    print("Loading realistic proteins for DQN Training...")
    # Load the 20 UniProt protein sequences used for benchmarking
    proteins = fetch_realistic_proteins()

    # Use all 20 proteins for training (train on all, generalize across all)
    test_proteins = proteins

    print(f"\nStarting DQN tests & Training...")

    # Open the output file for writing results
    txt_file = open("dqn_results.txt", "w")
    txt_file.write("=== DQN Optimization Results ===\n\n")

    # ── Shared agent across all proteins ──────────────────────────────────
    # A single agent is initialized once and used for ALL proteins.
    # This implements a form of CURRICULUM LEARNING / TRANSFER LEARNING:
    # patterns learned on one protein can be reused on the next.
    shared_agent = DQNAgent(num_actions=len(MOVE_TYPES), lr=5e-4)

    # Sequential processing retains network weights as transferable knowledge.
    for name, seq in test_proteins:
        # ── Clear the replay buffer between proteins ───────────────────────
        # This is necessary because different proteins produce different grid
        # content (different lengths → different spatial patterns), and mixing
        # them in the buffer without clearing can cause inconsistency.
        # IMPORTANT: only the BUFFER is cleared — the neural network WEIGHTS
        # (i.e., the "knowledge") are preserved across proteins!
        shared_agent.memory.buffer.clear()

        # Convert amino acid sequence to HP binary string
        hp_str = sequence_to_hp(seq)
        h_count = hp_str.count('H')  # number of H residues

        # Log protein info to both console and results file
        header = f"\nEvaluating {name} - Length: {len(seq)} AA - H count: {h_count}"
        print(header)
        txt_file.write(header + "\n")

        start = time.time()
        # Train the shared agent on this protein
        best_positions, best_energy = generate_2d_structure_dqn(
            hp_str,
            agent=shared_agent,
            episodes=50,                 # 50 episodes per protein (fast training)
            max_steps_per_episode=100,   # 100 moves per episode
            batch_size=64,               # sample 64 transitions per training step
            target_update_freq=20        # sync target network every 20 episodes
        )
        elapsed = time.time() - start

        # Log results
        res = f"  -> DQN Best Energy: {best_energy} | Time: {elapsed:.2f}s"
        print(res)

        # ── Save historical move-based weights to disk ────────────────────
        # Save after EVERY protein so an interrupted experiment still leaves
        # recoverable weights. The current Django site instead loads the
        # constructive checkpoint through agent_dqn_model.py.
        torch.save(shared_agent.policy_net.state_dict(), 'dqn_weights.pth')
        print(f"  [+] Weights updated and saved to tests/dqn_weights.pth after {name}\n")

    txt_file.close()
    print("\nDQN Tests completed! Results saved to dqn_results.txt")


# ── Script entry point ────────────────────────────────────────────────────────
if __name__ == '__main__':
    main()
