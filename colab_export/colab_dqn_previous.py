"""
Deep Q-Learning (DQN) for 2D HP Protein Folding
[Constructive Approach with Relative Directions]

This script is designed to run seamlessly in Google Colab.
Unlike standard perturbative approaches, which start from a straight chain
and then bend it, this agent builds the protein one amino acid at a time.

At each step, the agent decides where to place the next amino acid using
relative directions (0: Forward, 1: Left, 2: Right). This technique makes the
state invariant to translation and rotation, reducing the curse of
dimensionality and allowing the neural network to generalize more effectively.
"""

import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import numpy as np
import random
import time
import csv
from collections import deque
import matplotlib.pyplot as plt
from IPython.display import clear_output

# ==========================================
# 1. ENVIRONMENT (Constructive HP Model)
# ==========================================

class ConstructiveHPEnv:
    def __init__(self, hp_string, grid_size=21):
        self.hp_string = hp_string
        self.n = len(hp_string)
        self.grid_size = grid_size
        self.half_grid = grid_size // 2
        self.reset()
        
    def reset(self):
        # The environment starts by placing the first two amino acids.
        # This defines an initial direction without loss of generality.
        self.positions = [(0, 0), (1, 0)]
        self.current_step = 2
        self.current_dir = (1, 0) # Facing East (+X)
        self.done = False
        self.total_energy = 0
        return self._get_state()
        
    def step(self, action):
        """
        Actions: 
        0 = Forward
        1 = Left
        2 = Right
        """
        dx, dy = self.current_dir
        
        if action == 0:
            new_dir = (dx, dy)
        elif action == 1:
            new_dir = (-dy, dx) # 90-degree counterclockwise rotation
        elif action == 2:
            new_dir = (dy, -dx) # 90-degree clockwise rotation
            
        new_pos = (self.positions[-1][0] + new_dir[0], self.positions[-1][1] + new_dir[1])
        
        # Self-intersection check (collision)
        if new_pos in self.positions:
            self.done = True
            # Strong collision penalty; the episode terminates early.
            return self._get_state(), -1.0, self.done 
            
        # Valid move: place the amino acid.
        self.positions.append(new_pos)
        self.current_dir = new_dir
        
        # Calculate newly formed H-H contacts.
        reward = 0.0
        if self.hp_string[self.current_step] == 'H':
            # Check against all previous amino acids except the bonded neighbour.
            for i, p in enumerate(self.positions[:-2]):
                if self.hp_string[i] == 'H':
                    dist = abs(p[0] - new_pos[0]) + abs(p[1] - new_pos[1])
                    if dist == 1:
                        reward += 1.0 # +1 for each new hydrophobic contact
                        self.total_energy -= 1 # Real HP energy is negative
                        
        self.current_step += 1
        if self.current_step == self.n:
            self.done = True
            # Optional bonus for completing the chain without collision.
            reward += 0.5 
            
        return self._get_state(), reward, self.done

    def _get_state(self):
        """
        State representation:
        Local 2D grid centered on the chain head and oriented according to the
        current direction. This makes the state invariant to global position
        and rotation.
        
        Channels:
        0: Positions of already placed H residues
        1: Positions of already placed P residues
        """
        state = np.zeros((2, self.grid_size, self.grid_size), dtype=np.float32)
        head_x, head_y = self.positions[-1]
        dx, dy = self.current_dir
        
        for i, (x, y) in enumerate(self.positions):
            # Coordinates relative to the chain head.
            rel_x = x - head_x
            rel_y = y - head_y
            
            # Rotate so the current direction always points upward (+Y).
            # Rotation formula based on the current direction vector (dx, dy).
            rot_x = rel_x * dy - rel_y * dx
            rot_y = rel_x * dx + rel_y * dy
            
            # Map the rotated coordinates onto the centered grid.
            grid_x = rot_x + self.half_grid
            grid_y = rot_y + self.half_grid
            
            if 0 <= grid_x < self.grid_size and 0 <= grid_y < self.grid_size:
                channel = 0 if self.hp_string[i] == 'H' else 1
                state[channel, grid_y, grid_x] = 1.0
                
        return state

# ==========================================
# 2. NEURAL NETWORK (CNN)
# ==========================================

class DQN_CNN(nn.Module):
    def __init__(self, grid_size=21):
        super(DQN_CNN, self).__init__()
        # Input shape: (Batch, 2, grid_size, grid_size)
        self.conv1 = nn.Conv2d(2, 16, kernel_size=3, stride=1, padding=1)
        self.conv2 = nn.Conv2d(16, 32, kernel_size=3, stride=1, padding=1)
        
        # Flattened feature size after the convolutional layers.
        self.flat_size = 32 * grid_size * grid_size
        
        self.fc1 = nn.Linear(self.flat_size, 128)
        self.fc2 = nn.Linear(128, 3) # 3 Actions: Forward, Left, Right
        
    def forward(self, x):
        x = F.relu(self.conv1(x))
        x = F.relu(self.conv2(x))
        x = x.view(x.size(0), -1) # Flatten
        x = F.relu(self.fc1(x))
        return self.fc2(x)

# ==========================================
# 3. REPLAY BUFFER
# ==========================================

class ReplayBuffer:
    def __init__(self, capacity=10000):
        self.buffer = deque(maxlen=capacity)
        
    def push(self, state, action, reward, next_state, done):
        self.buffer.append((state, action, reward, next_state, done))
        
    def sample(self, batch_size):
        batch = random.sample(self.buffer, batch_size)
        states, actions, rewards, next_states, dones = zip(*batch)
        return (np.array(states), np.array(actions), 
                np.array(rewards, dtype=np.float32), 
                np.array(next_states), np.array(dones, dtype=np.float32))
        
    def __len__(self):
        return len(self.buffer)

# ==========================================
# 4. TRAINING LOOP
# ==========================================

def prepare_training_items(hp_sequences):
    """
    Normalize training input into a list of (name, hp_string) pairs.
    Accepts a single HP string, a list of HP strings, or a list of
    (name, hp_string) tuples.
    """
    if isinstance(hp_sequences, str):
        return [("Sequence_1", hp_sequences)]

    training_items = []
    for idx, item in enumerate(hp_sequences, start=1):
        if isinstance(item, (tuple, list)) and len(item) >= 2:
            name, hp_string = item[0], item[1]
        else:
            name, hp_string = f"Sequence_{idx}", item
        training_items.append((str(name), str(hp_string)))

    return training_items


def init_sequence_stats(name, hp_string):
    """Create the statistics container used for one training protein."""
    return {
        "protein": name,
        "length": len(hp_string),
        "h_count": hp_string.count("H"),
        "episodes_seen": 0,
        "complete_episodes": 0,
        "best_energy": None,
        "best_reward": None,
        "best_length": 0,
        "last_energy": None,
        "last_reward": None,
        "last_length": 0,
        "sum_energy": 0.0,
        "sum_reward": 0.0,
        "sum_length": 0.0,
    }


def train_dqn(hp_sequences, episodes=1500, batch_size=64):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Training on device: {device}")
    
    training_items = prepare_training_items(hp_sequences)
        
    # Initialize policy and target networks.
    policy_net = DQN_CNN(grid_size=21).to(device)
    target_net = DQN_CNN(grid_size=21).to(device)
    target_net.load_state_dict(policy_net.state_dict())
    target_net.eval()
    
    optimizer = optim.Adam(policy_net.parameters(), lr=0.001)
    memory = ReplayBuffer(10000)
    
    # Reinforcement learning hyperparameters.
    gamma = 0.99
    epsilon_start = 1.0
    epsilon_end = 0.05
    epsilon_decay = 0.995
    epsilon = epsilon_start
    target_update = 10 # Update the target network every X episodes
    
    best_energy = {}
    sequence_stats = {}
    episode_logs = []
    history_energy = []
    training_config = {
        "episodes": episodes,
        "batch_size": batch_size,
        "grid_size": 21,
        "replay_buffer_capacity": 10000,
        "learning_rate": 0.001,
        "gamma": gamma,
        "epsilon_start": epsilon_start,
        "epsilon_end": epsilon_end,
        "epsilon_decay": epsilon_decay,
        "target_update_episodes": target_update,
        "actions": "Forward,Left,Right",
        "device": str(device),
        "num_training_proteins": len(training_items),
    }
    
    start_time = time.time()
    
    for episode in range(episodes):
        protein_name, hp_string = random.choice(training_items)
        env = ConstructiveHPEnv(hp_string, grid_size=21)
        state = env.reset()
        total_reward = 0
        
        while not env.done:
            # Epsilon-Greedy Action Selection
            if random.random() < epsilon:
                action = random.randrange(3)
            else:
                with torch.no_grad():
                    state_tensor = torch.FloatTensor(state).unsqueeze(0).to(device)
                    q_values = policy_net(state_tensor)
                    action = q_values.max(1)[1].item()
                    
            next_state, reward, done = env.step(action)
            memory.push(state, action, reward, next_state, done)
            
            state = next_state
            total_reward += reward
            
            # --- Training Step ---
            if len(memory) > batch_size:
                states, actions, rewards, next_states, dones = memory.sample(batch_size)
                
                states = torch.FloatTensor(states).to(device)
                actions = torch.LongTensor(actions).unsqueeze(1).to(device)
                rewards = torch.FloatTensor(rewards).unsqueeze(1).to(device)
                next_states = torch.FloatTensor(next_states).to(device)
                dones = torch.FloatTensor(dones).unsqueeze(1).to(device)
                
                # Compute Q(s, a).
                q_values = policy_net(states).gather(1, actions)
                
                # Compute target = R + gamma * max Q(s', a').
                with torch.no_grad():
                    next_q_values = target_net(next_states).max(1)[0].unsqueeze(1)
                    target = rewards + gamma * next_q_values * (1 - dones)
                    
                loss = F.mse_loss(q_values, target)
                
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                
        # Update epsilon and target network.
        epsilon = max(epsilon_end, epsilon * epsilon_decay)
        if episode % target_update == 0:
            target_net.load_state_dict(policy_net.state_dict())
            
        # Logging
        current_energy = env.total_energy
        reached_length = len(env.positions)
        completed_chain = reached_length == env.n
        history_energy.append(current_energy)
        
        if protein_name not in best_energy or current_energy < best_energy[protein_name]:
            best_energy[protein_name] = current_energy

        if protein_name not in sequence_stats:
            sequence_stats[protein_name] = init_sequence_stats(protein_name, hp_string)

        stats = sequence_stats[protein_name]
        stats["episodes_seen"] += 1
        stats["complete_episodes"] += int(completed_chain)
        stats["sum_energy"] += current_energy
        stats["sum_reward"] += total_reward
        stats["sum_length"] += reached_length
        stats["last_energy"] = current_energy
        stats["last_reward"] = total_reward
        stats["last_length"] = reached_length
        stats["best_length"] = max(stats["best_length"], reached_length)

        if stats["best_energy"] is None or current_energy < stats["best_energy"]:
            stats["best_energy"] = current_energy
        if stats["best_reward"] is None or total_reward > stats["best_reward"]:
            stats["best_reward"] = total_reward

        episode_logs.append({
            "episode": episode + 1,
            "protein": protein_name,
            "length": env.n,
            "h_count": hp_string.count("H"),
            "epsilon": round(epsilon, 6),
            "total_reward": round(total_reward, 6),
            "final_energy": current_energy,
            "reached_length": reached_length,
            "completed_chain": int(completed_chain),
            "best_energy_for_protein": best_energy[protein_name],
            "elapsed_s": round(time.time() - start_time, 3),
        })
            
        if (episode + 1) % 50 == 0:
            elapsed = time.time() - start_time
            print(f"Episode {episode+1:4d} | Epsilon: {epsilon:.2f} | "
                  f"Protein: {protein_name} | Seq Len: {env.n} | "
                  f"Best for Protein: {best_energy[protein_name]} | "
                  f"Reward: {total_reward:.2f} | Last Energy: {current_energy} | "
                  f"Len: {reached_length}/{env.n} | Time: {elapsed:.1f}s")
                  
    training_config["training_time_s"] = round(time.time() - start_time, 3)
    return policy_net, history_energy, best_energy, sequence_stats, episode_logs, training_config

# ==========================================
# 5. EXECUTION & VISUALIZATION
# ==========================================

def evaluate_dqn(policy_net, hp_string):
    """Evaluate the agent without exploration (greedy policy)."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    policy_net.eval()
    env = ConstructiveHPEnv(hp_string, grid_size=21)
    state = env.reset()
    
    while not env.done:
        with torch.no_grad():
            state_tensor = torch.FloatTensor(state).unsqueeze(0).to(device)
            q_values = policy_net(state_tensor)
            action = q_values.max(1)[1].item()
            
        state, reward, done = env.step(action)
        
    return env.total_energy, len(env.positions)


def evaluate_training_set(policy_net, training_items):
    """Run greedy evaluation for every protein after training."""
    results = {}
    for name, hp_string in training_items:
        energy, reached_length = evaluate_dqn(policy_net, hp_string)
        results[name] = {
            "greedy_final_energy": energy,
            "greedy_reached_length": reached_length,
            "greedy_completed_chain": int(reached_length == len(hp_string)),
        }
    return results


def write_dqn_reports(sequence_stats, episode_logs, greedy_results, training_config,
                      summary_path="dqn_previous_summary.csv",
                      episode_log_path="dqn_previous_episode_log.csv",
                      report_path="dqn_previous_report.txt"):
    """Write DQN results in CSV/TXT formats suitable for report tables."""
    episode_fields = [
        "episode", "protein", "length", "h_count", "epsilon", "total_reward",
        "final_energy", "reached_length", "completed_chain",
        "best_energy_for_protein", "elapsed_s",
    ]
    with open(episode_log_path, "w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=episode_fields)
        writer.writeheader()
        writer.writerows(episode_logs)

    summary_fields = [
        "protein", "length", "h_count", "episodes_seen", "best_energy",
        "avg_final_energy", "last_final_energy", "best_reward", "avg_reward",
        "last_reward", "best_length", "avg_reached_length", "completion_rate",
        "greedy_final_energy", "greedy_reached_length", "greedy_completed_chain",
        "training_time_s", "episodes", "batch_size", "gamma", "epsilon_start",
        "epsilon_end", "epsilon_decay", "target_update_episodes",
        "replay_buffer_capacity", "learning_rate",
    ]
    summary_rows = []

    for name, stats in sequence_stats.items():
        episodes_seen = max(1, stats["episodes_seen"])
        greedy = greedy_results.get(name, {})
        summary_rows.append({
            "protein": name,
            "length": stats["length"],
            "h_count": stats["h_count"],
            "episodes_seen": stats["episodes_seen"],
            "best_energy": stats["best_energy"],
            "avg_final_energy": round(stats["sum_energy"] / episodes_seen, 3),
            "last_final_energy": stats["last_energy"],
            "best_reward": round(stats["best_reward"], 3),
            "avg_reward": round(stats["sum_reward"] / episodes_seen, 3),
            "last_reward": round(stats["last_reward"], 3),
            "best_length": stats["best_length"],
            "avg_reached_length": round(stats["sum_length"] / episodes_seen, 3),
            "completion_rate": round(stats["complete_episodes"] / episodes_seen, 3),
            "greedy_final_energy": greedy.get("greedy_final_energy"),
            "greedy_reached_length": greedy.get("greedy_reached_length"),
            "greedy_completed_chain": greedy.get("greedy_completed_chain"),
            "training_time_s": training_config["training_time_s"],
            "episodes": training_config["episodes"],
            "batch_size": training_config["batch_size"],
            "gamma": training_config["gamma"],
            "epsilon_start": training_config["epsilon_start"],
            "epsilon_end": training_config["epsilon_end"],
            "epsilon_decay": training_config["epsilon_decay"],
            "target_update_episodes": training_config["target_update_episodes"],
            "replay_buffer_capacity": training_config["replay_buffer_capacity"],
            "learning_rate": training_config["learning_rate"],
        })

    with open(summary_path, "w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    best_values = [row["best_energy"] for row in summary_rows if row["best_energy"] is not None]
    avg_values = [row["avg_final_energy"] for row in summary_rows]
    completion_rates = [row["completion_rate"] for row in summary_rows]

    with open(report_path, "w") as txt_file:
        txt_file.write("=== Constructive DQN Training Report ===\n")
        txt_file.write("Training configuration:\n")
        for key, value in training_config.items():
            txt_file.write(f"  {key}: {value}\n")

        txt_file.write("\nAggregate metrics:\n")
        txt_file.write(f"  Proteins evaluated: {len(summary_rows)}\n")
        txt_file.write(f"  Median best energy: {float(np.median(best_values)):.3f}\n")
        txt_file.write(f"  Mean best energy: {float(np.mean(best_values)):.3f}\n")
        txt_file.write(f"  Best achieved energy: {min(best_values)}\n")
        txt_file.write(f"  Mean final energy: {float(np.mean(avg_values)):.3f}\n")
        txt_file.write(f"  Mean completion rate: {float(np.mean(completion_rates)):.3f}\n")

        txt_file.write("\nPer-protein summary:\n")
        txt_file.write("protein,length,h_count,episodes_seen,best_energy,avg_final_energy,"
                       "best_reward,avg_reward,completion_rate,greedy_final_energy,"
                       "greedy_reached_length\n")
        for row in summary_rows:
            txt_file.write(
                f"{row['protein']},{row['length']},{row['h_count']},"
                f"{row['episodes_seen']},{row['best_energy']},"
                f"{row['avg_final_energy']},{row['best_reward']},"
                f"{row['avg_reward']},{row['completion_rate']},"
                f"{row['greedy_final_energy']},{row['greedy_reached_length']}\n"
            )

    print(f"DQN episode log saved to '{episode_log_path}'")
    print(f"DQN summary CSV saved to '{summary_path}'")
    print(f"DQN text report saved to '{report_path}'")

PROTEINS = [
    ("A1L190", "MDDADPEERNYDNMLKMLSDLNKDLEKLLEEMEKISVQATWMAYDMVVMRTNPTLAESMRRLEDAFVNCKEEMEKNWQELLHETKQRL"),
    ("C9JLW8", "MTSSPVSRVVYNGKRTSSPRSPPSSSEIFTPAHEENVRFIYEAWQGVERDLRGQVPGGERGLVEEYVEKVPNPSLKTFKPIDLSDLKRRSTQDAKKS"),
    ("O00168", "MASLGHILVFCVGLLTMAKAESPKEHDPFTYDYQSLQIGGLVIAGILFILGILIVLSRRCRCKFNQQQRTGEPDEEEGTFRSSIRRLSTRRR"),
    ("O00453", "MLSRNDDICIYGGLGLGGLLLLAVVLLSACLCWLHRRVKRLERSWAQGSSEQELHYASLQRLPVPSSEGPDLRGRDKRGTKEDPRADYACIAENKPT"),
    ("O14949", "MGREFGNLTRMRHVISYSLSPFEQRAYPHVFTKGIPNVLRRIRESFFRVVPQFVVFYLIYTWGTEEFERSKRKNPAAYENDK"),
    ("O14957", "MVTRFLGPRYRELVKNWVPTAYTWGAVGAVGLVWATDWRLILDWVPYINGKFKKDN"),
    ("O15263", "MRVLYLLFSFLFIFLMPLPGVFGGIGDPVTCLKSGAICHPVFCPRRYKQIGTCGLPGTKCCKKP"),
    ("O43715", "MNSVGEACTDMKREYDQCFNRWFAEKFLKGDSSGDPCTDLFKRYQQCVQKAIKEKEIPIEGLEFMGHGKEKPENSS"),
    ("O75438", "MVNLLQIVRDHWVHVLVPMGFVIGCYLDRKSDERLTAFRNKSMLFKRELQPSEEVTWK"),
    ("O95777", "MTSALENYINRTVAVITSDGRMIVGTLKGFDQTINLILDESHERVFSSSQGVEQVVLGLYIVRGDNVAVIGEIDEETDSALDLGNIRAEPLNSVAH"),
    ("P04080", "MMCGAPSATQPATAETQHIADQVRSQLEEKENKKFPVFKAVSFKSQVVAGTNYFIKVHVGDEDFVHLRVFQSLPHENKPLTLSNYQTNKAKHDELTYF"),
    ("P04155", "MATMENKVICALVLVSMLALGTLAEAQTETCTVAPRERQNCGFPGVTPSQCANKGCCFDDTVRGVPWCFYPNTIDVPPEEECEF"),
    ("P04731", "MDPNCSCATGGSCTCTGSCKCKECKCTSCKKSCCSCCPMSCAKCAQGCICKGASEKCSCCA"),
    ("P05204", "MPKRKAEGDAKGDKAKVKDEPQRRSARLSAKPAPPKPEPKPKKAPAKKGEKVPKGKKGKADAGKEGNNPAENGDAKTDQAQKAEGAGDAK"),
    ("P06028", "MYGKIIFVLLLSEIVSISALSTTEVAMHTSTSSSVTKSYISSQTNGETGQLVHRFTVPAPVVIILIILCVMAGIIGTILLISYSIRRLIKA"),
    ("P07438", "MDPNCSCTTGGSCACAGSCKCKECKCTSCKKCCCSCCPVGCAKCAQGCVCKGSSEKCRCCA"),
    ("P10176", "MSVLTPLLLRGLTGSARRLPVPRAKIHSLPPEGKLGIMELAVGLTSCFVTFLLPAGWILSHLETYRRPE"),
    ("P29034", "MMCSSLEQALAVLVTTFHKYSCQEGDKFKLSKGEMKELLHKELPSFVGEKVDEEGLKKLMGSLDENSDQQVDFQEYAVFLALITVMCNDFFQGCPDRP"),
    ("P33763", "METPLEKALTTMVTTFHKYSGREGSKLTLSRKELKELIKKELCLGEMKESSIDDLMKSLDKNSDQEIDFKEYSVFLTMLCMAYNDFFLEDNK"),
    ("P35321", "MNSQQQKQPCTPPPQPQQQQVKQPCQPPPQEPCIPKTKEPCHPKVPEPCHPKVPEPCQPKVPEPCQPKVPEPCPSTVTPAPAQQKTKQK")
]

def sequence_to_hp(sequence):
    hydrophobic = set('ACFILMVWY')
    return ''.join('H' if aa in hydrophobic else 'P' for aa in sequence)

if __name__ == "__main__":
    training_proteins = [(name, sequence_to_hp(seq)) for name, seq in PROTEINS]
    
    print(f"Starting Constructive DQN Training on a dataset of {len(training_proteins)} proteins...")
    
    # Train the agent on multiple proteins.
    trained_net, energy_hist, best_E_dict, sequence_stats, episode_logs, training_config = train_dqn(
        training_proteins,
        episodes=20000,
        batch_size=128,
    )
    
    print("\nTraining Completed!")
    print("Best energies found during training for each protein:")
    for protein_name, energy in best_E_dict.items():
        stats = sequence_stats[protein_name]
        avg_energy = stats["sum_energy"] / stats["episodes_seen"]
        avg_reward = stats["sum_reward"] / stats["episodes_seen"]
        completion_rate = stats["complete_episodes"] / stats["episodes_seen"]
        print(
            f"  {protein_name:6s} | Len: {stats['length']:3d} | H: {stats['h_count']:3d} | "
            f"Episodes: {stats['episodes_seen']:4d} | Best E: {energy:4d} | "
            f"Avg E: {avg_energy:7.3f} | Avg Reward: {avg_reward:7.3f} | "
            f"Completion: {completion_rate:.3f}"
        )
        
    torch.save(trained_net.state_dict(), 'dqn_previous_weights.pth')
    print("Model weights saved to 'dqn_previous_weights.pth'")

    greedy_results = evaluate_training_set(trained_net, training_proteins)
    write_dqn_reports(sequence_stats, episode_logs, greedy_results, training_config)
        
    print("\n--- Evaluation on Unseen/Large Protein ---")
    test_hp = "HPHHHPPHPPPPHPHPHPPPPPHPPHPHPHPPPPPHPPPPPPHHHHPPPHPPPPPHPPHPHPPPPPHPHHHPHPPP"
    print(f"Testing Constructive DQN on sequence length {len(test_hp)}...")
    
    test_energy, test_len = evaluate_dqn(trained_net, test_hp)
    print(f"Evaluation Result - Energy: {test_energy} | Length reached: {test_len}/{len(test_hp)}")
    
    # Plot the learning curve. Energy varies across proteins, so a moving
    # average is plotted to show the overall trend.
    plt.figure(figsize=(10, 5))
    plt.plot(energy_hist, alpha=0.3, color='blue', label='Energy per episode (Mixed lengths)')
    
    # Moving-average smoothing.
    window = 50
    if len(energy_hist) > window:
        smoothed = np.convolve(energy_hist, np.ones(window)/window, mode='valid')
        plt.plot(range(window-1, len(energy_hist)), smoothed, color='red', linewidth=2, label='Moving Average')
        
    plt.title('DQN Constructive Approach Learning Curve (Multi-Protein Dataset)')
    plt.xlabel('Episodes')
    plt.ylabel('Energy (Lower is better)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('dqn_learning_curve.png')
    print("Learning curve saved to 'dqn_learning_curve.png'")
