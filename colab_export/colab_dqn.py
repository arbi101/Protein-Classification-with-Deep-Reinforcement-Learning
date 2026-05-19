"""
colab_dqn.py
-----------
DQN agent implementation for 2D HP protein folding (standalone executable).

This module provides:
  - DQNCNN: Convolutional neural network for Q-function approximation
  - ReplayBuffer: Experience replay for stable training
  - DQNAgent: Main agent class managing policy/target networks and training
  - Core folding utilities: move application, energy calculation, grid encoding

ACTION SPACE DESIGN (Why these 4 moves?):
  The agent's action space consists of 4 established structural moves from
  protein folding literature, proven effective for 2D HP lattice models:
  
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
    - ε-greedy exploration (explore vs exploit via policy network)
    - Experience replay sampling (random mini-batch training)
  is essential for DQN exploration—the agent must sample diverse conformations
  and trajectories to learn effective folding policies and avoid local optima.
"""

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import random
import time
from collections import deque

# ─────────────────────────────────────────────────────────────────────────────
# CORE LOGIC (Proteins & Moves)
# ─────────────────────────────────────────────────────────────────────────────

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

MOVE_TYPES = ['end_flip', 'kink_jump', 'crankshaft', 'pivot']

def sequence_to_hp(sequence):
    hydrophobic = set('ACFILMVWY')
    return ''.join('H' if aa in hydrophobic else 'P' for aa in sequence)

def get_neighbors(x, y):
    return [(x+1, y), (x-1, y), (x, y+1), (x, y-1)]

def calculate_energy(pos_list, hp_string):
    e = 0
    pos_dict = {pos: i for i, pos in enumerate(pos_list)}
    for i, (x, y) in enumerate(pos_list):
        if hp_string[i] == 'H':
            for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                nb = (x + dx, y + dy)
                if nb in pos_dict:
                    j = pos_dict[nb]
                    if abs(i - j) > 1 and hp_string[j] == 'H':
                        e += 1
    return -(e // 2)

def apply_move(current_positions, move_type, hp_string):
    n = len(hp_string)
    new_positions = list(current_positions)

    if move_type == 'end_flip':
        idx = random.choice([0, n - 1])
        anchor_idx = 1 if idx == 0 else n - 2
        ax, ay = current_positions[anchor_idx]
        candidates = []
        pos_set = set(current_positions)
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            nx, ny = ax + dx, ay + dy
            if (nx, ny) not in pos_set and (nx, ny) != current_positions[idx]:
                candidates.append((nx, ny))
        if not candidates: return None
        new_positions[idx] = random.choice(candidates)
        return new_positions

    elif move_type == 'kink_jump':
        if n <= 2: return None
        idx = random.randint(1, n - 2)
        p_prev = current_positions[idx - 1]
        p_curr = current_positions[idx]
        p_next = current_positions[idx + 1]
        if abs(p_prev[0] - p_next[0]) == 1 and abs(p_prev[1] - p_next[1]) == 1:
            nx = p_prev[0] + p_next[0] - p_curr[0]
            ny = p_prev[1] + p_next[1] - p_curr[1]
            if (nx, ny) not in set(current_positions):
                new_positions[idx] = (nx, ny)
                return new_positions
        return None

    elif move_type == 'crankshaft':
        if n <= 3: return None
        idx = random.randint(1, n - 3)
        p_prev  = current_positions[idx - 1]
        p_curr  = current_positions[idx]
        p_next  = current_positions[idx + 1]
        p_next2 = current_positions[idx + 2]
        dist_sq = (p_prev[0] - p_next2[0])**2 + (p_prev[1] - p_next2[1])**2
        if dist_sq == 1:
            nx_curr = p_prev[0] + p_next2[0] - p_next[0]
            ny_curr = p_prev[1] + p_next2[1] - p_next[1]
            nx_next = p_prev[0] + p_next2[0] - p_curr[0]
            ny_next = p_prev[1] + p_next2[1] - p_curr[1]
            pos_set = set(current_positions)
            if (nx_curr, ny_curr) not in pos_set and (nx_next, ny_next) not in pos_set:
                new_positions[idx]     = (nx_curr, ny_curr)
                new_positions[idx + 1] = (nx_next, ny_next)
                return new_positions
        return None

    elif move_type == 'pivot':
        if n <= 2: return None
        pivot_idx = random.randint(1, n - 2)
        angle = random.choice([90, -90, 180])
        cx, cy = current_positions[pivot_idx]
        if angle == 90: cos_a, sin_a = 0, 1
        elif angle == -90: cos_a, sin_a = 0, -1
        else: cos_a, sin_a = -1, 0
        for i in range(pivot_idx + 1, n):
            x, y = current_positions[i]
            tx, ty = x - cx, y - cy
            new_positions[i] = (tx * cos_a - ty * sin_a + cx, tx * sin_a + ty * cos_a + cy)
        if len(set(new_positions)) == len(new_positions):
            return new_positions
        return None

    return None

# ─────────────────────────────────────────────────────────────────────────────
# NEURAL NETWORK ARCHITECTURE
# ─────────────────────────────────────────────────────────────────────────────

class DQNCNN(nn.Module):
    def __init__(self, num_actions=4):
        super(DQNCNN, self).__init__()
        self.conv_layers = nn.Sequential(
            nn.Conv2d(1, 16, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(2), # 200 -> 100
            nn.Conv2d(16, 32, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(2), # 100 -> 50
            nn.Conv2d(32, 64, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(2)  # 50 -> 25
        )
        self.fc_layers = nn.Sequential(
            nn.Linear(64 * 25 * 25, 256),
            nn.ReLU(),
            nn.Linear(256, num_actions)
        )
        
    def forward(self, x):
        x = self.conv_layers(x)
        x = x.view(x.size(0), -1)
        x = self.fc_layers(x)
        return x

class ReplayBuffer:
    def __init__(self, capacity):
        self.buffer = deque(maxlen=capacity)
    def push(self, state, action, reward, next_state):
        self.buffer.append((state, action, reward, next_state))
    def sample(self, batch_size):
        batch = random.sample(self.buffer, batch_size)
        states, actions, rewards, next_states = zip(*batch)
        return (np.stack(states), np.array(actions), np.array(rewards, dtype=np.float32), np.stack(next_states))
    def __len__(self):
        return len(self.buffer)

def positions_to_grid(positions, hp_string):
    # Griglia fissa a 200x200 per conservare la memoria spaziale!
    grid_size = 200
    grid = np.zeros((1, grid_size, grid_size), dtype=np.float32)
    if not positions: return grid
    ox, oy = positions[0]
    for i, (x, y) in enumerate(positions):
        tx = max(0, min(grid_size - 1, x - ox + grid_size//2))
        ty = max(0, min(grid_size - 1, y - oy + grid_size//2))
        grid[0, tx, ty] = 1.0 if hp_string[i] == 'H' else -1.0
    return grid

class DQNAgent:
    def __init__(self, num_actions=4, lr=1e-3, gamma=0.95, buffer_size=10000):
        # Auto-detect GPU for Colab
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        print(f"Using device: {self.device}")
        
        self.num_actions = num_actions
        self.gamma = gamma
        self.policy_net = DQNCNN(num_actions).to(self.device)
        self.target_net = DQNCNN(num_actions).to(self.device)
        self.target_net.load_state_dict(self.policy_net.state_dict())
        self.target_net.eval()
        self.optimizer = optim.Adam(self.policy_net.parameters(), lr=lr)
        self.memory = ReplayBuffer(buffer_size)
        
    def select_action(self, state_grid, epsilon):
        if random.random() < epsilon: return random.randrange(self.num_actions)
        with torch.no_grad():
            state_tensor = torch.FloatTensor(state_grid).unsqueeze(0).to(self.device)
            return self.policy_net(state_tensor).argmax().item()
            
    def train_step(self, batch_size):
        if len(self.memory) < batch_size: return 0.0
        states, actions, rewards, next_states = self.memory.sample(batch_size)
        states = torch.FloatTensor(states).to(self.device)
        actions = torch.LongTensor(actions).unsqueeze(1).to(self.device)
        rewards = torch.FloatTensor(rewards).unsqueeze(1).to(self.device)
        next_states = torch.FloatTensor(next_states).to(self.device)
        
        q_values = self.policy_net(states).gather(1, actions)
        with torch.no_grad():
            next_q_values = self.target_net(next_states).max(1)[0].unsqueeze(1)
        target_q_values = rewards + self.gamma * next_q_values
        loss = nn.MSELoss()(q_values, target_q_values)
        self.optimizer.zero_grad()
        loss.backward()
        self.optimizer.step()
        return loss.item()
        
    def update_target_network(self):
        self.target_net.load_state_dict(self.policy_net.state_dict())

# ─────────────────────────────────────────────────────────────────────────────
# TRAINING LOOP
# ─────────────────────────────────────────────────────────────────────────────

def generate_2d_structure_dqn(hp_string, agent, episodes=100, max_steps_per_episode=150,
                              batch_size=64, target_update_freq=20,
                              epsilon_start=1.0, epsilon_end=0.05, epsilon_decay_episodes=None):
    n = len(hp_string)
    if epsilon_decay_episodes is None:
        epsilon_decay_episodes = int(episodes * 0.7)
        
    best_positions = [(i, 0) for i in range(n)]
    best_energy = calculate_energy(best_positions, hp_string)
    losses = []
    
    for episode in range(episodes):
        current_positions = [(i, 0) for i in range(n)]
        current_energy = calculate_energy(current_positions, hp_string)
        current_state = positions_to_grid(current_positions, hp_string)
        epsilon = max(epsilon_end, epsilon_start - (epsilon_start - epsilon_end) * episode / epsilon_decay_episodes)
        
        for step in range(max_steps_per_episode):
            action_idx = agent.select_action(current_state, epsilon)
            action = MOVE_TYPES[action_idx]
            new_positions = apply_move(current_positions, action, hp_string)
            
            if new_positions is None:
                reward, next_state, next_energy = -2.0, current_state, current_energy
            else:
                next_energy = calculate_energy(new_positions, hp_string)
                reward = float(current_energy - next_energy)
                if next_energy < best_energy: reward += 5.0
                next_state = positions_to_grid(new_positions, hp_string)
            
            agent.memory.push(current_state, action_idx, reward, next_state)
            loss = agent.train_step(batch_size)
            if loss > 0: losses.append(loss)
            
            if new_positions is not None:
                current_positions, current_energy, current_state = new_positions, next_energy, next_state
                if current_energy < best_energy:
                    best_energy, best_positions = current_energy, list(current_positions)
                    
        if episode % target_update_freq == 0:
            agent.update_target_network()
            
        if (episode + 1) % 50 == 0:
            avg_loss = sum(losses[-50:]) / 50 if losses else 0
            print(f"  Episode {episode+1}/{episodes} | Best Energy: {best_energy} | Eps: {epsilon:.2f} | Avg Loss: {avg_loss:.4f}")
            
    return best_positions, best_energy

if __name__ == '__main__':
    print("=== Google Colab DQN Optimization ===")
    
    # Inizializziamo l'agente una sola volta per generalizzare
    shared_agent = DQNAgent(num_actions=len(MOVE_TYPES), lr=5e-4)

    for name, seq in PROTEINS:
        shared_agent.memory.buffer.clear() # Reset memory to avoid grid size mismatch crashes
        hp_str = sequence_to_hp(seq)
        
        print(f"\nEvaluating {name} - Length: {len(seq)} AA")
        start = time.time()
        
        best_positions, best_energy = generate_2d_structure_dqn(
            hp_str,
            agent=shared_agent,
            episodes=100, 
            max_steps_per_episode=150,
            batch_size=64, 
            target_update_freq=20
        )
        
        elapsed = time.time() - start
        print(f"  -> DQN Best Energy: {best_energy} | Time: {elapsed:.2f}s")
        
        torch.save(shared_agent.policy_net.state_dict(), 'dqn_weights.pth')
        print(f"  [+] Weights saved to dqn_weights.pth")
