"""
Advanced Deep Q-Learning (DQN) for 2D HP Protein Folding
[Constructive Approach with Relative Directions]

This script is designed to run seamlessly in Google Colab.
Unlike standard perturbative approaches (che partono da una catena dritta e la piegano),
questo agente "costruisce" la proteina un amminoacido alla volta. 

Ad ogni passo, l'agente decide dove posizionare il prossimo amminoacido
usando direzioni relative (0: Avanti, 1: Sinistra, 2: Destra).
Questa tecnica rende lo stato invariante a traslazioni e rotazioni,
risolvendo il problema della "curse of dimensionality" e permettendo
alla rete neurale di generalizzare molto meglio.
"""

import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import numpy as np
import random
import time
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
        # L'ambiente parte sempre posizionando i primi 2 amminoacidi
        # per definire la direzione iniziale (senza perdita di generalita')
        self.positions = [(0, 0), (1, 0)]
        self.current_step = 2
        self.current_dir = (1, 0) # Guardiamo verso Est (+X)
        self.done = False
        self.total_energy = 0
        return self._get_state()
        
    def step(self, action):
        """
        Actions: 
        0 = Forward (dritto)
        1 = Left (gira a sinistra)
        2 = Right (gira a destra)
        """
        dx, dy = self.current_dir
        
        if action == 0:
            new_dir = (dx, dy)
        elif action == 1:
            new_dir = (-dy, dx) # Rotazione 90 gradi antioraria
        elif action == 2:
            new_dir = (dy, -dx) # Rotazione 90 gradi oraria
            
        new_pos = (self.positions[-1][0] + new_dir[0], self.positions[-1][1] + new_dir[1])
        
        # Controllo Auto-intersezione (Collisione)
        if new_pos in self.positions:
            self.done = True
            # Forte penalita' per collisione, l'episodio finisce in anticipo
            return self._get_state(), -1.0, self.done 
            
        # Mossa valida: posiziona l'amminoacido
        self.positions.append(new_pos)
        self.current_dir = new_dir
        
        # Calcolo dei nuovi contatti H-H formati
        reward = 0.0
        if self.hp_string[self.current_step] == 'H':
            # Check con tutti gli amminoacidi tranne l'ultimo (con cui e' legato)
            for i, p in enumerate(self.positions[:-2]):
                if self.hp_string[i] == 'H':
                    dist = abs(p[0] - new_pos[0]) + abs(p[1] - new_pos[1])
                    if dist == 1:
                        reward += 1.0 # +1 per ogni nuovo legame idrofobico
                        self.total_energy -= 1 # L'energia reale e' negativa
                        
        self.current_step += 1
        if self.current_step == self.n:
            self.done = True
            # Extra reward opzionale per aver finito la proteina senza collidere
            reward += 0.5 
            
        return self._get_state(), reward, self.done

    def _get_state(self):
        """
        Rappresentazione dello stato:
        Griglia 2D locale centrata sulla "testa" della catena e orientata 
        nella direzione corrente. Questo e' il trucco magico che rende 
        lo stato invariante alla posizione e alla rotazione globale!
        
        Canali:
        0: Posizioni degli 'H' gia' piazzati
        1: Posizioni dei 'P' gia' piazzati
        """
        state = np.zeros((2, self.grid_size, self.grid_size), dtype=np.float32)
        head_x, head_y = self.positions[-1]
        dx, dy = self.current_dir
        
        for i, (x, y) in enumerate(self.positions):
            # Coordinate relative alla testa
            rel_x = x - head_x
            rel_y = y - head_y
            
            # Rotazione in modo che la direzione corrente punti sempre in alto (+Y)
            # Matematica della rotazione in base alla direzione corrente (dx, dy)
            rot_x = rel_x * dy - rel_y * dx
            rot_y = rel_x * dx + rel_y * dy
            
            # Mappa sulla griglia centrata
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
        
        # Calcolo dimensione flat
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

def train_dqn(hp_sequences, episodes=1500, batch_size=64):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Training on device: {device}")
    
    if isinstance(hp_sequences, str):
        hp_sequences = [hp_sequences]
        
    # Inizializza reti
    policy_net = DQN_CNN(grid_size=21).to(device)
    target_net = DQN_CNN(grid_size=21).to(device)
    target_net.load_state_dict(policy_net.state_dict())
    target_net.eval()
    
    optimizer = optim.Adam(policy_net.parameters(), lr=0.001)
    memory = ReplayBuffer(10000)
    
    # Iperparametri RL
    gamma = 0.99
    epsilon_start = 1.0
    epsilon_end = 0.05
    epsilon_decay = 0.995
    epsilon = epsilon_start
    target_update = 10 # Aggiorna la target network ogni X episodi
    
    best_energy = {}
    history_energy = []
    
    start_time = time.time()
    
    for episode in range(episodes):
        hp_string = random.choice(hp_sequences)
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
                
                # Calcolo Q(s, a)
                q_values = policy_net(states).gather(1, actions)
                
                # Calcolo Target = R + gamma * max Q(s', a')
                with torch.no_grad():
                    next_q_values = target_net(next_states).max(1)[0].unsqueeze(1)
                    target = rewards + gamma * next_q_values * (1 - dones)
                    
                loss = F.mse_loss(q_values, target)
                
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                
        # Aggiornamento Epsilon e Target Net
        epsilon = max(epsilon_end, epsilon * epsilon_decay)
        if episode % target_update == 0:
            target_net.load_state_dict(policy_net.state_dict())
            
        # Logging
        current_energy = env.total_energy
        history_energy.append(current_energy)
        
        if hp_string not in best_energy or current_energy < best_energy[hp_string]:
            best_energy[hp_string] = current_energy
            
        if (episode + 1) % 50 == 0:
            elapsed = time.time() - start_time
            print(f"Episode {episode+1:4d} | Epsilon: {epsilon:.2f} | "
                  f"Seq Len: {env.n} | Best for Seq: {best_energy[hp_string]} | "
                  f"Last Energy: {current_energy} | Len: {len(env.positions)}/{env.n} | Time: {elapsed:.1f}s")
                  
    return policy_net, history_energy, best_energy

# ==========================================
# 5. EXECUTION & VISUALIZATION
# ==========================================

def evaluate_dqn(policy_net, hp_string):
    """Valuta l'agente senza esplorazione (Greedy)."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    env = ConstructiveHPEnv(hp_string, grid_size=21)
    state = env.reset()
    
    while not env.done:
        with torch.no_grad():
            state_tensor = torch.FloatTensor(state).unsqueeze(0).to(device)
            q_values = policy_net(state_tensor)
            action = q_values.max(1)[1].item()
            
        state, reward, done = env.step(action)
        
    return env.total_energy, len(env.positions)

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
    training_proteins = [sequence_to_hp(seq) for name, seq in PROTEINS]
    
    print(f"Starting Constructive DQN Training on a dataset of {len(training_proteins)} proteins...")
    
    # Alleniamo l'agente su diverse proteine
    trained_net, energy_hist, best_E_dict = train_dqn(training_proteins, episodes=20000, batch_size=128)
    
    print("\nTraining Completed!")
    print("Best energies found during training for each sequence:")
    for seq, energy in best_E_dict.items():
        print(f"  Len {len(seq):2d}: {energy}")
        
    torch.save(trained_net.state_dict(), 'dqn_advanced_weights.pth')
    print("Model weights saved to 'dqn_advanced_weights.pth'")
        
    print("\n--- Evaluation on Unseen/Large Protein ---")
    test_hp = "HPHHHPPHPPPPHPHPHPPPPPHPPHPHPHPPPPPHPPPPPPHHHHPPPHPPPPPHPPHPHPPPPPHPHHHPHPPP"
    print(f"Testing Constructive DQN on sequence length {len(test_hp)}...")
    
    test_energy, test_len = evaluate_dqn(trained_net, test_hp)
    print(f"Evaluation Result - Energy: {test_energy} | Length reached: {test_len}/{len(test_hp)}")
    
    # Plotting learning curve (Energia varia tra proteine, quindi plotto la media mobile per trend)
    plt.figure(figsize=(10, 5))
    plt.plot(energy_hist, alpha=0.3, color='blue', label='Energy per episode (Mixed lengths)')
    
    # Smoothing con media mobile
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
