"""
agent_dqn.py
------------
Implements the Deep Q-Network (DQN) agent for 2D HP protein folding.

Instead of a lookup table (as in Tabular Q-Learning), a Convolutional
Neural Network (CNN) is used to approximate Q-values. The CNN takes a
200×200 spatial grid of the current protein conformation as input and
outputs one Q-value per possible move action.

Key components:
  - DQNCNN     : the CNN neural network architecture (policy + target net)
  - ReplayBuffer: experience replay memory buffer
  - DQNAgent   : wraps the two networks, optimizer, and memory
  - run_dqn_inference : hybrid inference (DRL exploration + greedy refinement)
"""

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import random
from collections import deque


# ─────────────────────────────────────────────────────────────────────────────
# Neural Network Architecture
# ─────────────────────────────────────────────────────────────────────────────

class DQNCNN(nn.Module):
    """
    Convolutional Neural Network that approximates Q(s, a) for all 4 actions.

    Architecture:
      Input  : (batch, 1, 200, 200) — single-channel 200×200 grid
      Conv1  : 1→16 channels, 3×3 kernel, padding=1 → (batch, 16, 200, 200)
               + ReLU + MaxPool(2)                  → (batch, 16, 100, 100)
      Conv2  : 16→32 channels, 3×3 kernel, padding=1 → (batch, 32, 100, 100)
               + ReLU + MaxPool(2)                  → (batch, 32,  50,  50)
      Conv3  : 32→64 channels, 3×3 kernel, padding=1 → (batch, 64,  50,  50)
               + ReLU + MaxPool(2)                  → (batch, 64,  25,  25)
      Flatten: 64×25×25 = 40,000 features
      FC1    : Linear(40000 → 256) + ReLU
      FC2    : Linear(256 → num_actions)  ← Q-values for each action

    The CNN is particularly suited to this problem because protein folding
    patterns are translation-invariant (a cluster of H residues provides the
    same energy regardless of its position on the grid).
    """

    def __init__(self, num_actions=4):
        """
        Parameters
        ----------
        num_actions : int  Number of output Q-values (= number of move types)
        """
        super(DQNCNN, self).__init__()

        # ── Convolutional feature extractor ───────────────────────────────
        # 3 blocks of: Conv2d → ReLU → MaxPool2d(2)
        # Each MaxPool halves the spatial dimensions (200→100→50→25)
        self.conv_layers = nn.Sequential(
            # Block 1: extract low-level spatial features (edges, small patterns)
            nn.Conv2d(1, 16, kernel_size=3, stride=1, padding=1),  # 1 channel in, 16 out
            nn.ReLU(),                  # non-linearity (zero out negatives)
            nn.MaxPool2d(2),            # halve spatial size: 200×200 → 100×100

            # Block 2: extract mid-level features (combinations of patterns)
            nn.Conv2d(16, 32, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(2),            # 100×100 → 50×50

            # Block 3: extract high-level features (H-H cluster representations)
            nn.Conv2d(32, 64, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(2)             # 50×50 → 25×25
        )

        # ── Fully connected head ──────────────────────────────────────────
        # Takes the flattened CNN output (64×25×25 = 40,000) and produces
        # one Q-value per action
        self.fc_layers = nn.Sequential(
            nn.Linear(64 * 25 * 25, 256),  # 40,000 features → 256 hidden units
            nn.ReLU(),                      # non-linearity
            nn.Linear(256, num_actions)     # 256 → Q-value for each action
        )

    def forward(self, x):
        """
        Forward pass: maps a batch of grid states to Q-value vectors.

        Parameters
        ----------
        x : torch.Tensor  Shape (batch, 1, 200, 200)

        Returns
        -------
        torch.Tensor  Shape (batch, num_actions) — Q-values for each action
        """
        x = self.conv_layers(x)       # extract spatial features
        x = x.view(x.size(0), -1)     # flatten: (batch, 64, 25, 25) → (batch, 40000)
        x = self.fc_layers(x)         # compute Q-values
        return x


# ─────────────────────────────────────────────────────────────────────────────
# Experience Replay Buffer
# ─────────────────────────────────────────────────────────────────────────────

class ReplayBuffer:
    """
    Fixed-size circular buffer that stores past (s, a, r, s') transitions.

    Experience Replay serves two purposes:
      1. BREAKS temporal correlations: consecutive transitions are highly
         correlated; training on a random mini-batch reduces this.
      2. DATA EFFICIENCY: each transition can be reused multiple times.

    When the buffer is full, the oldest transitions are automatically
    discarded (deque with maxlen).
    """

    def __init__(self, capacity):
        """
        Parameters
        ----------
        capacity : int  Maximum number of transitions to store
        """
        # deque with maxlen automatically discards the oldest entry when full
        self.buffer = deque(maxlen=capacity)

    def push(self, state, action, reward, next_state):
        """
        Stores a single transition in the buffer.

        Parameters
        ----------
        state      : np.ndarray  Current grid state  (1, 200, 200)
        action     : int         Action index taken (0–3)
        reward     : float       Reward received
        next_state : np.ndarray  Next grid state     (1, 200, 200)
        """
        self.buffer.append((state, action, reward, next_state))

    def sample(self, batch_size):
        """
        Randomly samples a mini-batch of transitions from the buffer.

        Parameters
        ----------
        batch_size : int  Number of transitions to sample

        Returns
        -------
        states      : np.ndarray  (batch, 1, 200, 200)
        actions     : np.ndarray  (batch,)
        rewards     : np.ndarray  (batch,)
        next_states : np.ndarray  (batch, 1, 200, 200)
        """
        batch = random.sample(self.buffer, batch_size)  # random sample (no replacement)
        # Unzip the list of tuples into 4 separate tuples
        states, actions, rewards, next_states = zip(*batch)
        return (
            np.stack(states),                              # stack into array
            np.array(actions),
            np.array(rewards, dtype=np.float32),
            np.stack(next_states)
        )

    def __len__(self):
        """Returns the current number of transitions stored in the buffer."""
        return len(self.buffer)


# ─────────────────────────────────────────────────────────────────────────────
# State encoding: protein positions → 2D grid
# ─────────────────────────────────────────────────────────────────────────────

def positions_to_grid(positions, hp_string):
    """
    Converts the current list of lattice positions into a 2D numpy grid
    that the CNN can process as an input image.

    Encoding:
      +1.0  →  Hydrophobic (H) residue at this cell
      −1.0  →  Polar (P) residue at this cell
       0.0  →  Empty cell (no residue)

    The grid is centred on the first residue: residue 0 is always at
    cell (100, 100) in a 200×200 grid. This ensures the protein is
    always roughly centred regardless of its position in global coordinates.

    Parameters
    ----------
    positions  : list of (int, int)  Lattice coordinates of all residues
    hp_string  : str                 HP binary string

    Returns
    -------
    np.ndarray  Shape (1, 200, 200), dtype float32
    """
    grid_size = 200                                     # grid is 200×200
    grid = np.zeros((1, grid_size, grid_size), dtype=np.float32)  # start with all zeros

    if not positions:
        return grid  # empty sequence → return blank grid

    # Translate coordinates so residue 0 is at the centre of the grid (100,100)
    ox, oy = positions[0]
    for i, (x, y) in enumerate(positions):
        tx = x - ox + (grid_size // 2)   # translate x to grid coordinates
        ty = y - oy + (grid_size // 2)   # translate y to grid coordinates

        # Clip to stay within grid bounds (prevent index out of range)
        tx = max(0, min(grid_size - 1, tx))
        ty = max(0, min(grid_size - 1, ty))

        # Place the residue value in the grid
        if hp_string[i] == 'H':
            grid[0, tx, ty] = 1.0    # H residue: positive value
        else:
            grid[0, tx, ty] = -1.0   # P residue: negative value

    return grid


# ─────────────────────────────────────────────────────────────────────────────
# DQN Agent
# ─────────────────────────────────────────────────────────────────────────────

class DQNAgent:
    """
    Deep Q-Network agent that manages:
      - policy_net  : the main CNN used to SELECT actions and compute Q(s,a)
      - target_net  : a slowly-updated copy of policy_net used to compute
                      STABLE target Q-values during training
      - memory      : the experience replay buffer
      - optimizer   : Adam optimizer for gradient updates

    Key Design: RANDOMNESS IN RL vs Deterministic Search
    ────────────────────────────────────────────────────────────────────────
    Unlike HC (greedy search), RL agents REQUIRE randomness for learning:
    
    1. ε-greedy EXPLORATION: With probability ε, take a random action instead
       of the "best" action according to the policy. This forces the agent to:
       - Discover whether other moves are better (exploration vs exploitation trade-off)
       - Sample diverse trajectories that the network learns from
       
    2. Experience REPLAY RANDOMNESS: Mini-batches sampled randomly from the buffer
       break temporal correlations and improve learning stability.
       
    3. MOVE PARAMETER RANDOMNESS: Even when applying a deterministic move type
       (e.g., end_flip), the target residue and destination are chosen randomly
       from valid options. This is a stochastic move set, not deterministic.
    
    WITHOUT this randomness, DQN would:
      - Get stuck in local minima like greedy search
      - Overfit to specific patterns
      - Fail to explore alternative folding pathways
    
    The TARGET NETWORK is the stabilization trick: it prevents the network
    from chasing a constantly-moving target during training, which would
    cause divergence or oscillations.
    """

    def __init__(self, num_actions=4, lr=1e-3, gamma=0.95, buffer_size=10000):
        """
        Parameters
        ----------
        num_actions  : int    Number of discrete actions (default: 4 move types)
        lr           : float  Adam learning rate (default: 0.001)
        gamma        : float  Discount factor for future rewards (default: 0.95)
        buffer_size  : int    Replay buffer capacity (default: 10,000)
        """
        # Use CPU for compatibility (no GPU required)
        self.device = torch.device("cpu")
        self.num_actions = num_actions
        self.gamma = gamma  # discount factor γ

        # ── Networks ──────────────────────────────────────────────────────
        # Policy network: updated every training step, used to select actions
        self.policy_net = DQNCNN(num_actions).to(self.device)
        # Target network: updated periodically (every N steps), used only
        # for computing stable TD targets (not for action selection)
        self.target_net = DQNCNN(num_actions).to(self.device)
        # Initialize target with the SAME weights as the policy network
        self.target_net.load_state_dict(self.policy_net.state_dict())
        # Target network is NEVER trained directly → set to eval mode
        self.target_net.eval()

        # ── Optimizer ─────────────────────────────────────────────────────
        # Adam adapts the learning rate per parameter → more stable than SGD
        self.optimizer = optim.Adam(self.policy_net.parameters(), lr=lr)

        # ── Replay Buffer ─────────────────────────────────────────────────
        self.memory = ReplayBuffer(buffer_size)

    def select_action(self, state_grid, epsilon):
        """
        Selects an action using the ε-greedy policy.
        
        This is the CORE EXPLORATION MECHANISM of DQN:
          - With probability ε : pick a RANDOM action (exploration)
          - With probability 1−ε: pick the action with the HIGHEST Q-value (exploitation)
        
        Why ε-greedy randomness is ESSENTIAL:
          Without it, the agent would always exploit (greedily choose the "best" move),
          which causes the same problems as Hill Climbing: getting stuck in local minima.
          By occasionally taking random actions, the agent:
            1. Discovers whether "suboptimal" moves lead to better long-term outcomes
            2. Samples diverse conformations and trajectory data for training
            3. Averages over different exploration paths, improving generalization
        
        ε is typically DECAYED over time: start high (100% random) to explore,
        gradually decrease to exploit the learned policy. This balances exploration
        early in training with exploitation of learned knowledge later.

        Parameters
        ----------
        state_grid : np.ndarray  Current grid state (1, 200, 200)
        epsilon    : float       Current exploration probability ε ∈ [0, 1]

        Returns
        -------
        int  Action index (0–3)
        """
        if random.random() < epsilon:
            # EXPLORE: random action (uniform over all moves)
            return random.randrange(self.num_actions)
        else:
            # EXPLOIT: use the policy network to find the best action
            with torch.no_grad():  # no gradient needed for inference
                # Add batch dimension: (1, 200, 200) → (1, 1, 200, 200)
                state_tensor = torch.FloatTensor(state_grid).unsqueeze(0).to(self.device)
                # Forward pass through the CNN to get Q-values
                q_values = self.policy_net(state_tensor)
                # Return the index of the highest Q-value action
                return q_values.argmax().item()

    def train_step(self, batch_size):
        """
        Performs one gradient descent step on the DQN loss.

        The loss is the Mean Squared Error between:
          - Predicted Q(s, a)  from the POLICY network
          - Target Q = r + γ · max_a' Q(s', a')  from the TARGET network

        Parameters
        ----------
        batch_size : int  Number of transitions to sample from the replay buffer

        Returns
        -------
        float  Loss value for this step (0.0 if buffer too small)
        """
        if len(self.memory) < batch_size:
            # Not enough transitions in buffer yet → skip training
            return 0.0

        # ── Sample a random mini-batch from the replay buffer ─────────────
        states, actions, rewards, next_states = self.memory.sample(batch_size)

        # Convert numpy arrays to PyTorch tensors on the device
        states      = torch.FloatTensor(states).to(self.device)        # (batch, 1, 200, 200)
        actions     = torch.LongTensor(actions).unsqueeze(1).to(self.device)  # (batch, 1)
        rewards     = torch.FloatTensor(rewards).unsqueeze(1).to(self.device) # (batch, 1)
        next_states = torch.FloatTensor(next_states).to(self.device)   # (batch, 1, 200, 200)

        # ── Predicted Q(s, a) from the policy network ─────────────────────
        # policy_net(states) → (batch, num_actions)
        # .gather(1, actions) selects the Q-value for the action actually taken
        q_values = self.policy_net(states).gather(1, actions)  # (batch, 1)

        # ── Target Q-values from the target network ────────────────────────
        with torch.no_grad():  # no gradient through the target network
            # max Q(s', a') → greedy value of the next state
            next_q_values = self.target_net(next_states).max(1)[0].unsqueeze(1)  # (batch, 1)

        # Bellman target: r + γ · max_a' Q(s', a'; θ⁻)
        target_q_values = rewards + self.gamma * next_q_values  # (batch, 1)

        # ── Compute MSE loss and backpropagate ─────────────────────────────
        loss = nn.MSELoss()(q_values, target_q_values)

        self.optimizer.zero_grad()  # clear previous gradients
        loss.backward()             # compute gradients via backpropagation
        self.optimizer.step()       # update policy network weights

        return loss.item()  # return scalar loss value for logging

    def update_target_network(self):
        """
        Synchronizes the target network weights with the policy network.
        Called periodically (e.g., every 10–20 training steps) to provide
        stable Q-value targets and prevent training oscillations.
        """
        self.target_net.load_state_dict(self.policy_net.state_dict())


# ─────────────────────────────────────────────────────────────────────────────
# Hybrid DRL inference
# ─────────────────────────────────────────────────────────────────────────────

def run_dqn_inference(hp_string, weights_path, max_steps=600):
    """
    Runs a HYBRID DRL inference to fold a protein using pre-trained weights.

    The inference has two phases:
      Phase 1 – DRL-guided exploration (first half of steps, ε=0.4):
                The trained policy network selects moves with some randomness.
                This avoids getting stuck in the straight-line starting point
                which would happen with pure greedy (ε=0) inference.
      Phase 2 – Greedy local search (second half of steps, ε=0):
                From the best structure found in Phase 1, try all 4 moves
                at each step and pick the best improving one. If none improve,
                apply a random move to escape local minima.

    Parameters
    ----------
    hp_string    : str   HP binary sequence to fold
    weights_path : str   Path to the .pth file with trained CNN weights
    max_steps    : int   Total inference steps (default: 600)

    Returns
    -------
    best_positions : list of (int, int)  Best conformation found
    best_energy    : int                 Minimum energy found
    """
    import os
    # Import the shared move functions from the Q-Learning module
    from test_ql import apply_move, calculate_energy, MOVE_TYPES

    n = len(hp_string)
    # Start from a straight-line conformation (same as training)
    initial_positions = [(i, 0) for i in range(n)]
    best_positions = list(initial_positions)
    best_energy = calculate_energy(initial_positions, hp_string)

    # If the weights file doesn't exist, fall back to the straight-line
    if not os.path.exists(weights_path):
        return best_positions, best_energy

    # Load the pre-trained DQN agent
    agent = DQNAgent(num_actions=len(MOVE_TYPES))
    try:
        # Load the saved CNN weights into the policy network
        agent.policy_net.load_state_dict(
            torch.load(weights_path, map_location=agent.device)
        )
        agent.policy_net.eval()  # set to evaluation mode (disables dropout, etc.)
    except Exception as e:
        print(f"Failed to load weights: {e}")
        return best_positions, best_energy  # fall back to straight-line

    # ── Phase 1: DRL-guided exploration ───────────────────────────────────
    # Use the trained network with ε=0.4 (40% random actions) to explore
    # the conformation space and find a good starting point for Phase 2.
    explore_steps = max_steps // 2     # use half the budget for exploration
    current_positions = list(initial_positions)
    current_energy = best_energy

    for _ in range(explore_steps):
        # Encode current conformation as a 200×200 grid
        state_grid = positions_to_grid(current_positions, hp_string)
        # Select action with partial exploration (ε=0.4)
        action_idx = agent.select_action(state_grid, epsilon=0.4)
        action = MOVE_TYPES[action_idx]

        # Apply the selected move
        new_positions = apply_move(current_positions, action, hp_string)
        if new_positions is None:
            continue  # Invalid move → skip and try next step

        new_energy = calculate_energy(new_positions, hp_string)
        current_positions = new_positions
        current_energy = new_energy

        # Track the best conformation seen during exploration
        if current_energy < best_energy:
            best_energy = current_energy
            best_positions = list(current_positions)

    # ── Phase 2: Greedy local search from the best Phase-1 structure ──────
    # Deterministic refinement: try all 4 moves and pick the best one.
    refine_steps = max_steps - explore_steps  # remaining budget
    current_positions = list(best_positions)  # start from Phase 1 best

    for _ in range(refine_steps):
        found_better = False
        # Try every possible move type and pick the one that improves energy most
        for action in MOVE_TYPES:
            new_positions = apply_move(current_positions, action, hp_string)
            if new_positions is None:
                continue
            new_energy = calculate_energy(new_positions, hp_string)
            if new_energy < best_energy:
                # Found an improving move → accept it
                best_energy = new_energy
                best_positions = list(new_positions)
                current_positions = new_positions
                found_better = True
                break  # restart from the better conformation immediately

        if not found_better:
            # No improving move found → stuck in local minimum.
            # Apply a random move to escape (accepts any valid move)
            import random
            action = random.choice(MOVE_TYPES)
            new_positions = apply_move(current_positions, action, hp_string)
            if new_positions:
                current_positions = new_positions  # move regardless of energy

    return best_positions, best_energy
