from google.colab import drive
# Mount Google Drive so long-running checkpoints survive Colab disconnections.
drive.mount('/content/drive')

"""
Deep Q-Learning (DQN) for 2D HP Protein Folding
[Constructive Approach with Relative Directions]

Training features:
  1. Deeper CNN: 3 convolutional layers with a wider receptive field
  2. Larger replay buffer: 100,000 transitions (previously 10,000)
  3. Target-network update every 1,000 optimization STEPS, not episodes
  4. Gradient clipping (max_norm=10) for numerical stability
  5. Double DQN to reduce Q-value overestimation
  6. Reward shaping for proteins with sparse hydrophobic residues
  7. Uniform experience-replay sampling from a large cross-protein memory
"""

import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import numpy as np
import random
import time
import csv
import os
from collections import deque
import matplotlib.pyplot as plt

# ==========================================
# 1. ENVIRONMENT (Constructive HP Model)
# ==========================================

class ConstructiveHPEnv:
    """Construct a self-avoiding HP chain one residue at a time.

    The environment supplies states, rewards and terminal conditions to DQN.
    Its three relative actions are Forward, Left and Right.
    """

    def __init__(self, hp_string, grid_size=21):
        """Store sequence metadata and initialize the first episode."""

        self.hp_string = hp_string
        # Total number of residues that must be placed before completion.
        self.n = len(hp_string)
        # An odd size provides one unambiguous center cell for the chain head.
        self.grid_size = grid_size
        self.half_grid = grid_size // 2
        # Precompute H density; sparse-H proteins receive extra shaping guidance.
        self.h_density = hp_string.count('H') / len(hp_string)
        self.reset()

    def reset(self):
        """Reset geometry and return the initial normalized state grid."""

        # The first two residues establish a fixed initial direction.
        self.positions = [(0, 0), (1, 0)]
        # A set provides average O(1) collision lookup instead of O(n) list scan.
        self.position_set = {(0, 0), (1, 0)}
        # current_step is both the placed-residue count and next residue index.
        self.current_step = 2
        # Unit vector pointing from residue 0 to residue 1 (east/right).
        self.current_dir = (1, 0)
        # Training proteins contain more than two residues, so construction starts active.
        self.done = False
        # A straight chain has no non-bonded H-H contacts initially.
        self.total_energy = 0
        return self._get_state()

    def step(self, action):
        """
        Apply one relative action and return ``(next_state, reward, done)``.

        Actions:
        0 = Forward
        1 = Left (counterclockwise)
        2 = Right (clockwise)
        """
        dx, dy = self.current_dir

        if action == 0:
            # Forward preserves the current unit direction.
            new_dir = (dx, dy)
        elif action == 1:
            # A +90-degree rotation maps (dx,dy) to (-dy,dx).
            new_dir = (-dy, dx)
        elif action == 2:
            # A -90-degree rotation maps (dx,dy) to (dy,-dx).
            new_dir = (dy, -dx)

        new_pos = (self.positions[-1][0] + new_dir[0],
                   self.positions[-1][1] + new_dir[1])

        # A repeated lattice coordinate violates the self-avoiding-walk constraint.
        if new_pos in self.position_set:
            self.done = True
            return self._get_state(), -1.0, self.done

        self.positions.append(new_pos)
        # Keep list and set synchronized: list preserves residue order, while the
        # set accelerates collision and neighborhood membership tests.
        self.position_set.add(new_pos)
        self.current_dir = new_dir

        reward = 0.0

        # Only a newly placed H can create new hydrophobic contacts.
        if self.hp_string[self.current_step] == 'H':
            new_hh = 0
            # Exclude the new residue and its immediate backbone predecessor.
            for i, p in enumerate(self.positions[:-2]):
                if self.hp_string[i] == 'H':
                    dist = abs(p[0] - new_pos[0]) + abs(p[1] - new_pos[1])
                    if dist == 1:
                        # Every newly formed non-bonded H-H contact contributes -1 energy.
                        new_hh += 1
                        self.total_energy -= 1

            reward += float(new_hh)

            # ── REWARD SHAPING ──────────────────────────────────────────────
            # Low-H-density proteins produce sparse rewards. Add a small bonus
            # for proximity to already placed H residues. This changes the
            # learning signal but never changes the reported physical HP energy.
            if self.h_density < 0.35:
                h_neighbors = 0
                for dxn, dyn in [(0,1),(0,-1),(1,0),(-1,0)]:
                    nb = (new_pos[0]+dxn, new_pos[1]+dyn)
                    if nb in self.position_set:
                        # Recover the sequence index so the neighbor type is known.
                        idx = self.positions.index(nb)
                        if self.hp_string[idx] == 'H':
                            h_neighbors += 1
                reward += 0.1 * h_neighbors
            # ────────────────────────────────────────────────────────────────

        self.current_step += 1
        if self.current_step == self.n:
            self.done = True
            reward += 0.5  # Completion bonus favors full collision-free folds.

        return self._get_state(), reward, self.done

    def _get_state(self):
        """
        Build a 2x21x21 grid centered on the chain head.

        Coordinates are rotation-normalized so the current direction always
        has the same canonical orientation from the CNN's perspective.
        """
        # Channel 0 stores H occupancy; channel 1 stores P occupancy.
        state = np.zeros((2, self.grid_size, self.grid_size), dtype=np.float32)
        head_x, head_y = self.positions[-1]
        dx, dy = self.current_dir

        for i, (x, y) in enumerate(self.positions):
            # First translate every residue relative to the current chain head.
            rel_x = x - head_x
            rel_y = y - head_y

            # Then rotate by the current direction to remove global orientation.
            rot_x = rel_x * dy - rel_y * dx
            rot_y = rel_x * dx + rel_y * dy

            # Shift local coordinates into non-negative NumPy array indices.
            grid_x = rot_x + self.half_grid
            grid_y = rot_y + self.half_grid

            if 0 <= grid_x < self.grid_size and 0 <= grid_y < self.grid_size:
                channel = 0 if self.hp_string[i] == 'H' else 1
                state[channel, int(grid_y), int(grid_x)] = 1.0

        return state


# ==========================================
# 2. NEURAL NETWORK — Deeper CNN
# ==========================================

class DQN_CNN(nn.Module):
    """
    Improvements over the previous architecture:
    - 3 convolutional layers instead of 2: 7x7 receptive field instead of 5x5
    - More feature maps in intermediate layers
    - Larger fully connected representation: 256 units instead of 128
    """
    def __init__(self, grid_size=21):
        super(DQN_CNN, self).__init__()

        # Padding=1 preserves 21x21 spatial dimensions after each 3x3 convolution.
        self.conv1 = nn.Conv2d(2,  32, kernel_size=3, stride=1, padding=1)  # → 32×21×21
        self.conv2 = nn.Conv2d(32, 64, kernel_size=3, stride=1, padding=1)  # → 64×21×21
        self.conv3 = nn.Conv2d(64, 64, kernel_size=3, stride=1, padding=1)  # → 64×21×21

        self.flat_size = 64 * grid_size * grid_size  # 28,224

        self.fc1 = nn.Linear(self.flat_size, 256)
        self.fc2 = nn.Linear(256, 3)  # Three Q-values: Forward, Left, Right.

    def forward(self, x):
        """Map a batch of 2-channel state grids to three Q-values each."""

        x = F.relu(self.conv1(x))
        x = F.relu(self.conv2(x))
        x = F.relu(self.conv3(x))
        x = x.view(x.size(0), -1)
        x = F.relu(self.fc1(x))
        return self.fc2(x)


# ==========================================
# 3. REPLAY BUFFER — 100k capacity
# ==========================================

class ReplayBuffer:
    """
    Replay memory expanded to 100,000 transitions.

    With 20 proteins and chains of 56-98 residues, the larger memory delays
    overwriting old experience and helps reduce catastrophic forgetting.
    """
    def __init__(self, capacity=100_000):
        # deque automatically removes the oldest item when capacity is reached.
        self.buffer = deque(maxlen=capacity)

    def push(self, state, action, reward, next_state, done):
        """Store one complete transition tuple for later reuse."""
        self.buffer.append((state, action, reward, next_state, done))

    def sample(self, batch_size):
        """Return a uniformly sampled mini-batch as NumPy arrays."""

        # Random replay breaks temporal correlations between adjacent steps.
        batch = random.sample(self.buffer, batch_size)
        states, actions, rewards, next_states, dones = zip(*batch)
        return (np.array(states),
                np.array(actions),
                np.array(rewards, dtype=np.float32),
                np.array(next_states),
                np.array(dones, dtype=np.float32))

    def __len__(self):
        return len(self.buffer)


# ==========================================
# 4. UTILITY
# ==========================================

def prepare_training_items(hp_sequences):
    """Normalize supported inputs into ``[(name, hp_string), ...]``."""

    if isinstance(hp_sequences, str):
        # A single anonymous sequence receives a generated display name.
        return [("Sequence_1", hp_sequences)]
    training_items = []
    for idx, item in enumerate(hp_sequences, start=1):
        if isinstance(item, (tuple, list)) and len(item) >= 2:
            # Preserve an explicitly supplied protein name and sequence.
            name, hp_string = item[0], item[1]
        else:
            # Bare sequence elements receive deterministic generated names.
            name, hp_string = f"Sequence_{idx}", item
        training_items.append((str(name), str(hp_string)))
    return training_items


def init_sequence_stats(name, hp_string):
    """Create zeroed training accumulators for one protein."""

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


def atomic_torch_save(obj, path):
    """Save a torch object through a temporary file, then replace atomically."""
    tmp_path = f"{path}.tmp"
    # Write completely to a sibling file before replacing the public checkpoint.
    torch.save(obj, tmp_path)
    # os.replace is atomic on the same filesystem, preventing half-written files.
    os.replace(tmp_path, path)


def atomic_write_dqn_reports(sequence_stats, episode_logs, greedy_results,
                             training_config, summary_path, episode_log_path,
                             report_path):
    """Write reports safely so an interrupted write does not corrupt old files."""
    tmp_summary = f"{summary_path}.tmp"
    tmp_episode = f"{episode_log_path}.tmp"
    tmp_report = f"{report_path}.tmp"
    # Generate all three files under temporary names first.
    write_dqn_reports(
        sequence_stats,
        episode_logs,
        greedy_results,
        training_config,
        summary_path=tmp_summary,
        episode_log_path=tmp_episode,
        report_path=tmp_report,
    )
    # Publish only successfully completed report files.
    os.replace(tmp_summary, summary_path)
    os.replace(tmp_episode, episode_log_path)
    os.replace(tmp_report, report_path)


def save_training_artifacts(policy_net, target_net, optimizer, memory,
                            sequence_stats, episode_logs, history_energy,
                            best_energy, training_config, checkpoint_dir,
                            episode, epsilon, total_opt_steps, start_time,
                            save_replay_buffer=False,
                            save_numbered_checkpoints=False,
                            final=False):
    """
    Save enough state to keep useful results if Colab disconnects.

    The lightweight checkpoint stores model/optimizer state and logs. The replay
    buffer is optional because with 100k transitions it can become very large.
    """
    os.makedirs(checkpoint_dir, exist_ok=True)

    # Include time from earlier resumed sessions in the accumulated duration.
    elapsed_s = round(time.time() - start_time, 3)
    checkpoint_config = dict(training_config)
    checkpoint_config["training_time_s"] = elapsed_s
    checkpoint_config["total_opt_steps"] = total_opt_steps
    checkpoint_config["last_saved_episode"] = episode
    checkpoint_config["checkpoint_complete"] = bool(final)
    checkpoint_config["replay_buffer_saved"] = bool(save_replay_buffer)

    checkpoint = {
        # Progress markers needed to continue from the same training point.
        "episode": episode,
        "epsilon": epsilon,
        "total_opt_steps": total_opt_steps,
        # Learnable network and optimizer states.
        "policy_state_dict": policy_net.state_dict(),
        "target_state_dict": target_net.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "best_energy": best_energy,
        "sequence_stats": sequence_stats,
        "episode_logs": episode_logs,
        "history_energy": history_energy,
        "training_config": checkpoint_config,
        # Random-generator states make resumed stochastic behavior reproducible.
        "python_random_state": random.getstate(),
        "numpy_random_state": np.random.get_state(),
        "torch_random_state": torch.get_rng_state(),
    }
    if torch.cuda.is_available():
        # CUDA maintains separate generator state for each visible GPU.
        checkpoint["torch_cuda_random_state_all"] = torch.cuda.get_rng_state_all()
    if save_replay_buffer:
        # Optional because serializing 100k high-dimensional states is expensive.
        checkpoint["replay_buffer"] = list(memory.buffer)

    latest_checkpoint = os.path.join(checkpoint_dir, "dqn_checkpoint_latest.pth")
    latest_weights = os.path.join(checkpoint_dir, "dqn_weights_latest.pth")
    # Always overwrite stable "latest" names for convenient recovery.
    atomic_torch_save(checkpoint, latest_checkpoint)
    atomic_torch_save(policy_net.state_dict(), latest_weights)

    if save_numbered_checkpoints or final:
        # Numbered snapshots preserve historical milestones instead of replacing them.
        episode_checkpoint = os.path.join(
            checkpoint_dir,
            f"dqn_checkpoint_ep{episode:05d}.pth",
        )
        episode_weights = os.path.join(
            checkpoint_dir,
            f"dqn_weights_ep{episode:05d}.pth",
        )
        atomic_torch_save(checkpoint, episode_checkpoint)
        atomic_torch_save(policy_net.state_dict(), episode_weights)

    atomic_write_dqn_reports(
        sequence_stats,
        episode_logs,
        greedy_results={},
        training_config=checkpoint_config,
        summary_path=os.path.join(checkpoint_dir, "dqn_summary_latest.csv"),
        episode_log_path=os.path.join(checkpoint_dir, "dqn_episode_log_latest.csv"),
        report_path=os.path.join(checkpoint_dir, "dqn_report_latest.txt"),
    )

    status = "final" if final else "checkpoint"
    print(f"\n[{status}] Saved episode {episode} artifacts in '{checkpoint_dir}'")
    print(f"  checkpoint: {latest_checkpoint}")
    print(f"  weights:    {latest_weights}\n")


def load_training_checkpoint(checkpoint_path, policy_net, target_net, optimizer,
                             memory, device):
    """Resume model, optimizer, logs, and optional replay buffer from checkpoint."""
    try:
        checkpoint = torch.load(
            checkpoint_path,
            map_location=device,
            weights_only=False,
        )
    except TypeError:
        checkpoint = torch.load(checkpoint_path, map_location=device)
    # Restore learned parameters and Adam's running moment estimates.
    policy_net.load_state_dict(checkpoint["policy_state_dict"])
    target_net.load_state_dict(checkpoint["target_state_dict"])
    optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
    for state in optimizer.state.values():
        # map_location moves model tensors, but optimizer internals may still need
        # explicit transfer to the active device.
        for key, value in state.items():
            if torch.is_tensor(value):
                state[key] = value.to(device)

    if "replay_buffer" in checkpoint:
        # Replace current memory only when the large optional buffer was saved.
        memory.buffer.clear()
        memory.buffer.extend(checkpoint["replay_buffer"])

    if "python_random_state" in checkpoint:
        # Restoring all random sources keeps exploration sequences reproducible.
        random.setstate(checkpoint["python_random_state"])
    if "numpy_random_state" in checkpoint:
        np.random.set_state(checkpoint["numpy_random_state"])
    if "torch_random_state" in checkpoint:
        torch.set_rng_state(checkpoint["torch_random_state"].cpu())
    if torch.cuda.is_available() and "torch_cuda_random_state_all" in checkpoint:
        cuda_states = [
            state.cpu() for state in checkpoint["torch_cuda_random_state_all"]
        ]
        torch.cuda.set_rng_state_all(cuda_states)

    return checkpoint


# ==========================================
# 5. TRAINING LOOP — complete stabilized version
# ==========================================

def train_dqn(hp_sequences, episodes=20000, batch_size=128,
              checkpoint_interval=500,
              checkpoint_dir="dqn_checkpoints",
              resume_from=None,
              save_replay_buffer=False,
              save_numbered_checkpoints=False):
    """Train one shared constructive Double-DQN over multiple HP proteins.

    Each episode samples one protein, constructs it until completion/collision,
    stores transitions in replay memory and performs online mini-batch updates.
    The function supports atomic checkpoints and exact random-state recovery.
    """

    # Prefer a CUDA GPU in Colab, while preserving CPU compatibility.
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Training on device: {device}")

    training_items = prepare_training_items(hp_sequences)

    # ── Policy and target networks ─────────────────────────────────────────
    # policy_net is optimized; target_net provides slowly changing Bellman targets.
    policy_net = DQN_CNN(grid_size=21).to(device)
    target_net = DQN_CNN(grid_size=21).to(device)
    target_net.load_state_dict(policy_net.state_dict())
    # Evaluation mode prevents training-specific behavior in the target network.
    target_net.eval()

    # Adam adapts per-parameter learning rates using first/second moments.
    optimizer = optim.Adam(policy_net.parameters(), lr=0.001)

    # ── Large replay buffer ────────────────────────────────────────────────
    memory = ReplayBuffer(capacity=100_000)

    # ── Hyperparameters ────────────────────────────────────────────────────
    # gamma controls how strongly future Q-values influence the current target.
    gamma          = 0.99
    # Epsilon starts fully exploratory and decays multiplicatively per episode.
    epsilon_start  = 1.0
    epsilon_end    = 0.05
    epsilon_decay  = 0.995
    epsilon        = epsilon_start

    # Update target by optimization-step count, independent of episode length.
    target_update_steps = 1_000
    total_opt_steps     = 0

    # Per-protein minima, aggregates and detailed logs are maintained separately.
    best_energy    = {}
    sequence_stats = {}
    episode_logs   = []
    history_energy = []
    start_episode  = 0
    previous_training_time_s = 0.0

    # Store enough metadata to interpret reports and reproduce the experiment.
    training_config = {
        "episodes": episodes,
        "batch_size": batch_size,
        "grid_size": 21,
        "replay_buffer_capacity": 100_000,
        "learning_rate": 0.001,
        "gamma": gamma,
        "epsilon_start": epsilon_start,
        "epsilon_end": epsilon_end,
        "epsilon_decay": epsilon_decay,
        "target_update_steps": target_update_steps,
        "actions": "Forward,Left,Right",
        "device": str(device),
        "num_training_proteins": len(training_items),
        "double_dqn": True,
        "gradient_clipping": 10,
        "reward_shaping": True,
        "checkpoint_interval": checkpoint_interval,
        "checkpoint_dir": checkpoint_dir,
        "resume_from": resume_from,
        "save_replay_buffer": save_replay_buffer,
        "save_numbered_checkpoints": save_numbered_checkpoints,
    }

    if resume_from:
        # Resuming is explicit: a missing requested checkpoint is a hard error.
        if not os.path.exists(resume_from):
            raise FileNotFoundError(f"Checkpoint not found: {resume_from}")

        checkpoint = load_training_checkpoint(
            resume_from,
            policy_net,
            target_net,
            optimizer,
            memory,
            device,
        )
        # .get provides backward compatibility with older checkpoint formats.
        start_episode = int(checkpoint.get("episode", 0))
        epsilon = float(checkpoint.get("epsilon", epsilon))
        total_opt_steps = int(checkpoint.get("total_opt_steps", total_opt_steps))
        best_energy = checkpoint.get("best_energy", best_energy)
        sequence_stats = checkpoint.get("sequence_stats", sequence_stats)
        episode_logs = checkpoint.get("episode_logs", episode_logs)
        history_energy = checkpoint.get("history_energy", history_energy)
        previous_config = checkpoint.get("training_config", {})
        previous_training_time_s = float(previous_config.get("training_time_s", 0.0))
        print(f"Resumed from '{resume_from}' at episode {start_episode}.")
        print(f"Previous optimisation steps: {total_opt_steps:,}")
        print(f"Replay buffer restored: {len(memory):,} transitions")

    # Shift the new timer backwards so reported duration includes earlier sessions.
    start_time = time.time() - previous_training_time_s

    try:
        for episode in range(start_episode, episodes):
            # Uniform random protein sampling interleaves experience and trains
            # one shared policy rather than one separate model per sequence.
            protein_name, hp_string = random.choice(training_items)
            env = ConstructiveHPEnv(hp_string, grid_size=21)
            state = env.reset()
            total_reward = 0.0

            while not env.done:
                # ── Epsilon-Greedy ──────────────────────────────────────────
                if random.random() < epsilon:
                    # Exploration samples Forward/Left/Right uniformly.
                    action = random.randrange(3)
                else:
                    # Exploitation uses the current policy's largest Q-value.
                    with torch.no_grad():
                        state_t = torch.FloatTensor(state).unsqueeze(0).to(device)
                        q_values = policy_net(state_t)
                        action = q_values.max(1)[1].item()

                next_state, reward, done = env.step(action)
                # Save the transition before advancing the local state variable.
                memory.push(state, action, reward, next_state, done)

                state = next_state
                total_reward += reward

                # ── Optimization ───────────────────────────────────────────
                # Delay learning until one complete mini-batch can be sampled.
                if len(memory) > batch_size:
                    states_b, actions_b, rewards_b, next_states_b, dones_b = \
                        memory.sample(batch_size)

                    # Convert NumPy batches to device tensors with shapes:
                    # states (B,2,21,21), scalar columns (B,1).
                    states_b      = torch.FloatTensor(states_b).to(device)
                    actions_b     = torch.LongTensor(actions_b).unsqueeze(1).to(device)
                    rewards_b     = torch.FloatTensor(rewards_b).unsqueeze(1).to(device)
                    next_states_b = torch.FloatTensor(next_states_b).to(device)
                    dones_b       = torch.FloatTensor(dones_b).unsqueeze(1).to(device)

                    # Current Q(s,a): gather selects the output corresponding to
                    # the action actually stored in each replay transition.
                    q_values_b = policy_net(states_b).gather(1, actions_b)

                    # ── DOUBLE DQN ──────────────────────────────────────────
                    # The policy network SELECTS the next action, while the target
                    # network EVALUATES it. This reduces max-Q overestimation.
                    with torch.no_grad():
                        next_actions = policy_net(next_states_b).max(1)[1].unsqueeze(1)
                        next_q = target_net(next_states_b).gather(1, next_actions)
                        # Terminal transitions multiply future value by zero.
                        target = rewards_b + gamma * next_q * (1 - dones_b)

                    # Regress current predictions toward Double-DQN Bellman targets.
                    loss = F.mse_loss(q_values_b, target)

                    # Standard PyTorch update: clear, backpropagate, clip, step.
                    optimizer.zero_grad()
                    loss.backward() # Compute gradients for all learnable parameters.
                    # Clipping prevents occasional large gradients from destabilizing training.
                    torch.nn.utils.clip_grad_norm_(policy_net.parameters(),
                                                   max_norm=10)
                    optimizer.step()

                    total_opt_steps += 1

                    # ── TARGET NETWORK UPDATE every 1,000 optimization steps ─
                    if total_opt_steps % target_update_steps == 0:
                        # Hard update: target becomes an exact policy copy.
                        target_net.load_state_dict(policy_net.state_dict())

            # ── Epsilon decay ──────────────────────────────────────────────
            # max enforces a permanent 5% exploration floor.
            epsilon = max(epsilon_end, epsilon * epsilon_decay)

            # ── Logging ────────────────────────────────────────────────────
            current_energy = env.total_energy
            reached_length = len(env.positions)
            completed_chain = reached_length == env.n
            history_energy.append(current_energy)

            # Update the best physical HP energy observed for this protein.
            if protein_name not in best_energy or \
                    current_energy < best_energy[protein_name]:
                best_energy[protein_name] = current_energy

            # Lazily initialize aggregates the first time a protein is sampled.
            if protein_name not in sequence_stats:
                sequence_stats[protein_name] = \
                    init_sequence_stats(protein_name, hp_string)

            # Accumulate values now; averages are computed only while reporting.
            stats = sequence_stats[protein_name]
            stats["episodes_seen"]    += 1
            stats["complete_episodes"] += int(completed_chain)
            stats["sum_energy"]        += current_energy
            stats["sum_reward"]        += total_reward
            stats["sum_length"]        += reached_length
            stats["last_energy"]        = current_energy
            stats["last_reward"]        = total_reward
            stats["last_length"]        = reached_length
            stats["best_length"]        = max(stats["best_length"], reached_length)

            if stats["best_energy"] is None or \
                    current_energy < stats["best_energy"]:
                stats["best_energy"] = current_energy
            if stats["best_reward"] is None or \
                    total_reward > stats["best_reward"]:
                stats["best_reward"] = total_reward

            # Preserve one detailed row per episode for later CSV analysis.
            episode_logs.append({
                "episode":                episode + 1,
                "protein":                protein_name,
                "length":                 env.n,
                "h_count":                hp_string.count("H"),
                "epsilon":                round(epsilon, 6),
                "total_reward":           round(total_reward, 6),
                "final_energy":           current_energy,
                "reached_length":         reached_length,
                "completed_chain":        int(completed_chain),
                "best_energy_for_protein": best_energy[protein_name],
                "total_opt_steps":        total_opt_steps,
                "elapsed_s":              round(time.time() - start_time, 3),
            })

            # Print compact progress without flooding the Colab output cell.
            if (episode + 1) % 50 == 0:
                elapsed = time.time() - start_time
                print(f"Episode {episode+1:5d} | ε: {epsilon:.3f} | "
                      f"Protein: {protein_name} | Len: {env.n} | "
                      f"Best: {best_energy[protein_name]:4d} | "
                      f"Reward: {total_reward:.2f} | "
                      f"Energy: {current_energy:4d} | "
                      f"Steps: {reached_length}/{env.n} | "
                      f"OptSteps: {total_opt_steps} | "
                      f"Time: {elapsed:.1f}s")

            # Periodically persist recoverable state to survive runtime resets.
            if checkpoint_interval and (episode + 1) % checkpoint_interval == 0:
                save_training_artifacts(
                    policy_net,
                    target_net,
                    optimizer,
                    memory,
                    sequence_stats,
                    episode_logs,
                    history_energy,
                    best_energy,
                    training_config,
                    checkpoint_dir,
                    episode=episode + 1,
                    epsilon=epsilon,
                    total_opt_steps=total_opt_steps,
                    start_time=start_time,
                    save_replay_buffer=save_replay_buffer,
                    save_numbered_checkpoints=save_numbered_checkpoints,
                    final=False,
                )
    except KeyboardInterrupt:
        # Ctrl+C or manual Colab interruption still produces usable artifacts.
        interrupted_episode = len(history_energy)
        save_training_artifacts(
            policy_net,
            target_net,
            optimizer,
            memory,
            sequence_stats,
            episode_logs,
            history_energy,
            best_energy,
            training_config,
            checkpoint_dir,
            episode=interrupted_episode,
            epsilon=epsilon,
            total_opt_steps=total_opt_steps,
            start_time=start_time,
            save_replay_buffer=save_replay_buffer,
            save_numbered_checkpoints=save_numbered_checkpoints,
            final=False,
        )
        print("Training interrupted. Latest checkpoint and partial reports were saved.")
        raise

    # Finalize metadata and save one checkpoint marked as complete.
    training_config["training_time_s"] = round(time.time() - start_time, 3)
    training_config["total_opt_steps"] = total_opt_steps
    training_config["completed_episodes"] = len(history_energy)
    save_training_artifacts(
        policy_net,
        target_net,
        optimizer,
        memory,
        sequence_stats,
        episode_logs,
        history_energy,
        best_energy,
        training_config,
        checkpoint_dir,
        episode=episodes,
        epsilon=epsilon,
        total_opt_steps=total_opt_steps,
        start_time=start_time,
        save_replay_buffer=save_replay_buffer,
        save_numbered_checkpoints=save_numbered_checkpoints,
        final=True,
    )
    return (policy_net, history_energy, best_energy,
            sequence_stats, episode_logs, training_config)


# ==========================================
# 6. EVALUATION
# ==========================================

def evaluate_dqn(policy_net, hp_string, n_rollouts=10):
    """
    Evaluate the trained policy with greedy constructive rollouts.

    The function returns the lowest energy observed and its reached length.
    The current policy/environment are deterministic under greedy argmax, so
    repeated rollouts normally reproduce the same result; the parameter remains
    useful if stochastic evaluation is introduced later.
    """
    # Use the same device-selection rule as training.
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # Evaluation mode disables training-specific behavior such as dropout.
    policy_net.eval()

    best_e    = 0
    best_len  = 0

    for _ in range(n_rollouts):
        # Every rollout starts from the canonical two-residue initial chain.
        env   = ConstructiveHPEnv(hp_string, grid_size=21)
        state = env.reset()

        while not env.done:
            # Greedy evaluation uses no epsilon exploration and no gradients.
            with torch.no_grad():
                state_t = torch.FloatTensor(state).unsqueeze(0).to(device)
                q_vals  = policy_net(state_t)
                action  = q_vals.max(1)[1].item()
            state, _, _ = env.step(action)

        # Lower (more negative) energy is considered better here. Note that the
        # current comparison does not prioritize completion over partial length.
        if env.total_energy < best_e:
            best_e   = env.total_energy
            best_len = len(env.positions)

    if best_len == 0:
        # Energy zero never beats the initial best_e=0, so preserve reached length.
        best_len = len(env.positions)

    return best_e, best_len


def evaluate_training_set(policy_net, training_items, n_rollouts=10):
    """Evaluate every training protein and return report-ready metrics."""

    results = {}
    for name, hp_string in training_items:
        energy, reached = evaluate_dqn(policy_net, hp_string,
                                       n_rollouts=n_rollouts)
        # Keep keys aligned with columns consumed by write_dqn_reports().
        results[name] = {
            "greedy_final_energy":    energy,
            "greedy_reached_length":  reached,
            "greedy_completed_chain": int(reached == len(hp_string)),
        }
    return results


# ==========================================
# 7. REPORT
# ==========================================

def write_dqn_reports(sequence_stats, episode_logs, greedy_results,
                      training_config,
                      summary_path="dqn_summary.csv",
                      episode_log_path="dqn_episode_log.csv",
                      report_path="dqn_report.txt"):
    """Write detailed episodes, per-protein summaries and a text report."""

    # The episode log preserves raw chronological training observations.
    episode_fields = [
        "episode", "protein", "length", "h_count", "epsilon",
        "total_reward", "final_energy", "reached_length",
        "completed_chain", "best_energy_for_protein",
        "total_opt_steps", "elapsed_s",
    ]
    with open(episode_log_path, "w", newline="") as f:
        # newline="" prevents blank rows on platforms with CSV newline handling.
        writer = csv.DictWriter(f, fieldnames=episode_fields)
        writer.writeheader()
        writer.writerows(episode_logs)

    # The summary combines training aggregates, greedy evaluation and config.
    summary_fields = [
        "protein", "length", "h_count", "episodes_seen", "best_energy",
        "avg_final_energy", "last_final_energy", "best_reward",
        "avg_reward", "last_reward", "best_length", "avg_reached_length",
        "completion_rate", "greedy_final_energy", "greedy_reached_length",
        "greedy_completed_chain", "training_time_s", "episodes",
        "batch_size", "gamma", "epsilon_start", "epsilon_end",
        "epsilon_decay", "target_update_steps", "replay_buffer_capacity",
        "learning_rate",
    ]
    summary_rows = []
    for name, stats in sequence_stats.items():
        # max(1, ...) prevents division by zero for partially initialized stats.
        n = max(1, stats["episodes_seen"])
        g = greedy_results.get(name, {})
        # Convert running sums into averages only at report-generation time.
        summary_rows.append({
            "protein":              name,
            "length":               stats["length"],
            "h_count":              stats["h_count"],
            "episodes_seen":        stats["episodes_seen"],
            "best_energy":          stats["best_energy"],
            "avg_final_energy":     round(stats["sum_energy"]  / n, 3),
            "last_final_energy":    stats["last_energy"],
            "best_reward":          round(stats["best_reward"], 3),
            "avg_reward":           round(stats["sum_reward"]  / n, 3),
            "last_reward":          round(stats["last_reward"], 3),
            "best_length":          stats["best_length"],
            "avg_reached_length":   round(stats["sum_length"]  / n, 3),
            "completion_rate":      round(stats["complete_episodes"] / n, 3),
            "greedy_final_energy":  g.get("greedy_final_energy"),
            "greedy_reached_length":g.get("greedy_reached_length"),
            "greedy_completed_chain":g.get("greedy_completed_chain"),
            "training_time_s":      training_config["training_time_s"],
            "episodes":             training_config["episodes"],
            "batch_size":           training_config["batch_size"],
            "gamma":                training_config["gamma"],
            "epsilon_start":        training_config["epsilon_start"],
            "epsilon_end":          training_config["epsilon_end"],
            "epsilon_decay":        training_config["epsilon_decay"],
            "target_update_steps":  training_config["target_update_steps"],
            "replay_buffer_capacity":training_config["replay_buffer_capacity"],
            "learning_rate":        training_config["learning_rate"],
        })

    with open(summary_path, "w", newline="") as f:
        # One row per protein supports spreadsheet and plotting workflows.
        writer = csv.DictWriter(f, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    # Build clean metric lists before aggregate NumPy calculations.
    best_vals  = [r["best_energy"]        for r in summary_rows
                  if r["best_energy"] is not None]
    avg_vals   = [r["avg_final_energy"]   for r in summary_rows]
    comp_rates = [r["completion_rate"]    for r in summary_rows]

    with open(report_path, "w") as f:
        # The text report is intentionally human-readable for thesis inspection.
        f.write("=== Constructive DQN Training Report ===\n")
        f.write("Training configuration:\n")
        for k, v in training_config.items():
            f.write(f"  {k}: {v}\n")
        f.write("\nAggregate metrics:\n")
        f.write(f"  Proteins evaluated:  {len(summary_rows)}\n")
        if best_vals:
            f.write(f"  Median best energy:  {float(np.median(best_vals)):.3f}\n")
            f.write(f"  Mean best energy:    {float(np.mean(best_vals)):.3f}\n")
            f.write(f"  Best achieved energy:{min(best_vals)}\n")
        else:
            f.write("  Median best energy:  N/A\n")
            f.write("  Mean best energy:    N/A\n")
            f.write("  Best achieved energy:N/A\n")
        if avg_vals:
            f.write(f"  Mean final energy:   {float(np.mean(avg_vals)):.3f}\n")
        else:
            f.write("  Mean final energy:   N/A\n")
        if comp_rates:
            f.write(f"  Mean completion rate:{float(np.mean(comp_rates)):.3f}\n")
        else:
            f.write("  Mean completion rate:N/A\n")
        f.write("\nPer-protein summary:\n")
        f.write("protein,length,h_count,episodes_seen,best_energy,"
                "avg_final_energy,best_reward,avg_reward,completion_rate,"
                "greedy_final_energy,greedy_reached_length\n")
        for r in summary_rows:
            f.write(f"{r['protein']},{r['length']},{r['h_count']},"
                    f"{r['episodes_seen']},{r['best_energy']},"
                    f"{r['avg_final_energy']},{r['best_reward']},"
                    f"{r['avg_reward']},{r['completion_rate']},"
                    f"{r['greedy_final_energy']},{r['greedy_reached_length']}\n")

    print(f"Episode log  → '{episode_log_path}'")
    print(f"Summary CSV  → '{summary_path}'")
    print(f"Text report  → '{report_path}'")


# ==========================================
# 8. MAIN SCRIPT
# ==========================================

PROTEINS = [
    # Fixed multi-protein dataset used to train one shared policy network.
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
    ("P35321", "MNSQQQKQPCTPPPQPQQQQVKQPCQPPPQEPCIPKTKEPCHPKVPEPCHPKVPEPCQPKVPEPCQPKVPEPCPSTVTPAPAQQKTKQK"),
]

def sequence_to_hp(sequence):
    """Reduce a one-letter amino-acid sequence to the binary HP alphabet."""

    hydrophobic = set('ACFILMVWY')
    return ''.join('H' if aa in hydrophobic else 'P' for aa in sequence)


if __name__ == "__main__":
    # Convert all raw proteins once before entering the training loop.
    training_proteins = [(name, sequence_to_hp(seq)) for name, seq in PROTEINS]

    # ── Colab-safe checkpoint settings ─────────────────────────────────────
    # If Colab disconnects, use the latest checkpoint below to resume:
    #   RESUME_FROM = "dqn_checkpoints/dqn_checkpoint_latest.pth"
    #
    # Tip: if you mount Google Drive, set CHECKPOINT_DIR to a Drive path, e.g.:
    #   CHECKPOINT_DIR = "/content/drive/MyDrive/dqn_checkpoints"
    CHECKPOINT_DIR = "/content/drive/MyDrive/dqn_checkpoints"
    # Save recoverable state every 500 completed episodes.
    CHECKPOINT_INTERVAL = 500
    # Set this to None for a fresh run; otherwise the file must already exist.
    RESUME_FROM = "/content/drive/MyDrive/dqn_checkpoints/dqn_checkpoint_latest.pth"
    SAVE_REPLAY_BUFFER = False  # True resumes more exactly but creates huge files.
    SAVE_NUMBERED_CHECKPOINTS = False  # True keeps ep00500, ep01000, ... snapshots.

    print(f"Starting Constructive DQN on {len(training_proteins)} proteins...")
    print("Improvements active:")
    print("  ✓ 3-layer CNN (7×7 receptive field)")
    print("  ✓ Replay buffer 100k")
    print("  ✓ Target network update every 1,000 steps (not episodes)")
    print("  ✓ Gradient clipping (max_norm=10)")
    print("  ✓ Double DQN")
    print("  ✓ Reward shaping for low-H-density proteins")
    print()

    trained_net, energy_hist, best_E_dict, seq_stats, ep_logs, tr_config = \
        train_dqn(
            training_proteins,
            episodes=20000,
            batch_size=128,
            checkpoint_interval=CHECKPOINT_INTERVAL,
            checkpoint_dir=CHECKPOINT_DIR,
            resume_from=RESUME_FROM,
            save_replay_buffer=SAVE_REPLAY_BUFFER,
            save_numbered_checkpoints=SAVE_NUMBERED_CHECKPOINTS,
        )

    print("\n=== Training Completed ===")
    print(f"Total optimisation steps: {tr_config['total_opt_steps']:,}")
    print(f"Training time: {tr_config['training_time_s']:.1f}s")
    print()
    print("Best energies per protein:")
    for pname, energy in best_E_dict.items():
        stats = seq_stats[pname]
        avg_e = stats["sum_energy"] / stats["episodes_seen"]
        cr    = stats["complete_episodes"] / stats["episodes_seen"]
        print(f"  {pname:8s} | Len:{stats['length']:3d} | H:{stats['h_count']:3d} | "
              f"Best E:{energy:4d} | Avg E:{avg_e:7.2f} | "
              f"Completion:{cr:.2f}")

    torch.save(trained_net.state_dict(), 'dqn_weights.pth')
    print("\nWeights saved → 'dqn_weights.pth'")

    # Greedy evaluation over the complete training set.
    print("\nRunning greedy evaluation (10 rollouts per protein)...")
    greedy_res = evaluate_training_set(trained_net, training_proteins,
                                       n_rollouts=10)
    write_dqn_reports(seq_stats, ep_logs, greedy_res, tr_config)

    # Test generalization on one HP sequence unseen during training.
    print("\n--- Test on a sequence unseen during training ---")
    test_hp = "HPHHHPPHPPPPHPHPHPPPPPHPPHPHPHPPPPPHPPPPPPHHHHPPPHPPPPPHPPHPHPPPPPHPHHHPHPPP"
    print(f"Sequence length {len(test_hp)}, {test_hp.count('H')} H residues")
    te, tl = evaluate_dqn(trained_net, test_hp, n_rollouts=10)
    print(f"Greedy result (best over 10 rollouts) → Energy: {te} | "
          f"Length: {tl}/{len(test_hp)}")

    # Plot raw episode energy plus a 100-episode moving average.
    plt.figure(figsize=(12, 5))
    plt.plot(energy_hist, alpha=0.2, color='steelblue',
             label='Energy per episode')
    window = 100
    if len(energy_hist) > window:
        smoothed = np.convolve(energy_hist,
                               np.ones(window) / window, mode='valid')
        plt.plot(range(window - 1, len(energy_hist)), smoothed,
                 color='crimson', linewidth=2,
                 label=f'Moving average (w={window})')
    plt.title('DQN — Learning Curve (Multi-Protein Dataset)')
    plt.xlabel('Episodes')
    plt.ylabel('Energy (lower is better)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('dqn_learning_curve.png', dpi=150)
    print("\nLearning curve saved → 'dqn_learning_curve.png'")
