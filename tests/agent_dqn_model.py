"""
Inference utilities for the Colab-trained DQN model.

This model is the constructive DQN exported from the Colab DQN notebooks. It is
not compatible with the older DQNAgent in agent_dqn.py: the constructive model
uses a 2-channel 21x21 local grid and 3 relative actions.

This file is intentionally inference-only. Training is done offline in Colab;
the Django site imports this module to load the saved .pth weights and generate
2D HP folds quickly during a web request.
"""

import os
import random

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

try:
    # A web worker benefits from predictable CPU usage. The environment value
    # lets deployment choose a suitable thread count without changing code.
    torch.set_num_threads(int(os.environ.get("TORCH_NUM_THREADS", "1")))
except (TypeError, ValueError, RuntimeError):
    # Invalid or unsupported thread settings should not prevent inference.
    pass


class ConstructiveHPEnv:
    """Minimal constructive HP environment used only for DQN inference."""

    def __init__(self, hp_string, grid_size=21):
        """Store sequence/grid metadata and initialize the first episode."""

        # The HP sequence determines both chain length and residue channels.
        self.hp_string = hp_string
        self.n = len(hp_string)
        # The odd grid size gives a unique center cell for the chain head.
        self.grid_size = grid_size
        self.half_grid = grid_size // 2
        self.reset()

    def reset(self):
        """Reset the constructive episode and return its initial state."""

        # Start with one or two residues already placed. This gives the chain
        # an initial direction, matching the constructive policy trained in Colab.
        if self.n == 0:
            # Empty input is already terminal and has no possible energy.
            self.positions = []
            # No residue has been placed yet.
            self.current_step = 0
            # Keep a valid default direction so state creation remains defined.
            self.current_dir = (1, 0)
            # There is nothing to construct, therefore the episode is finished.
            self.done = True
            self.total_energy = 0
            return self._get_state()

        # Place residue zero at the origin. A second residue, when present,
        # establishes the initial orientation used by relative actions.
        self.positions = [(0, 0)]
        if self.n > 1:
            # The bond from residue 0 to residue 1 points along positive x.
            self.positions.append((1, 0))

        # current_step equals the number of residues successfully placed.
        self.current_step = len(self.positions)
        # Directions are unit vectors on the square lattice: this means right.
        self.current_dir = (1, 0)
        # Chains of length one or two are complete immediately after reset.
        self.done = self.current_step >= self.n
        self.total_energy = calculate_hp_energy(self.positions, self.hp_string)
        return self._get_state()

    def step(self, action):
        """Apply one relative action and return ``(state, reward, done)``."""

        if self.done:
            # Calling step after termination is harmless and changes no state.
            return self._get_state(), 0.0, True

        # Actions are relative to the current chain direction:
        # 0 = forward, 1 = left turn, 2 = right turn.
        dx, dy = self.current_dir
        if action == 0:
            # FORWARD: preserve the current unit direction vector.
            new_dir = (dx, dy)
        elif action == 1:
            # LEFT: rotate (dx, dy) by +90 degrees -> (-dy, dx).
            new_dir = (-dy, dx)
        elif action == 2:
            # RIGHT: rotate (dx, dy) by -90 degrees -> (dy, -dx).
            new_dir = (dy, -dx)
        else:
            raise ValueError(f"Unknown action index: {action}")

        new_pos = (
            # Extend the chain by one lattice unit in the selected direction.
            self.positions[-1][0] + new_dir[0],
            self.positions[-1][1] + new_dir[1],
        )

        if new_pos in self.positions:
            # A collision violates the self-avoiding walk constraint, so this
            # rollout stops immediately.
            self.done = True
            return self._get_state(), -1.0, True

        # Keep the old energy so the immediate improvement can become reward.
        before_energy = self.total_energy
        self.positions.append(new_pos)
        # Future relative turns must start from the direction just selected.
        self.current_dir = new_dir
        # Advancing this count also identifies the next HP residue to place.
        self.current_step += 1
        self.total_energy = calculate_hp_energy(self.positions, self.hp_string)

        # Lower HP energy is better: -2 -> -3 therefore gives reward +1.
        reward = float(before_energy - self.total_energy)
        if self.current_step == self.n:
            self.done = True
            # A small terminal bonus favors complete, collision-free folds.
            reward += 0.5

        return self._get_state(), reward, self.done

    def _get_state(self):
        """Encode placed residues as a head-centered, rotation-normalized grid."""

        # The state is local and rotation-invariant: the grid is centered on
        # the chain head, and the current direction is always mapped upward.
        # Axis order is (channel, row/y, column/x), as expected by Conv2d once
        # a leading batch dimension is added.
        state = np.zeros((2, self.grid_size, self.grid_size), dtype=np.float32)
        if not self.positions:
            return state

        # Every residue is represented relative to the current chain head.
        head_x, head_y = self.positions[-1]
        dx, dy = self.current_dir

        for i, (x, y) in enumerate(self.positions):
            # ``i`` links this coordinate to residue type hp_string[i].
            # Translate global coordinates into head-relative coordinates.
            rel_x = x - head_x
            rel_y = y - head_y
            # Rotate the local coordinate system so equivalent orientations
            # produce the same network input.
            rot_x = rel_x * dy - rel_y * dx
            rot_y = rel_x * dx + rel_y * dy
            # Move centered mathematical coordinates into valid array indices;
            # with size 21, local (0,0) becomes cell (10,10).
            grid_x = rot_x + self.half_grid
            grid_y = rot_y + self.half_grid

            # Very distant residues may fall outside the finite local window.
            if 0 <= grid_x < self.grid_size and 0 <= grid_y < self.grid_size:
                # Separate channels prevent H and P residues from being confused.
                channel = 0 if self.hp_string[i] == "H" else 1
                # Occupied cells are 1; untouched/empty cells remain 0.
                state[channel, grid_y, grid_x] = 1.0

        return state


class DQN_CNN(nn.Module):
    """CNN used by the final constructive DQN trained in Colab.

    Input shape is ``(batch, 2, 21, 21)`` and output shape is ``(batch, 3)``.
    Each output is the predicted Q-value of Forward, Left or Right.
    """

    def __init__(self, grid_size=21):
        """Create three spatial feature layers and a Q-value head."""

        super().__init__()
        # Padding preserves the 21x21 spatial dimensions after every 3x3 conv.
        # Shape: (batch, 2, 21, 21) -> (batch, 32, 21, 21).
        self.conv1 = nn.Conv2d(2, 32, kernel_size=3, stride=1, padding=1)
        # Shape: (batch, 32, 21, 21) -> (batch, 64, 21, 21).
        self.conv2 = nn.Conv2d(32, 64, kernel_size=3, stride=1, padding=1)
        # Shape remains (batch, 64, 21, 21) while features are refined.
        self.conv3 = nn.Conv2d(64, 64, kernel_size=3, stride=1, padding=1)
        # No pooling is used, so all 64 feature maps retain the full grid size.
        self.flat_size = 64 * grid_size * grid_size #matrix converted in a vector of size 64*21*21 = 28224
        self.fc1 = nn.Linear(self.flat_size, 256) #fully connected layer that takes the flattened output of the convolutional layers and maps it to a 256-dimensional hidden representation
        self.fc2 = nn.Linear(256, 3) #final fully connected layer that maps the hidden representation to three output values (Forward, Left, Right).

    def forward(self, x):
        """Map a batch of state grids to three unnormalized Q-values."""

        # ReLU introduces non-linearity after every learned hidden layer.
        x = F.relu(self.conv1(x)) #Relu activation function is applied to the output of the first convolutional layer, introducing non-linearity and allowing the network to learn complex patterns in the input data.
        x = F.relu(self.conv2(x))
        x = F.relu(self.conv3(x))
        # Keep the first dimension as batch size and flatten everything else.
        x = x.view(x.size(0), -1) #The output of the last convolutional layer is flattened into a 1D vector for each sample in the batch, preparing it for the fully connected layers.
        x = F.relu(self.fc1(x))
        # The output is intentionally linear: Q-values are not probabilities.
        return self.fc2(x) #The final fully connected layer produces three unnormalized Q-values corresponding to the three possible actions (Forward, Left, Right).


class LegacyDQN_CNN(nn.Module):
    """Older CNN kept only as a compatibility fallback."""

    def __init__(self, grid_size=21):
        """Recreate the smaller architecture expected by older checkpoints."""

        super().__init__()
        # This architecture must exactly match old layer names and tensor sizes.
        self.conv1 = nn.Conv2d(2, 16, kernel_size=3, stride=1, padding=1)
        self.conv2 = nn.Conv2d(16, 32, kernel_size=3, stride=1, padding=1)
        self.flat_size = 32 * grid_size * grid_size
        self.fc1 = nn.Linear(self.flat_size, 128)
        self.fc2 = nn.Linear(128, 3)

    def forward(self, x):
        """Compute legacy checkpoint Q-values for the three actions."""

        x = F.relu(self.conv1(x))
        x = F.relu(self.conv2(x))
        x = x.view(x.size(0), -1)
        x = F.relu(self.fc1(x))
        return self.fc2(x)


# Process-local cache: one Django worker loads each checkpoint only once.
_MODEL_CACHE = {}


def _torch_load_weights(weights_path, device):
    """Load a checkpoint safely across multiple supported PyTorch versions."""

    # weights_only=True is available in newer PyTorch versions. The fallback
    # keeps the loader compatible with older environments.
    try:
        # map_location allows GPU-trained checkpoints to load on a CPU server.
        return torch.load(weights_path, map_location=device, weights_only=True)
    except TypeError:
        return torch.load(weights_path, map_location=device)


def _extract_state_dict(weights_obj):
    """Extract raw network parameters from common checkpoint wrappers."""

    # Training scripts may save either a bare state_dict or a larger dictionary
    # containing optimizer/logging data under one of these conventional keys.
    if isinstance(weights_obj, dict):
        for key in ("policy_net_state_dict", "model_state_dict", "state_dict"):
            # Stop at the first recognized wrapper containing a parameter map.
            if key in weights_obj and isinstance(weights_obj[key], dict):
                return weights_obj[key]
    return weights_obj


def _build_policy_net_for_weights(state_dict, device):
    """Select the architecture whose layer names match the checkpoint."""

    # Final checkpoints contain conv3; its absence identifies the legacy CNN.
    if isinstance(state_dict, dict) and "conv3.weight" not in state_dict:
        return LegacyDQN_CNN(grid_size=21).to(device)
    return DQN_CNN(grid_size=21).to(device)


def load_dqn_policy(weights_path):
    """Load and cache the trained DQN policy network from disk."""
    # Normalize the path so different relative spellings share one cache entry.
    abs_path = os.path.abspath(weights_path)
    try:
        stat = os.stat(abs_path)
        # Modification time and size invalidate the cache when weights change.
        cache_key = (abs_path, stat.st_mtime_ns, stat.st_size)
    except OSError:
        # Keep a deterministic key; the caller will handle load failure later.
        cache_key = (abs_path, None, None)

    if cache_key in _MODEL_CACHE:
        # Return both network and device because rollout tensors need the same
        # CPU/GPU placement as the model parameters.
        return _MODEL_CACHE[cache_key]

    # Recreate the exact architecture, load the Colab-trained weights, and set
    # the model to evaluation mode because Django only performs inference.
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # A state_dict maps layer names such as conv1.weight to numeric tensors.
    state_dict = _extract_state_dict(_torch_load_weights(weights_path, device))
    policy_net = _build_policy_net_for_weights(state_dict, device)
    # Strict loading catches architecture/checkpoint mismatches immediately.
    policy_net.load_state_dict(state_dict)
    policy_net.eval()

    stale_keys = [
        # Select older cached versions of this same checkpoint path.
        key for key in _MODEL_CACHE
        if (key[0] if isinstance(key, tuple) else key) == abs_path
    ]
    # If the .pth file changed, remove older cached versions of the same path.
    for key in stale_keys:
        del _MODEL_CACHE[key]

    _MODEL_CACHE[cache_key] = (policy_net, device)
    # Inputs must later be allocated on the same device as the network.
    return policy_net, device


def calculate_hp_energy(positions, hp_string):
    """Compute standard 2D HP energy from non-bonded adjacent H-H contacts."""
    energy = 0
    # Partial constructive rollouts contain fewer positions than residues.
    limit = min(len(positions), len(hp_string))

    for i in range(limit):
        # Polar residues do not contribute to this simplified energy function.
        if hp_string[i] != "H":
            continue
        x1, y1 = positions[i]
        # Start at i+2 to exclude covalently adjacent backbone residues.
        for j in range(i + 2, limit):
            if hp_string[j] != "H":
                continue
            x2, y2 = positions[j]
            # Manhattan distance one means orthogonal lattice adjacency.
            if abs(x1 - x2) + abs(y1 - y2) == 1:
                energy -= 1

    return energy


def straight_line_fallback(hp_string):
    """Return a valid unoptimized fold when inference is unavailable."""

    # Used when the model weights are missing or cannot be loaded. This keeps
    # the web request from crashing, even though the fold is not optimized.
    positions = [(i, 0) for i in range(len(hp_string))]
    # A line is connected and self-avoiding, but normally has no non-bonded
    # contacts and is therefore only safe, not optimized.
    return positions, calculate_hp_energy(positions, hp_string)


def _rank_candidate(candidate):
    """Sort incomplete candidates by longest chain, then lowest energy."""

    # In failed rollouts, prefer the one that placed more residues; if tied,
    # prefer the lower HP energy.
    positions, energy = candidate[:2]
    # min() prefers smaller tuples: negative length favors longer partial folds,
    # then the second tuple element favors lower (more negative) energy.
    return (-len(positions), energy)


def _run_constructive_rollout(policy_net, device, hp_string, epsilon,
                              return_trace=False, trace_interval=1,
                              max_trace_frames=200):
    """Run one complete DQN-guided construction attempt."""
    env = ConstructiveHPEnv(hp_string, grid_size=21)
    # Construction starts from the same two-residue orientation used in training.
    state = env.reset()
    trace_frames = []

    def add_trace_frame(step):
        """Snapshot the partial fold for the website's step animation."""

        if not return_trace:
            return
        if len(trace_frames) >= max_trace_frames:
            return
        trace_frames.append({
            "step": step,
            "energy": calculate_hp_energy(env.positions, hp_string),
            # Copy coordinates so later environment updates cannot alter history.
            "positions": list(env.positions),
        })

    trace_interval = max(1, int(trace_interval))
    max_trace_frames = max(2, int(max_trace_frames))
    add_trace_frame(len(env.positions))

    while not env.done:
        # A small epsilon keeps some rollouts stochastic. This can escape
        # deterministic mistakes made by the greedy policy.
        if random.random() < epsilon:
            # Exploration samples uniformly from action IDs 0, 1 and 2.
            action = random.randrange(3)
        else:
            # inference_mode disables gradients and version tracking, reducing
            # latency and memory usage during web inference.
            with torch.inference_mode():
                #converts numpy array to a PyTorch tensor
                state_tensor = torch.as_tensor(
                    state, dtype=torch.float32, device=device
                ).unsqueeze(0) #Add a batch dimension to the state tensor, making its shape (1, 2, 21, 21) for processing by the CNN.
                # Greedy action = index of the largest predicted Q-value.
                action = policy_net(state_tensor).argmax(dim=1).item()

        # Reward and returned done are ignored during inference: no Bellman
        # update occurs, and the while condition reads env.done directly.
        state, _, _ = env.step(action)
        # Trace frames are consumed by the site's step explorer animation.
        if len(env.positions) % trace_interval == 0:
            add_trace_frame(len(env.positions))

    # Collision may terminate early, so evaluate exactly the placed prefix.
    final_energy = calculate_hp_energy(env.positions, hp_string)
    if return_trace and (not trace_frames or trace_frames[-1]["step"] != len(env.positions)):
        add_trace_frame(len(env.positions))

    if return_trace:
        # Django expects a third value when step visualization is requested.
        return list(env.positions), final_energy, trace_frames
    # Non-visual callers keep the simpler (positions, energy) interface.
    return list(env.positions), final_energy


def run_dqn_inference(hp_string, weights_path, max_rollouts=64,
                      return_trace=False, trace_interval=1,
                      max_trace_frames=200):
    """
    Fold an HP string using the constructive DQN weights.

    max_rollouts controls how many greedy/stochastic constructions are sampled.
    The returned structure prefers complete folds first, then lower energy.
    """
    if not hp_string:
        # Preserve the expected return tuple shape even for empty input.
        if return_trace:
            return [], 0, []
        return [], 0

    if not os.path.exists(weights_path):
        # The model file is optional for robustness: without it, the web app
        # still returns a valid straight-line structure instead of an error.
        positions, energy = straight_line_fallback(hp_string)
        if return_trace:
            # Produce progressive synthetic frames so the UI remains usable
            # even without a neural-network checkpoint.
            trace_frames = [
                {"step": idx + 1, "energy": calculate_hp_energy(positions[:idx + 1], hp_string), "positions": positions[:idx + 1]}
                for idx in range(0, len(positions), max(1, int(trace_interval)))
            ]
            if not trace_frames or trace_frames[-1]["step"] != len(positions):
                trace_frames.append({"step": len(positions), "energy": energy, "positions": positions})
            return positions, energy, trace_frames[:max_trace_frames]
        return positions, energy

    try:
        policy_net, device = load_dqn_policy(weights_path)
    except Exception as exc:
        # Log technical details server-side and return a safe fold to the user.
        print(f"Failed to load DQN weights: {exc}")
        positions, energy = straight_line_fallback(hp_string)
        if return_trace:
            return positions, energy, [{"step": len(positions), "energy": energy, "positions": positions}]
        return positions, energy

    # Protect the web request from zero, negative or excessive rollout counts.
    rollout_count = max(1, min(int(max_rollouts), 300))
    # Try a mix of greedy and slightly stochastic constructions. More rollouts
    # increase the chance of a complete, lower-energy fold but cost more time.
    epsilons = [0.0, 0.05, 0.10, 0.20]
    # Entries are (positions, energy), plus trace_frames when requested.
    candidates = []

    for idx in range(rollout_count):
        # Modulo cycles repeatedly through all four exploration levels.
        epsilon = epsilons[idx % len(epsilons)]
        candidates.append(_run_constructive_rollout(
            policy_net,
            device,
            hp_string,
            epsilon,
            return_trace=return_trace,
            trace_interval=trace_interval,
            max_trace_frames=max_trace_frames,
        ))

    # A collision can terminate a rollout early, so completion is considered
    # more important than energy when ranking candidates.
    complete = [candidate for candidate in candidates if len(candidate[0]) == len(hp_string)]
    # Prefer complete folds first. Among complete folds, choose the lowest
    # energy. If no rollout completes, choose the longest partial structure.
    # candidate[1] is energy, so min() selects the most negative complete fold.
    # If none completed, _rank_candidate prioritizes the longest placed prefix.
    selected = min(complete, key=lambda candidate: candidate[1]) if complete else min(candidates, key=_rank_candidate)

    if return_trace:
        positions, energy, trace_frames = selected
        return positions, energy, trace_frames
    return selected
