# utils.py

import numpy as np

def initialize_nuclei(N, domain_size):
    # Generate random positions within the domain
    positions = np.random.uniform(0, domain_size, (N, 2))
    # Initialize velocities to zero
    velocities = np.zeros_like(positions)
    return positions, velocities
