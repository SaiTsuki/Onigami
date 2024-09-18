# main.py

from simulation import Simulation

if __name__ == '__main__':
    grid_size = 5  # Adjust grid size as needed
    num_agents = 10  # Initial number of agents
    max_steps = 100  # Maximum number of simulation steps

    sim = Simulation(grid_size, num_agents)
    sim.run(max_steps)
