# Onigami
# **Technical Specification (TK): Voronoi Origami Simulation**

---

## **1. Overview**

The Voronoi Origami Simulation is a dynamic Python application that models a 2D system of interacting particles (nuclei) whose positions evolve over time under physical forces. The system generates a Voronoi diagram from these nuclei, simulates physical interactions using spring dynamics and Brownian motion, and incorporates a folding mechanism where Voronoi cells change state when certain conditions are met. The simulation includes an interactive user interface with controls and real-time visualization.

---

## **2. Technical Stack**

- **Programming Language**: Python 3.x
- **Primary Libraries**:
  - **NumPy**: Numerical computations and array operations.
  - **SciPy**: Computing Voronoi and Delaunay diagrams.
  - **Matplotlib**: Visualization and plotting, including interactive features.
  - **Shapely**: Handling geometric objects and spatial queries.
  - **Numba** (Optional): Just-In-Time compilation for performance optimization.
- **Additional Libraries**:
  - **Matplotlib Widgets**: For UI elements like buttons.
  - **Tkinter** or **PyQt5** (Optional): For more advanced GUI features, if required.

---

## **3. System Functionality**

### **3.1. Initial Setup**

#### **3.1.1. Simulation Domain**

- **Domain**: Define a square domain ranging from \([-1, -1]\) to \([1, 1]\).
- **Boundary Conditions**: Implement periodic boundary conditions to ensure seamless interaction at the edges.

#### **3.1.2. Nuclei Generation**

- **Number of Nuclei**: \( N = 2048 \).
- **Position Initialization**: Generate random positions within the domain using a uniform distribution.

```python
import numpy as np

# Number of nuclei
N = 2048

# Generate random positions
nuclei_positions = np.random.uniform(-1, 1, (N, 2))
```

---

### **3.2. Voronoi and Delaunay Tessellations**

#### **3.2.1. Voronoi Diagram**

- **Computation**: Use `scipy.spatial.Voronoi` to compute the Voronoi diagram.
- **Boundary Handling**: Since SciPy does not support periodic boundaries, manually handle edge cases or consider alternative libraries if necessary.

#### **3.2.2. Delaunay Triangulation**

- **Computation**: Use `scipy.spatial.Delaunay` to compute the Delaunay triangulation, which is the dual of the Voronoi diagram.
- **Neighbor Identification**: Extract neighbor relationships from the Delaunay triangulation for force calculations.

```python
from scipy.spatial import Voronoi, Delaunay

# Compute Delaunay triangulation
delaunay = Delaunay(nuclei_positions)

# Compute Voronoi diagram
voronoi = Voronoi(nuclei_positions)
```

---

### **3.3. Spring Dynamics Between Nuclei**

#### **3.3.1. Physical Model**

- **Springs**: Model the connections between neighboring nuclei as springs with damping.
- **Neighbors**: Only consider immediate neighbors identified from the Delaunay triangulation.

#### **3.3.2. Force Calculations**

- **Hooke's Law**:
  \[
  \mathbf{F}_{\text{spring}} = -k (r - r_0) \hat{\mathbf{u}}
  \]
  - \( k \): Spring constant.
  - \( r \): Current distance between nuclei.
  - \( r_0 \): Rest length (initial distance).
  - \( \hat{\mathbf{u}} \): Unit vector from one nucleus to another.

- **Damping Force**:
  \[
  \mathbf{F}_{\text{damping}} = -b (\mathbf{v}_i - \mathbf{v}_j)
  \]
  - \( b \): Damping coefficient.
  - \( \mathbf{v}_i, \mathbf{v}_j \): Velocities of the nuclei.

#### **3.3.3. Implementation Steps**

1. **Neighbor Extraction**:
   - Build a neighbor list for each nucleus from the Delaunay triangulation.
   - Use an adjacency list or set for efficient lookup.

2. **Rest Length Calculation**:
   - Compute and store the rest lengths between each pair of neighboring nuclei.

3. **Force Computation**:
   - For each nucleus, calculate the net force from all connected springs and damping effects.
   - Use vectorized operations where possible.

4. **Position and Velocity Updates**:
   - Update velocities and positions using an integration method (e.g., Euler or Verlet).
   - Apply periodic boundary conditions after each update.

```python
# Initialize velocities
velocities = np.zeros_like(nuclei_positions)

# Spring constant and damping coefficient
k = 1.0
b = 0.1
dt = 0.01  # Time step

# Extract neighbors from Delaunay triangulation
neighbors = [set() for _ in range(N)]
for simplex in delaunay.simplices:
    for idx in simplex:
        others = simplex[simplex != idx]
        neighbors[idx].update(others)

# Compute rest lengths
rest_lengths = {}
for i, nbrs in enumerate(neighbors):
    for j in nbrs:
        key = tuple(sorted((i, j)))
        if key not in rest_lengths:
            delta = nuclei_positions[i] - nuclei_positions[j]
            rest_lengths[key] = np.linalg.norm(delta)

# Force computation function
def compute_forces(nuclei_positions, velocities):
    forces = np.zeros_like(nuclei_positions)
    for i, nbrs in enumerate(neighbors):
        for j in nbrs:
            key = tuple(sorted((i, j)))
            delta_pos = nuclei_positions[i] - nuclei_positions[j]
            delta_vel = velocities[i] - velocities[j]
            distance = np.linalg.norm(delta_pos)
            if distance == 0:
                continue  # Avoid division by zero
            direction = delta_pos / distance
            # Spring force
            rest_length = rest_lengths[key]
            spring_force = -k * (distance - rest_length) * direction
            # Damping force
            damping_force = -b * delta_vel
            # Total force
            forces[i] += spring_force + damping_force
    return forces
```

---

### **3.4. Brownian Motion**

#### **3.4.1. Implementation**

- **Random Forces**: Apply small random forces to simulate Brownian motion.
- **Noise Scale**: Define a parameter to control the magnitude of these random forces.

```python
noise_scale = 0.01

def apply_brownian_motion(forces):
    random_forces = np.random.normal(scale=noise_scale, size=forces.shape)
    forces += random_forces
    return forces
```

---

### **3.5. Folding Mechanism**

#### **3.5.1. Folding Logic**

- **Event Trigger**: A folding event occurs when a nucleus exits its own Voronoi cell.
- **State Management**: Maintain a folding state for each Voronoi cell (e.g., folded or unfolded).
- **Visual Representation**: Change the color of the Voronoi cell to reflect its folding state.

#### **3.5.2. Implementation Steps**

1. **Voronoi Cell Extraction**:
   - For each nucleus, extract its corresponding Voronoi cell as a polygon.

2. **Point-in-Polygon Test**:
   - Use `Shapely` to check if the nucleus is inside its Voronoi cell.
   - If not, toggle the folding state.

3. **State Update Function**:
   - Create a function to update the folding states after each simulation step.

```python
from shapely.geometry import Point, Polygon

# Initialize folding states
folding_states = [False] * N

def update_folding_states(nuclei_positions, voronoi):
    for i, pos in enumerate(nuclei_positions):
        region_index = voronoi.point_region[i]
        region = voronoi.regions[region_index]
        if -1 in region or len(region) == 0:
            continue  # Skip regions with infinite vertices
        vertices = voronoi.vertices[region]
        cell_polygon = Polygon(vertices)
        nucleus_point = Point(pos)
        if not cell_polygon.contains(nucleus_point):
            folding_states[i] = not folding_states[i]  # Toggle state
```

---

### **3.6. Visualization with Matplotlib**

#### **3.6.1. Plotting Voronoi Cells**

- **Drawing Cells**: For each Voronoi cell, draw a polygon with a color based on its folding state.
- **Axes Configuration**: Set the plot limits to match the domain and enable equal aspect ratio.

#### **3.6.2. Interactive Controls**

- **Start/Pause Button**: Implement a button to start or pause the simulation using Matplotlib widgets.
- **Animation**: Use `FuncAnimation` to update the visualization in real-time.

```python
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.2)  # Adjust to make room for buttons

# Function to draw the Voronoi diagram
def draw_voronoi():
    ax.clear()
    for i, region_index in enumerate(voronoi.point_region):
        region = voronoi.regions[region_index]
        if -1 in region or len(region) == 0:
            continue
        vertices = voronoi.vertices[region]
        polygon = Polygon(vertices, True,
                          color='red' if folding_states[i] else 'blue')
        ax.add_patch(polygon)
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_aspect('equal')
    plt.draw()

# Animation function
def animate(frame):
    global nuclei_positions, velocities
    if is_paused:
        return
    # Compute forces
    forces = compute_forces(nuclei_positions, velocities)
    # Apply Brownian motion
    forces = apply_brownian_motion(forces)
    # Update positions and velocities
    velocities += forces * dt
    nuclei_positions += velocities * dt
    # Apply periodic boundaries
    nuclei_positions = (nuclei_positions + 1) % 2 - 1
    # Recompute Voronoi diagram
    voronoi = Voronoi(nuclei_positions)
    # Update folding states
    update_folding_states(nuclei_positions, voronoi)
    # Redraw
    draw_voronoi()

# Start/pause functionality
is_paused = False

def toggle_simulation(event):
    global is_paused
    is_paused = not is_paused

ax_button = plt.axes([0.8, 0.05, 0.1, 0.075])
button = Button(ax_button, 'Start/Pause')
button.on_clicked(toggle_simulation)

# Initialize plot
draw_voronoi()

# Start animation
anim = FuncAnimation(fig, animate, interval=50)
plt.show()
```

---

### **3.7. Real-Time Updates and Interaction**

#### **3.7.1. User Interaction**

- **Mouse Events**: Optionally, implement mouse event handlers to allow users to interact with the nuclei (e.g., drag and drop).
- **Matplotlib Event Handling**: Use Matplotlib's event system to capture mouse clicks and movements.

#### **3.7.2. Implementation**

```python
def on_press(event):
    # Logic to select and drag nuclei
    pass

def on_release(event):
    # Logic to release nuclei
    pass

def on_motion(event):
    # Logic to move selected nuclei
    pass

fig.canvas.mpl_connect('button_press_event', on_press)
fig.canvas.mpl_connect('button_release_event', on_release)
fig.canvas.mpl_connect('motion_notify_event', on_motion)
```

---

### **3.8. Boundary Conditions**

#### **3.8.1. Periodic Boundary Conditions**

- **Position Wrapping**: Adjust positions to wrap around when they cross the domain boundaries.
- **Distance Calculations**: Modify distance calculations to account for the shortest distance considering the periodicity.

#### **3.8.2. Implementation**

```python
def apply_periodic_boundary(positions):
    return (positions + 1) % 2 - 1  # Wrap positions to [-1, 1]

def compute_periodic_distance(delta):
    delta = delta - np.round(delta / 2) * 2  # Adjust delta for periodicity
    return delta
```

---

## **4. Computational Optimization**

### **4.1. Complexity Reduction**

#### **4.1.1. Neighbor-Based Calculations**

- **Limiting Interactions**: Only compute forces between neighboring nuclei to reduce computational complexity from O(N²) to O(N).

### **4.2. Vectorization with NumPy**

#### **4.2.1. Efficient Operations**

- **Batch Calculations**: Use NumPy arrays to perform operations on all nuclei simultaneously.
- **Avoid Loops**: Replace explicit Python loops with vectorized operations where possible.

### **4.3. Just-In-Time Compilation with Numba**

#### **4.3.1. Performance Enhancement**

- **Critical Functions**: Apply `@njit` decorator to functions like `compute_forces` to speed up execution.
- **Installation**: Ensure Numba is installed in the environment.

```python
from numba import njit

@njit
def compute_forces_numba(...):
    # Optimized force computation
    pass
```

### **4.4. Efficient Voronoi Updates**

#### **4.4.1. Incremental Updates**

- **Partial Recomputations**: Investigate algorithms or libraries that support incremental updates to the Voronoi diagram to avoid full recomputation each time.

#### **4.4.2. Alternative Libraries**

- **PyVoronoi or Others**: Consider using libraries optimized for dynamic Voronoi diagrams.

---

## **5. Testing and Validation**

### **5.1. Unit Testing**

#### **5.1.1. Physics Calculations**

- **Force Functions**: Write tests to verify the correctness of force calculations.
- **Boundary Conditions**: Test position wrapping and distance calculations under periodic boundaries.

### **5.2. Integration Testing**

#### **5.2.1. Folding Mechanism**

- **Folding Logic**: Simulate nuclei moving in and out of cells and check if the folding state updates correctly.

### **5.3. Visualization Checks**

#### **5.3.1. Visual Validation**

- **Consistency**: Ensure that the visual representation matches the internal state (e.g., colors reflect folding states).
- **Smoothness**: Verify that the animation runs smoothly without significant lag.

---

## **6. Code Documentation and Maintainability**

### **6.1. Code Comments and Docstrings**

- **Function Documentation**: Provide clear docstrings for all functions, explaining their purpose, inputs, and outputs.
- **Inline Comments**: Add comments within functions to explain complex logic or important steps.

### **6.2. Modular Code Structure**

- **Organization**: Separate the code into modules (e.g., `physics.py`, `visualization.py`, `main.py`) for better maintainability.
- **Reusability**: Write reusable functions and classes where appropriate.

### **6.3. Version Control**

- **Git Usage**: Use Git for version control, allowing for tracking changes and collaboration.
- **Commit Messages**: Write clear and descriptive commit messages.

---

## **7. Additional Recommendations**

### **7.1. Parameter Tuning**

- **Experimentation**: Adjust parameters like the spring constant \( k \), damping coefficient \( b \), noise scale, and time step \( dt \) to achieve desired simulation behavior.
- **Stability**: Ensure that parameter values lead to a stable simulation (e.g., avoid numerical instability).

### **7.2. Performance Profiling**

- **Profiling Tools**: Use tools like `cProfile` or `line_profiler` to identify performance bottlenecks.
- **Optimization Focus**: Prioritize optimizing functions that consume the most time.

### **7.3. Advanced Features**

#### **7.3.1. Improved Integration Methods**

- **Verlet Integration**: Implement Verlet or Velocity-Verlet integration for better numerical stability in physics simulations.

#### **7.3.2. User Interactivity**

- **Custom Nuclei Placement**: Allow users to add or remove nuclei during the simulation.
- **Parameter Controls**: Provide sliders or input fields to adjust simulation parameters in real-time.

#### **7.3.3. Data Saving and Loading**

- **State Serialization**: Implement functionality to save the simulation state to a file and load it later.

---

## **8. Conclusion**

By following this comprehensive technical specification, the Voronoi Origami Simulation will be accurate, efficient, and interactive. The corrected folding logic ensures that the simulation behaves realistically, and optimizations allow it to handle a large number of nuclei smoothly. Proper documentation and modular code structure facilitate maintenance and future enhancements.

---

# **Appendix**

## **A. Dependencies Installation**

- **NumPy**: `pip install numpy`
- **SciPy**: `pip install scipy`
- **Matplotlib**: `pip install matplotlib`
- **Shapely**: `pip install shapely`
- **Numba**: `pip install numba` (if using)
- **PyQt5**: `pip install pyqt5` (if using for advanced GUI)

## **B. References**

- **NumPy Documentation**: [https://numpy.org/doc/](https://numpy.org/doc/)
- **SciPy Spatial**: [https://docs.scipy.org/doc/scipy/reference/spatial.html](https://docs.scipy.org/doc/scipy/reference/spatial.html)
- **Matplotlib**: [https://matplotlib.org/](https://matplotlib.org/)
- **Shapely**: [https://shapely.readthedocs.io/en/stable/](https://shapely.readthedocs.io/en/stable/)
- **Numba**: [http://numba.pydata.org/](http://numba.pydata.org/)

---

By implementing the simulation according to this specification, you will create a robust and efficient Voronoi Origami Simulation that provides valuable insights into dynamic systems and computational geometry.
