# visualization.py

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon  # Alias to avoid conflict
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button, Slider
import numpy as np
from shapely.geometry import Polygon as ShapelyPolygon

class Visualization:
    def __init__(self, physics_engine, folding_mechanism, params):
        self.physics_engine = physics_engine
        self.folding_mechanism = folding_mechanism
        self.params = params
        self.is_paused = False
        self.dragging = False
        self.drag_index = None
        self.fig, self.ax = plt.subplots()
        plt.subplots_adjust(bottom=0.25)
        self.setup_ui()
        self.initialize_plot()

    def setup_ui(self):
        # Start/Pause Button
        ax_button = plt.axes([0.7, 0.05, 0.1, 0.04])
        self.button = Button(ax_button, 'Start/Pause')
        self.button.on_clicked(self.toggle_simulation)

        # Parameter Sliders
        ax_slider_k = plt.axes([0.15, 0.1, 0.4, 0.02])
        self.slider_k = Slider(ax_slider_k, 'Spring Constant', 0.1, 5.0, valinit=self.params.spring_constant)
        self.slider_k.on_changed(self.update_parameters)

        ax_slider_b = plt.axes([0.15, 0.06, 0.4, 0.02])
        self.slider_b = Slider(ax_slider_b, 'Damping Coeff.', 0.01, 1.0, valinit=self.params.damping_coefficient)
        self.slider_b.on_changed(self.update_parameters)

        # Connect event handlers
        self.cid_press = self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_release = self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def initialize_plot(self):
        self.ax.set_xlim(0, self.params.domain_size)
        self.ax.set_ylim(0, self.params.domain_size)
        self.ax.set_aspect('equal')
        self.polygons = []

    def toggle_simulation(self, event):
        self.is_paused = not self.is_paused

    def update_parameters(self, val):
        self.params.spring_constant = self.slider_k.val
        self.params.damping_coefficient = self.slider_b.val

    def draw_voronoi(self):
        self.ax.clear()
        voronoi = self.folding_mechanism.voronoi
        folding_states = self.folding_mechanism.folding_states
        positions = self.physics_engine.positions
        domain_size = self.params.domain_size

        # Draw polygons
        for i, region_index in enumerate(voronoi.point_region):
            region = voronoi.regions[region_index]
            if -1 in region or len(region) == 0:
                continue
            vertices = voronoi.vertices[region]
            # Apply periodic boundary corrections
            vertices = vertices % domain_size
            # Use MplPolygon for visualization
            polygon = MplPolygon(vertices,
                                 closed=True,
                                 color='red' if folding_states[i] else 'blue',
                                 alpha=0.6)
            self.ax.add_patch(polygon)

        # Draw nuclei
        self.ax.plot(positions[:, 0], positions[:, 1], 'ko', markersize=2)

        self.ax.set_xlim(0, domain_size)
        self.ax.set_ylim(0, domain_size)
        self.ax.set_aspect('equal')
        self.ax.set_title('Voronoi Origami Simulation')
        self.ax.axis('off')
        plt.draw()

    def on_press(self, event):
        if event.inaxes != self.ax:
            return
        positions = self.physics_engine.positions
        distances = np.hypot(positions[:, 0] - event.xdata, positions[:, 1] - event.ydata)
        if np.min(distances) < 0.02 * self.params.domain_size:
            self.dragging = True
            self.drag_index = np.argmin(distances)

    def on_release(self, event):
        self.dragging = False
        self.drag_index = None

    def on_motion(self, event):
        if self.dragging and event.inaxes == self.ax:
            self.physics_engine.positions[self.drag_index] = np.array([event.xdata, event.ydata]) % self.params.domain_size

    def animate(self, frame):
        if not self.is_paused:
            # Update physics
            self.physics_engine.update()
            # Update Voronoi diagram and folding states
            self.folding_mechanism.update_voronoi(self.physics_engine.positions)
            self.folding_mechanism.update_folding_states(self.physics_engine.positions)
            # Redraw
            self.draw_voronoi()

    def run(self):
        # Initialize Voronoi diagram and folding states
        self.folding_mechanism.update_voronoi(self.physics_engine.positions)
        # Start animation
        self.anim = FuncAnimation(self.fig, self.animate, interval=50, cache_frame_data=False)
        plt.show()
