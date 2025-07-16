import re
import numpy as np
from datetime import datetime, timedelta
from typing import List, Dict, Any
from txt_parsing import strip_invisibles
import re
from datetime import datetime
from pyquaternion import Quaternion
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import pyvista as pv
import time
import vtk

def parse_head_rotation_data(participant_id: str, filename: str):
    with open(filename, 'r', encoding='utf-8') as file:
        text=strip_invisibles(file.read()).replace("â€‹", "")

    time_lists = []
    quaternion_lists = []
    current_times = []
    current_quats = []
    inside_correct_participant = False

    # Pattern to match the desired participant header
    participant_pattern = re.compile(rf'Participant ID:\s*{re.escape(participant_id)}\s*head rotation data', re.IGNORECASE)
    # Pattern to match quaternion and timestamp entries
    rotation_pattern = re.compile(r'at (.+?):\s+\(([-\d., ]+)\)')

    for line in text.splitlines():
        if participant_pattern.search(line):
            if current_times and current_quats:
                time_lists.append(current_times)
                quaternion_lists.append(current_quats)
                current_times = []
                current_quats = []
            inside_correct_participant = True
            continue
        if inside_correct_participant:
            if line.strip().startswith("Participant ID:"):
                # Encountering a new participant ends the current one's block
                inside_correct_participant = False
                if current_times and current_quats:
                    time_lists.append(current_times)
                    quaternion_lists.append(current_quats)
                    current_times = []
                    current_quats = []
                continue
            match = rotation_pattern.search(line)
            if match:
                time_str, quat_str = match.groups()
                try:
                    timestamp = datetime.strptime(time_str.strip(), "%A, %B %d, %Y %I:%M:%S %p")
                    quaternion = Quaternion(*map(float, quat_str.split(',')))
                    current_times.append(timestamp)
                    current_quats.append(quaternion)
                except ValueError:
                    continue  # Skip malformed lines

    if current_times and current_quats:
        time_lists.append(current_times)
        quaternion_lists.append(current_quats)
    return time_lists, quaternion_lists


def analyse_head_movements(quaternion_lists):
    """Analyse head movements to calculate the difference in head rotation at the start vs at the end"""
    # arbitrarily going to choose the fifth and fifth-last seconds as samples
    initial_head_rotations = []
    final_head_rotations = []
    for trial in quaternion_lists:
        initial_head_rotations.append(trial[4])
        final_head_rotations.append(trial[-5])
    euler_angle_list = []
    for q1, q2 in zip(initial_head_rotations, final_head_rotations):
        delta_q = q1.inverse * q2
        yaw, pitch, roll = map(float, np.degrees(delta_q.yaw_pitch_roll))
        # Store as a tuple with labels for readability
        euler_angle_list.append({
            'yaw': yaw,
            'pitch': pitch,
            'roll': roll
        })
    return euler_angle_list

def quaternion_to_vtk_transform(q):
    """Convert a scipy-style quaternion into a vtkTransform."""
    # Build a 4×4 VTK matrix
    mat = vtk.vtkMatrix4x4()
    R = q.rotation_matrix  # shape (3,3)
    # copy the 3×3
    for i in range(3):
        for j in range(3):
            mat.SetElement(i, j, float(R[i, j]))
    # bottom row & right column for affine
    mat.SetElement(3, 3, 1.0)
    # wrap it in a vtkTransform
    transform = vtk.vtkTransform()
    transform.SetMatrix(mat)
    return transform


def animate_head_rotations(quaternion_list, pause=0.5):
    plotter = pv.Plotter()
    plotter.set_background("white")

    # Add a single cone actor
    cone  = pv.Cone(center=(0, 0, 0),
                    direction=(1, 0, 0),
                    height=1.0,
                    radius=0.3)
    actor = plotter.add_mesh(cone, color="lightblue")

    # Camera
    plotter.camera.position    = (5, 0, 2)
    plotter.camera.up          = (0, 0, 1)
    plotter.camera.focal_point = (0, 0, 0)

    # Initialize the interactor
    plotter.show(auto_close=False, interactive_update=True)

    for q in quaternion_list:
        # Build & apply a full 3×3 transform
        tf = quaternion_to_vtk_transform(q)
        actor.SetUserTransform(tf)

        # Push the new pose to screen
        plotter.update()
        time.sleep(pause)

    plotter.close()



