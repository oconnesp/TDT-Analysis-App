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




def animate_head_rotations (quaternion_list):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    def update(i):
        ax.cla()  # Clear previous arrows and reset axis

        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        ax.set_zlim([-1, 1])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        R = quaternion_list[i].rotation_matrix
        x_axis = R[:, 0]
        y_axis = R[:, 1]
        z_axis = R[:, 2]

        ax.quiver(0, 0, 0, *x_axis, color='r', length=0.8)
        ax.quiver(0, 0, 0, *y_axis, color='g', length=0.8)
        ax.quiver(0, 0, 0, *z_axis, color='b', length=0.8)
        ax.set_title(f"Frame {i}")
    ani = FuncAnimation(fig, update, frames=len(quaternion_list), interval=500)
    plt.show()
    return ani
