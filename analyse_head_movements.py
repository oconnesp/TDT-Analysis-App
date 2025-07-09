import re
import numpy as np
from datetime import datetime, timedelta
from typing import List, Dict, Any
from txt_parsing import strip_invisibles
import re
from datetime import datetime
from pyquaternion import Quaternion
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
    #arbitrarily going to choose the fifth and fifth-last seconds as samples
    initial_head_rotations = []
    final_head_rotations = []
    for trial in quaternion_lists:
        initial_head_rotations.append(trial[4])
        final_head_rotations.append(trial[-5])
    euler_angle_list = []
    for i, (q1, q2) in enumerate(zip(initial_head_rotations, final_head_rotations)):
        delta_q = q1.inverse() * q2
        yaw, pitch, roll = delta_q.yaw_pitch_roll()
        euler_angle_list.append((yaw, pitch, roll))
    return euler_angle_list


    