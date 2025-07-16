from analyse_head_movements import parse_head_rotation_data, analyse_head_movements, animate_head_rotations

time_list,quat_list, = parse_head_rotation_data("1202507151601", r"C:\Users\OCONNESP\Downloads\rotation  test.txt")
print(analyse_head_movements(quat_list))
animate_head_rotations(quat_list[1])