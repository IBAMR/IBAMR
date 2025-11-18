# Read the vertex file
with open('eel2d.vertex', 'r') as f:
    lines = f.readlines()

# First line is the number of vertices
num_vertices = lines[0].strip()

# Shift values
x_shift = 2.0
y_shift = -0.5

# Process each vertex line
new_lines = [num_vertices + '\n']
for i in range(1, len(lines)):
    parts = lines[i].split()
    if len(parts) == 2:
        x = float(parts[0]) + x_shift
        y = float(parts[1]) + y_shift
        new_lines.append(f'{x:.6f}\t{y:.6f}\n')

# Write to new file
with open('eel2d_shifted.vertex', 'w') as f:
    f.writelines(new_lines)