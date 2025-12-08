import sys
import math


def is_point_in_cone(point_x, point_y, point_z, drone_x, drone_y, drone_z, 
                     direction_x, direction_y, direction_z, 
                     direction_length, cos_half_angle):
    vector_x = point_x - drone_x
    vector_y = point_y - drone_y
    vector_z = point_z - drone_z
    
    if vector_x == 0 and vector_y == 0 and vector_z == 0:
        return True
    
    dot_product = vector_x * direction_x + vector_y * direction_y + vector_z * direction_z
    
    if dot_product <= 0:
        return False
    
    vector_length_squared = vector_x * vector_x + vector_y * vector_y + vector_z * vector_z
    vector_length = math.sqrt(vector_length_squared)
    cos_angle = dot_product / (vector_length * direction_length)
    
    return cos_angle + 1e-12 >= cos_half_angle


def main():
    first_line = sys.stdin.readline().strip()
    if not first_line:
        return
    
    tokens = first_line.split()
    if len(tokens) < 7:
        return
    
    drone_x, drone_y, drone_z = int(tokens[0]), int(tokens[1]), int(tokens[2])
    direction_x, direction_y, direction_z = int(tokens[3]), int(tokens[4]), int(tokens[5])
    half_angle_degrees = int(tokens[6])
    
    second_line = sys.stdin.readline().strip()
    while not second_line:
        second_line = sys.stdin.readline().strip()
    
    num_modules = int(second_line)
    
    direction_length_squared = direction_x * direction_x + direction_y * direction_y + direction_z * direction_z
    
    if direction_length_squared == 0:
        for _ in range(num_modules):
            sys.stdin.readline()
        print(0)
        return
    
    direction_length = math.sqrt(direction_length_squared)
    cos_half_angle = math.cos(math.radians(half_angle_degrees))
    
    occupied_cells = {}
    
    for _ in range(num_modules):
        line = sys.stdin.readline().strip()
        while not line:
            line = sys.stdin.readline().strip()
        
        tokens = line.split()
        while len(tokens) < 13:
            additional_line = sys.stdin.readline().strip()
            if not additional_line:
                break
            tokens.extend(additional_line.split())
        
        if len(tokens) < 13:
            continue
        
        values = list(map(int, tokens))
        building_id = values[0]
        coordinates = values[1:]
        
        for section_idx in range(0, 12, 3):
            section_x = coordinates[section_idx]
            section_y = coordinates[section_idx + 1]
            section_z = coordinates[section_idx + 2]
            
            if is_point_in_cone(section_x, section_y, section_z,
                               drone_x, drone_y, drone_z,
                               direction_x, direction_y, direction_z,
                               direction_length, cos_half_angle):
                cell_key = (section_x, section_y, section_z)
                if cell_key not in occupied_cells:
                    occupied_cells[cell_key] = set()
                occupied_cells[cell_key].add(building_id)
    
    conflict_cells = []
    
    for (cell_x, cell_y, cell_z), building_ids in occupied_cells.items():
        if len(building_ids) >= 2:
            sorted_building_ids = sorted(building_ids)
            conflict_cells.append((cell_x, cell_y, cell_z, sorted_building_ids))
    
    conflict_cells.sort()
    
    output_lines = [str(len(conflict_cells))]
    for cell_x, cell_y, cell_z, building_ids in conflict_cells:
        line_parts = [str(cell_x), str(cell_y), str(cell_z)] + [str(bid) for bid in building_ids]
        output_lines.append(" ".join(line_parts))
    
    sys.stdout.write("\n".join(output_lines))


if __name__ == "__main__":
    main()
