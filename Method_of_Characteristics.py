import numpy as np
from math import asin, tan
from tablegenerator import interpolate
import pandas as pd
import matplotlib.pyplot as plt
import csv
import matplotlib.tri as mtri
from matplotlib.path import Path

def main():
    A, N, x = get_matrix()
    Mach_number_at_end = float(input("What is Mdes?"))
    exp_start = float(input("What is the start angle of expansion corner?"))
    A = first_step(A, Mach_number_at_end, N, x, exp_start)
    A = second_step(A, x)
    A = third_step(A, x)
    A = fourth_step(A, x)
    A = fifth_step(A, x)
    A = sixth_step(A, N, x)
    A = seventh_step(A, N, x)
    np.savetxt("MoC.csv", A, fmt='%10.5f', delimiter = ',')
    plot_moc_data_specific_lines_v2('MoC.csv', x)
    
def first_step(A, M, N, x, exp_start):
    nu = interpolate(M, 'Mach', 'v')
    A[N-1][4] = M
    A[N-2][4] = M
    A[N-1][3] = nu
    A[N-2][3] = nu
    for i in range(x):
        A[i][8] = 0.0
        A[i][9] = 1.0
    thetamax = nu/2
    steps = (thetamax-exp_start)/(x-1)
    print(steps)
    while thetamax >= exp_start:
        A[x-1][2] = thetamax
        A[x-1][3] = thetamax
        thetamax -= steps
        x -= 1
    A[0][2] = exp_start
    A[0][3] = exp_start
    print("Step 1/7 completed")
    return A

def second_step(A, x):
    for i in range(x):
        A[i][1] = A[i][2] + A[i][3]
    print("Step 2/7 completed")
    return A

def third_step(A, x):
    for i in range(x):
        A[i][4] = interpolate(A[i][3], 'v', 'Mach')
        A[i][5] = np.degrees(asin(1/A[i][4]))
        A[i][6] = A[i][2] + A[i][5]
        A[i][7] = A[i][2] - A[i][5]
    print("Step 3/7 completed")
    return A

def fourth_step(A, x):
    lines = x
    pos = x
    times = 0
    while lines > 0:
        for i in range(times, x):
            A[pos+i][1] = A[i][2] + A[i][3]
        pos += lines
        lines -= 1
        times +=1
    print("Step 4/7 completed")
    return A

def fifth_step(A, x):
    lines = x
    pos = x
    times = 0
    while lines > 0:
        for i in range(times, x+1):
            A[pos+i][0] = A[times][1]
        pos += lines
        lines -= 1
        times += 1
    print("Step 5/7 completed")
    return A

def sixth_step(A, N, x):
    for i in range(x, N):
        if A[i][1] == 0:
            A[i][2] = A[i-1][2]
            A[i][3] = A[i-1][3]
        else:
            A[i][2] = (A[i][1]-A[i][0])/2
            A[i][3] = (A[i][1]+A[i][0])/2
        A[i][4] = interpolate(A[i][3], 'v', 'Mach')
        A[i][5] = np.degrees(asin(1/A[i][4]))
        A[i][6] = A[i][2] + A[i][5]
        A[i][7] = A[i][2] - A[i][5]
    print("Step 6/7 completed")   
    return A
    
def seventh_step(A, N, x):
    pos = x
    step = x
    delta = x
    linedifference = x+1
    factor = 1
    line = 0
    while step > 0:
        for i in range(pos+1, pos+step):
            A[pos][8] = A[pos-delta][8] -A[pos-delta][9]/tan(np.radians(0.5*(A[pos-delta][7]+A[pos][7])))
            A[pos][9] = 0
            A[i][8] = (A[i-delta][8]* tan(np.radians(0.5*(A[i-delta][7]+A[i][7])))- A[i-1][8]*tan(np.radians(0.5*(A[i-1][6]+A[i][6]))) -\
                      (A[i-delta][9]-A[i-1][9]))/(tan(np.radians(0.5*(A[i-delta][7]+A[i][7])))-tan(np.radians(0.5*(A[i-1][6]+A[i][6]))))
            A[i][9] = A[i-1][9] + (A[i][8]-A[i-1][8])*tan(np.radians(0.5*(A[i-1][6]+A[i][6])))
            A[pos+step][8] = (A[pos+step - linedifference][8]*tan(np.radians(0.5*(A[pos+step - linedifference][2]+A[pos+step][2])))-\
                              A[pos+step-1][8]*tan(np.radians(A[pos+step-1][6]))-(A[pos+step - linedifference][9]-A[pos+step-1][9]))/\
                (tan(np.radians(0.5*(A[pos+step - linedifference][2]+A[pos+step][2])))-tan(np.radians(A[pos+step-1][6])))
            A[pos+step][9] = A[pos+step-1][9] + (A[pos+step][8]-A[pos+step-1][8])*tan(np.radians(A[pos+step-1][6]))
        if delta == x and factor == 1:
            factor -= 1 
        else:
            delta -= 1
        linedifference -= 1
        line += 1
        step -= 1 
        pos += x+1 
        x -= 1
    A[N-2][8] = A[N-4][8]- A[N-4][9]/tan(np.radians(0.5*(A[N-4][7]+A[N-2][7])))
    A[N-2][9] = 0
    A[N-1][8] = (A[N-3][8]*tan(np.radians(0.5*(A[N-3][2]+A[N-1][2])))-A[N-2][8]*\
                 tan(np.radians(A[N-2][6]))-(A[N-3][9]-A[N-2][9]))/(tan(np.radians(0.5*(A[N-3][2]+A[N-1][2])))-tan(np.radians(A[N-2][6])))
    A[N-1][9] = A[N-2][9] + (A[N-1][8]-A[N-2][8])*tan(np.radians(A[N-2][6]))
    print("Step 7/7 completed")
    return A


def plot_moc_data_specific_lines_v2(file_path="MoC.csv", num_initial_lines=5 ):
    """
    Plots the Method of Characteristics nozzle design with smooth,
    filled contours for the Mach number.
    """
    try:
        df = pd.read_csv(file_path, header=None)
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found. Please ensure your main script is run first to generate it.")
        return
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return

    # Column indices from your MoC.csv
    x_col = 8
    y_col = 9
    M_col = 4

    # Extract data as numpy arrays
    x = df.iloc[:, x_col].to_numpy()
    y = df.iloc[:, y_col].to_numpy()
    M = df.iloc[:, M_col].to_numpy()
    
    # Check if we have the number of lines, which is crucial for plotting
    if num_initial_lines is None:
        # Try to infer from data if not provided (optional, but good practice)
        num_initial_lines = len(df[df[x_col] == 0])
        if num_initial_lines == 0:
             print("Error: 'num_initial_lines' was not provided and could not be inferred.")
             return
        print(f"Inferred num_initial_lines = {num_initial_lines}")

    plt.figure(figsize=(12, 7))
    ax = plt.gca()

    # --- 1. Define the Nozzle Boundaries (Wall and Centerline) ---
    # This logic is derived from your `get_matrix` and `seventh_step` functions
    # to find the indices of all wall and centerline points.
    
    wall_points_idx = []
    center_points_idx = []
    
    n = num_initial_lines
    
    # Add initial fan points (all are at (0, 1))
    # We only need the first and last for the boundary polygon
    wall_points_idx.append(0) # (0, 1)
    
    # Loop variables to trace point indices
    n_local = n
    pos = n
    
    # This loop finds the first (centerline) and last (wall)
    # point of each characteristic segment.
    while n_local > 0:
        center_points_idx.append(pos)           # e.g., 5, 11, 16...
        wall_points_idx.append(pos + n_local)   # e.g., 10, 15, 19...
        
        pos += n_local + 1
        n_local -= 1
        
    # Add the final two points (from your `seventh_step` end-calc)
    N = len(x)
    if N > 0:
        # N-2 is the final centerline point
        # N-1 is the final wall point
        if (N-2) not in center_points_idx:
            center_points_idx.append(N - 2)
        if (N-1) not in wall_points_idx:
            wall_points_idx.append(N - 1)

    # Get coordinates for these boundary points
    wall_x = x[wall_points_idx]
    wall_y = y[wall_points_idx]
    center_x = x[center_points_idx]
    center_y = y[center_points_idx] # Should all be ~0

    # Sort them by x-position to create a continuous line
    wall_x_sorted = wall_x[np.argsort(wall_x)]
    wall_y_sorted = wall_y[np.argsort(wall_x)]
    center_x_sorted = center_x[np.argsort(center_x)]
    center_y_sorted = center_y[np.argsort(center_x)] # all ~0
    
    # --- 2. Create the Contour Plot ---

    # Identify exit x-location
    x_exit = wall_x_sorted[-1]

    # Take final Mach value at exit (from last wall point)
    M_exit = M[wall_points_idx[-1]]

    # Create artificial exit line between centerline and wall
    y_exit = np.linspace(0, wall_y_sorted[-1], 30)
    x_exit_arr = np.full_like(y_exit, x_exit)
    M_exit_arr = np.full_like(y_exit, M_exit)

    # Append to existing data
    x = np.concatenate([x, x_exit_arr])
    y = np.concatenate([y, y_exit])
    M = np.concatenate([M, M_exit_arr])

    # Create a triangulation of all the (x, y) points
    triang = mtri.Triangulation(x, y)
    
    # Define contour levels for a smooth gradient
    levels = np.linspace(M.min(), M.max(), 40)
    
    # Plot the filled contours
    contour = ax.tricontourf(triang, M, levels=levels)

    # --- 3. Create a Clipping Path ---
    # This ensures the contours only appear *inside* the nozzle shape
    
    # Create a polygon from the wall and centerline vertices
    # 1. Start at (0,0)
    verts = [(0, 0)]
    # 2. Add all centerline points
    verts.extend(list(zip(center_x_sorted, center_y_sorted)))
    # 3. Add all wall points in reverse order
    verts.extend(list(zip(wall_x_sorted[::-1], wall_y_sorted[::-1])))
    # 4. Close the polygon back at (0,0)
    verts.append((0, 0))
    
    # Create the Path object and set it as the clip path
    clip_path = Path(verts)
    patch = plt.Polygon(verts, facecolor='none', edgecolor='none')
    ax.add_patch(patch)
    
    for collection in contour.collections:
        collection.set_clip_path(patch)

    # --- 4. Plot the Boundaries on Top ---
    ax.plot(wall_x_sorted, wall_y_sorted, 'k-', linewidth=2, label='Nozzle Wall')
    ax.plot(center_x_sorted, center_y_sorted, 'k-', linewidth=2, label='Centerline')

    # --- 5. Final Plot Settings ---
    ax.set_xlabel('Axial Position (x/r_throat)')
    ax.set_ylabel('Radial Position (y/r_throat)')
    ax.set_title('Mach Number Contours in MoC Nozzle')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.set_aspect('equal', adjustable='box')
    ax.legend(loc='best', fontsize='small')
    
    # Add a colorbar
    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('Mach Number (M)')
    
    # Save the plot
    plot_file_name = "MoC_Contour_Plot.png"
    plt.savefig(plot_file_name, dpi=300)
    plt.show()
    plt.close()

    print(f"Contour plot successfully saved as {plot_file_name}")
        
def get_matrix():
    x = int(input("How many characteristics lines?"))
    N = x
    lines = x
    point = 0
    while x > 0:
        point += (x+1)
        x -= 1
    N += point
    return np.zeros((N, 10)), N, lines
           
if __name__ == '__main__':
    main()