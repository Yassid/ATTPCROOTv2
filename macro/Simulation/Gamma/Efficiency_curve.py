import os
import csv
import re
from pushbullet import Pushbullet
import numpy as np



array2= np.linspace(0.0001, 0.25, 30)
array3 = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9])
array4 = np.linspace(1, 5, 25)
#array5 = np.arange(3, 10, 0.5)
#energy_values = np.concatenate((array2,array3,array4))
energy_values = [0.06903793103448276]

no_events = [100000]

#pb = Pushbullet("o.iWhnN6tpdFI9kkPnu9961eqJjYMIniMC")

output_lines = []
output_lines2 = []

# Column headers
headers = ["no_event","Energy","PhotopeakCount","Efficiency","Error"]

# File path
file_path = "./EfficiencyCurves/efficiency_curve_PxCT_Peek.csv"


# Check if the file already exists
if not os.path.isfile(file_path):
    # Open the file in write mode and write the headers
    with open(file_path, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)


for energy in energy_values:
    for no_event in no_events:
        os.system(f"root -l -b -q 'gamma_sim.C({energy},{no_event})'")
        
        output = os.popen(f"root -l -b -q 'Simp_gamma_analysis.C({energy},{no_event})'").read()

        # Extract all relevant values from the terminal output
        num_ev_match = re.search(r"Total number of events : (\d+)", output)
        PhotopeakCount_match = re.search(r"Number of events in photopeak : (\d+)", output)
        efficiency_match = re.search(r"Photopeak Efficency : ([\d.]+)%", output)
        error_match = re.search(r"Error: ([\d.]+)", output)

        num_ev = num_ev_match.group(1) if num_ev_match else "N/A"
        PhotopeakCount = PhotopeakCount_match.group(1) if PhotopeakCount_match else "N/A"
        efficiency = efficiency_match.group(1) if efficiency_match else "N/A"
        error = error_match.group(1) if error_match else "N/A"

        # Print values
        print(f"Energy: {energy}")
        print(f"Total number of events: {num_ev}")
        print(f"Number of events in photopeak: {PhotopeakCount}")
        print(f"Photopeak Efficiency: {efficiency}%")
        print(f"Error: {error}\n")

        # Store output line for CSV
        output_lines.append([no_event, energy, PhotopeakCount, efficiency, error])

        # Open the CSV file in append mode and write the output line

        with open(file_path, "a", newline='') as f:
            writer = csv.writer(f)
            writer.writerow([no_event, energy, PhotopeakCount, efficiency, error])


#pb.push_note("Simulation Done", "The simulation and analysis are complete!")

os.system("more " + file_path)
