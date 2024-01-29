import subprocess
import sys

# Number of times to run the script
num_runs = 1000

# Path to your target script
script_path = "circuit.py"

# File to log errors
error_log_file = "error_log.txt"

# Get the path to the Python interpreter
python_executable = sys.executable

# Open the error log file for writing
with open(error_log_file, "w") as log_file:
    for i in range(num_runs):
        print(f"Running iteration {i + 1}/{num_runs}...")
        result = subprocess.run([python_executable, script_path], capture_output=True, text=True)
        
        # Check if there was an error
        if result.returncode != 0:
            error_message = f"Error in iteration {i + 1}:\n{result.stderr}\n"
            log_file.write(error_message)
        else:
            print(f"Iteration {i + 1} completed successfully.")

print("All iterations completed.")
