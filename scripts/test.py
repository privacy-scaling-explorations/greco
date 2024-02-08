import argparse
import subprocess
import sys
from random import getrandbits
from bfv.utils import find_odd_pairwise_coprimes

# Set up command line argument parsing
parser = argparse.ArgumentParser(description="Run script multiple times with dynamic input")
parser.add_argument('runs', type=int, help='Number of times to run the script')
args = parser.parse_args()

# Number of times to run the script, taken from CLI
runs = args.runs

# Command and arguments for the Python script
script_to_run = "scripts/circuit_sk.py"
args_base = [
    "-n", "1024",
    "-t", "65537",
    "-output_input", "./src/data/sk_enc_input.json",
    "-output_constant", "./src/constants/sk_enc.rs",
]

# Command for running cargo test
cargo_test_command = ["cargo", "test", "--release", "--", "--nocapture", "test_sk_enc_valid"]

error_log_file = "scripts/error_log.txt"

# Open the error log file
with open(error_log_file, "w") as error_log:
    # Loop to run the script multiple times
    for i in range(runs):
        print(f"Running iteration: {i + 1}/{runs}")  # Print the current iteration
        
        # Generate qis randomly
        start = getrandbits(59)  # Start with a random 59-bit number
        qis = find_odd_pairwise_coprimes(start, 15)  # Find 15 pairwise coprimes starting from `start`
        
        # Convert qis list to string format expected by the script
        qis_str = f"[{', '.join(map(str, qis))}]"
        
        # Complete the args with dynamically generated qis
        args = args_base + ["-qis", qis_str]
        try:
            # Run the script using the same Python interpreter that's running this script
            # and include the required arguments
            subprocess.run([sys.executable, script_to_run] + args, check=True)
            
            # After running the script, run cargo test
            print(f"Running cargo test for iteration: {i + 1}")
            result = subprocess.run(cargo_test_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            # Capture and write stderr and stdout to the error log file
            error_log.write(f"Error occurred in run {i + 1}:\n")
            error_log.write(f"STDOUT: {e.stdout.decode()}\n")
            error_log.write(f"STDERR: {e.stderr.decode()}\n")

print(f"Script executed {runs} times with errors logged to {error_log_file}")
