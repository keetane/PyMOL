import subprocess
from pymol import cmd

def sub(command: str):
    try:
        result = subprocess.run(command, capture_output=True, text=True, shell=True)
        if result.returncode != 0:
            print(f"Command failed with return code {result.returncode}")
            print(result.stderr)
        else:
            print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")