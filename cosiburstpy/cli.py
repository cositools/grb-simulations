import sys
import subprocess

def main():
    if len(sys.argv) < 2:
        print("Usage: cosiburstpy [setup]")
        sys.exit(1)

    command = sys.argv[1]

    if command == "setup":
        run_setup()
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)


def run_setup():
    try:
        subprocess.run(["gdt-data", "init"], check=True)
    except FileNotFoundError:
        print("Error: gdt-data is not installed.")
        print("Install with: pip install cosiburstpy[fermi-gbm]")
        sys.exit(1)