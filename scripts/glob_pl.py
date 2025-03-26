from pathlib import Path

def glob():
    # Find all .pl files in the current directory
    pl_files_i = Path('.').glob('*.pl')
    
    # Print the results in Meson-compatible format
    for i in pl_files_i:
        print(i)


if __name__ == "__main__":
    glob()