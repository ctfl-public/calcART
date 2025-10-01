import os
import subprocess

def check_mydir():
    calcART_data_DIR = os.getenv('calcART_data_DIR')
    if calcART_data_DIR:
        mydir = calcART_data_DIR
    else:
        print("Environment variable 'calcART_data_DIR' is not set. Use local dump directory, instead.")
        mydir = os.path.join(os.getcwd(), "dump")
    if not os.path.exists(mydir):
        os.makedirs(mydir)
    return mydir

# Get current git commits
def get_git_commit(repo_path):
    try:
        result = subprocess.run(['git', 'rev-parse', 'HEAD'], 
                                cwd=repo_path, capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        return None

def check_version_file(mydir):
	"""
	Check version.txt file in mydir for git commit numbers of RMCRT and calcART.
	If commits don't match current commits, warn user and provide choice to replace or append.
	"""

	import shutil
	from datetime import datetime
	
	version_file = os.path.join(mydir, "version.txt")
	
	
	# Get current commits
	calcART_dir = os.getenv('calcART_DIR')
	if calcART_dir:
		calcART_commit = get_git_commit(calcART_dir)
	else:
		raise EnvironmentError("Environment variable 'calcART_DIR' is not set.")
	RMCRT_DIR = os.getenv('RMCRT_DIR')
	if RMCRT_DIR:
		RMCRT_commit = get_git_commit(RMCRT_DIR)
	else:
		raise EnvironmentError("Environment variable 'RMCRT_DIR' is not set.")
	if not calcART_commit or not RMCRT_commit:
		print("Warning: Could not retrieve git commit information")
		return
	
	current_info = f"calcART: {calcART_commit}\nRMCRT: {RMCRT_commit}\nTimestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
	
	# Check if version.txt exists
	if not os.path.exists(version_file):
		print(f"Creating version.txt in {mydir}")
		os.makedirs(mydir, exist_ok=True)
		with open(version_file, 'w') as f:
			f.write(current_info)
		return
	
	# Read existing version file
	with open(version_file, 'r') as f:
		existing_content = f.read()
	
	# Check if current commits are already in the file
	if calcART_commit in existing_content and RMCRT_commit in existing_content:
		# Current commits are already recorded in version.txt
		return
	
	# Commits don't match - warn user and get choice
	print("WARNING: Current git commits differ from those in version.txt")
	print(f"\nCurrent commits:")
	print(f"calcART: {calcART_commit}")
	print(f"RMCRT:  {RMCRT_commit}")
	print(f"\nExisting version.txt content:")
	print(existing_content)
	
	while True:
		choice = input("\nChoose action:\n[r] Replace - Delete all files in mydir and create new version.txt\n[a] Append - Add current commits to existing version.txt\n[c] Cancel\nEnter choice (r/a/c): ").lower().strip()
		
		if choice == 'r':
			# Replace: delete all files and create new version.txt
			print(f"Deleting all files in {mydir}...")
			for filename in os.listdir(mydir):
				file_path = os.path.join(mydir, filename)
				try:
					if os.path.isfile(file_path) or os.path.islink(file_path):
						os.unlink(file_path)
					elif os.path.isdir(file_path):
						shutil.rmtree(file_path)
				except Exception as e:
					print(f"Failed to delete {file_path}: {e}")
			
			# Create new version.txt
			with open(version_file, 'w') as f:
				f.write(current_info)
			print("All files deleted. New version.txt created with current commits.")
			break
			
		elif choice == 'a':
			# Append: add current commits to existing file
			with open(version_file, 'a') as f:
				f.write(current_info)
			print("Current commits appended to version.txt")
			break
			
		elif choice == 'c':
			raise SystemExit("Failed to initialize calcART data management.")
			
		else:
			print("Invalid choice. Please enter 'r', 'a', or 'c'.")


def version_check(mydir):
	check_version_file(mydir)
