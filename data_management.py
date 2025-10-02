import os
import subprocess
import shutil
from datetime import datetime

calcART_dump_DIR_default = "dump"


# Get current git commits
def get_git_commit(repo_path):
    try:
        result = subprocess.run(['git', 'rev-parse', 'HEAD'], 
                                cwd=repo_path, capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        return None

def has_uncommitted_changes(repo_path):
	"""Return True if the repository at repo_path has staged or unstaged changes."""
	try:
		# --porcelain gives a stable, machine-readable output; any non-empty output means changes
		result = subprocess.run(['git', 'status', '--porcelain'], cwd=repo_path,
								capture_output=True, text=True, check=True)
		return bool(result.stdout.strip())
	except subprocess.CalledProcessError:
		return False

def get_git_diff(repo_path):
	"""Return combined diff (staged + unstaged) for repo relative to HEAD. Empty string if none or error."""
	try:
		# git diff HEAD covers both staged and unstaged relative to HEAD
		result = subprocess.run(['git', 'diff', 'HEAD'], cwd=repo_path,
								capture_output=True, text=True, check=True)
		return result.stdout
	except subprocess.CalledProcessError:
		return ""

def check_version_file():
	"""
	Check and manage version tracking for calcART and RMCRT repositories.
	
	Workflow:
	1. Retrieve current git commits from calcART_DIR and RMCRT_DIR environment variables
	2. Check for uncommitted changes in both repositories
	3. Determine output directory:
	   - If uncommitted changes exist: use calcART_dump_DIR or fallback to local 'dump'
	   - If clean repositories: use calcART_data_DIR or fallback to local 'dump'
	4. Handle version.txt file:
	   - If doesn't exist: create new with current commits
	   - If exists and commits match: continue with existing
	   - If exists but commits differ: prompt user to replace or append
	5. Save git diffs if uncommitted changes are detected
	6. Return the determined output directory path
	
	Returns:
		str: Path to the data directory (mydir) where version.txt and outputs should be stored
	
	Raises:
		EnvironmentError: If required environment variables are not set
		SystemExit: If user chooses to quit or directory creation fails
	"""
	
	# Get current commits
	calcART_DIR = os.getenv('calcART_DIR')
	if calcART_DIR:
		calcART_commit = get_git_commit(calcART_DIR)
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

	timeStamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	current_info = f"\nTimestamp: {timeStamp}\ncalcART: {calcART_commit}\nRMCRT: {RMCRT_commit}\n"

	# Early: handle uncommitted changes scenario (before normal commit mismatch flow)
	calcART_dirty = has_uncommitted_changes(calcART_DIR)
	RMCRT_dirty = has_uncommitted_changes(RMCRT_DIR)

	# create mydir
	if calcART_dirty or RMCRT_dirty:
		# If either repo is dirty, change to dump directory
		print("WARNING: Detected uncommitted changes in the following repositories:")
		if calcART_dirty:
			print(" - calcART repository has uncommitted changes")
		if RMCRT_dirty:
			print(" - RMCRT repository has uncommitted changes")

		# Determine mydir for dirty state handling
		calcART_dump_dir = os.getenv('calcART_dump_DIR') # relative to cwd
		if calcART_dump_dir:
			print("Redirecting output to 'calcART_dump_DIR' directory.")
			mydir = os.path.join(os.getcwd(), calcART_dump_dir)
		else:
			print("Environment variable 'calcART_dump_DIR' is not set. Using local 'dump' directory instead.")
			mydir = os.path.join(os.getcwd(), calcART_dump_DIR_default)
	else:
		calcART_data_DIR = os.getenv('calcART_data_DIR')
		if calcART_data_DIR:
			mydir = calcART_data_DIR
		else:
			print("Environment variable 'calcART_data_DIR' is not set. Use local dump directory, instead.")
			mydir = os.path.join(os.getcwd(), calcART_dump_DIR_default)
		if not os.path.exists(mydir):
			os.makedirs(mydir)


	# Create mydir if it doesn't exist
	try:
		os.makedirs(mydir, exist_ok=True)
	except Exception as e:
		print(f"Could not create directory '{mydir}': {e}")
		raise SystemExit("Failed to initialize calcART data management.")

	# Create version.txt if it doesn't exist
	# if does not exist, create it with current commits and return.
	version_file = os.path.join(mydir, "version.txt")
	if not os.path.exists(version_file):
		print(f"Creating new version.txt in {mydir}")
		note = append_dirty_notes(current_info, calcART_dirty, RMCRT_dirty)
		update_version_file(version_file, note, 'w')
		if calcART_dirty or RMCRT_dirty:
			save_diff(mydir, calcART_DIR, RMCRT_DIR, calcART_dirty, RMCRT_dirty)
		return mydir
	
	# if exists, read version file
	with open(version_file, 'r') as f:
		existing_content = f.read()
	
	# Check if current commits are already in the file
	if calcART_commit in existing_content and RMCRT_commit in existing_content:
		# Current commits are already recorded in version.txt

		# If dirty, still save diff snapshot
		if calcART_dirty or RMCRT_dirty:
			save_diff(mydir, calcART_DIR, RMCRT_DIR, calcART_dirty, RMCRT_dirty)
		return mydir
	
	# Commits don't match - warn user and get choice
	print(f"WARNING: Current git commits differ from those in {version_file}")
	print(f"\nCurrent commits:")
	print(   "----------------")
	print(f"calcART: {calcART_commit}")
	print(f"RMCRT:  {RMCRT_commit}")
	print(f"\nExisting version.txt content:")
	print(   "-----------------------------")
	print(existing_content)
	
	# Prompt user for action
	while True:
		choice = input("\nChoose action:\n[r] Replace - Delete all files in mydir and create new version.txt\n[a] Append - Add current commits to existing version.txt\n[q] Quit\nEnter choice (r/a/q): ").lower().strip()
		
		if choice == 'r':
			# Replace: delete all files and create new version.txt
			print(f"Deleting all files and folders in {mydir}...")
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
			note = append_dirty_notes(current_info, calcART_dirty, RMCRT_dirty)
			update_version_file(version_file, note, 'w')

			print("All files deleted. New version.txt created with current commits.")
			break
			
		elif choice == 'a':
			# Append: add current commits to existing file
			note = append_dirty_notes(current_info, calcART_dirty, RMCRT_dirty)
			update_version_file(version_file, note, 'a')

			print("Current commits appended to version.txt")
			break
			
		elif choice == 'q':
			raise SystemExit("Failed to initialize calcART data management.")
			
		else:
			print("Invalid choice. Please enter 'r', 'a', or 'c'.")

	if calcART_dirty or RMCRT_dirty:
		save_diff(mydir, calcART_DIR, RMCRT_DIR, calcART_dirty, RMCRT_dirty)
	return mydir


def update_version_file(version_file, note, mode):
	with open(version_file, mode) as f:
		f.write(note)

def append_dirty_notes(note, calcART_dirty, RMCRT_dirty):
	if calcART_dirty:
		note += "calcART repository had uncommitted changes when recorded.\n"
	if RMCRT_dirty:
		note += "RMCRT repository had uncommitted changes when recorded.\n"
	return note

def save_diff(mydir, calcART_DIR, RMCRT_DIR, calcART_dirty, RMCRT_dirty):
	timeStamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	diff_dir = os.path.join(mydir, f'git_diff {timeStamp}')
	os.makedirs(diff_dir, exist_ok=True)
	if calcART_dirty:
		with open(os.path.join(diff_dir, 'calcART.diff'), 'w') as df:
			df.write(get_git_diff(calcART_DIR))
	if RMCRT_dirty:
		with open(os.path.join(diff_dir, 'RMCRT.diff'), 'w') as df:
			df.write(get_git_diff(RMCRT_DIR))
	print(f"Snapshot with diffs stored in {mydir}")