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

def check_version_file(mydir):
	"""
	Check version.txt file in mydir for git commit numbers of RMCRT and calcART.
	If commits don't match current commits, warn user and provide choice to replace or append.
	Additionally, if either repository has uncommitted (staged or unstaged) changes, warn user and prompt to:
	 - Append current (dirty) state info to existing version.txt in same directory, OR
	 - Use a custom folder to snapshot the current commits and store diffs of uncommitted changes.
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
	
	# Early: handle uncommitted changes scenario (before normal mismatch flow)
	calcART_dirty = has_uncommitted_changes(calcART_dir)
	RMCRT_dirty = has_uncommitted_changes(RMCRT_DIR)

	if calcART_dirty or RMCRT_dirty:
		print("WARNING: Detected uncommitted changes in the following repositories:")
		if calcART_dirty:
			print(" - calcART repository has uncommitted changes")
		if RMCRT_dirty:
			print(" - RMCRT repository has uncommitted changes")

		# Ensure base dir exists for potential append action
		os.makedirs(mydir, exist_ok=True)

		# Prompt user for action
		while True:
			choice_dirty = input("\nChoose action for dirty working tree(s):\n[a] Append dirty state info to existing version.txt in current directory\n[c] Use custom folder (snapshot)\n[x] Cancel\nEnter choice (a/c/x): ").lower().strip()
			if choice_dirty == 'a':
				# Append (or create) version.txt in current dir with dirty note
				note = current_info
				if calcART_dirty:
					note += "calcART repository had uncommitted changes when recorded.\n"
				if RMCRT_dirty:
					note += "RMCRT repository had uncommitted changes when recorded.\n"
				note += "(Diffs not archived; choose custom folder option next time to snapshot diffs)\n\n"
				mode = 'a' if os.path.exists(version_file) else 'w'
				with open(version_file, mode) as f:
					f.write(note)
				print("Dirty state info appended to version.txt")
				return
			elif choice_dirty == 'c':
				# Custom folder snapshot
				while True:
					custom_path = input("Enter custom folder path relative to cwd (or blank to cancel): ").strip()
					if not custom_path:
						print("Cancelled custom folder entry. Aborting initialization.")
						raise SystemExit("Failed to initialize calcART data management.")
					try:
						custom_path = os.path.join(os.getcwd(), custom_path)
						os.makedirs(custom_path, exist_ok=True)
						break
					except Exception as e:
						print(f"Could not create directory '{custom_path}': {e}")

				# Write version.txt in custom folder
				custom_version_file = os.path.join(custom_path, 'version.txt')
				with open(custom_version_file, 'w') as f:
					f.write(current_info)
					if calcART_dirty:
						f.write("calcART repository had uncommitted changes. See git_diff/calcART.diff.\n")
					if RMCRT_dirty:
						f.write("RMCRT repository had uncommitted changes. See git_diff/RMCRT.diff.\n")
					f.write("\n")

				# Collect diffs
				diff_dir = os.path.join(custom_path, 'git_diff')
				os.makedirs(diff_dir, exist_ok=True)
				if calcART_dirty:
					with open(os.path.join(diff_dir, 'calcART.diff'), 'w') as df:
						df.write(get_git_diff(calcART_dir))
				if RMCRT_dirty:
					with open(os.path.join(diff_dir, 'RMCRT.diff'), 'w') as df:
						df.write(get_git_diff(RMCRT_DIR))
				print(f"Snapshot with diffs stored in {custom_path}")
				return custom_path
			elif choice_dirty == 'x':
				raise SystemExit("Failed to initialize calcART data management.")
			else:
				print("Invalid choice. Please enter 'a', 'c', or 'x'.")

	# Check if version.txt exists (normal path when no uncommitted changes and proceeding)
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
	
	# Prompt user for action
	while True:
		choice = input("\nChoose action:\n[r] Replace - Delete all files in mydir and create new version.txt\n[a] Append - Add current commits to existing version.txt\n[x] Exit\nEnter choice (r/a/x): ").lower().strip()
		
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
			
		elif choice == 'x':
			raise SystemExit("Failed to initialize calcART data management.")
			
		else:
			print("Invalid choice. Please enter 'r', 'a', or 'c'.")

