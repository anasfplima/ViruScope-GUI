import tempfile
import string
import random
import os
import shutil
import time
import threading

# Function to create a safe temporary directory as BLAST does not accept paths with spaces, underscores, etc
def safe_tempdir(base_dir):
    dirname = ''.join(random.choices(string.ascii_letters + string.digits, k=8))
    path = os.path.join(base_dir, dirname)
    os.makedirs(path, exist_ok=True)
    return path

def delayed_delete(path, delay=2):
    """Delete a folder after a short delay to allow file handles to close."""
    def _delete():
        time.sleep(delay)
        try:
            shutil.rmtree(path, ignore_errors=True)
            print(f"Deleted temp directory: {path}")
        except Exception as e:
            print(f"Could not delete {path}: {e}")
    threading.Thread(target=_delete, daemon=True).start()

def move_and_delete_later(path):
    """Move the directory to a temp location, then delete in background."""
    if not os.path.exists(path):
        return
    tmp_delete_dir = tempfile.mkdtemp(prefix="to_delete_")
    new_path = os.path.join(tmp_delete_dir, os.path.basename(path))
    try:
        shutil.move(path, new_path)
        delayed_delete(tmp_delete_dir)
    except Exception as e:
        print(f"Could not move {path}: {e}")

def safe_delete_folder(path, retries=5, delay=0.5):
    for i in range(retries):
        try:
            shutil.rmtree(path)
            print(f"[TempDir] Deleted: {path}")
            return
        except PermissionError:
            print(f"[TempDir] PermissionError deleting {path}, retry {i+1}/{retries}")
            time.sleep(delay)
    print(f"[TempDir] Could not delete {path} after {retries} retries.")