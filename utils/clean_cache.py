"""
Utility script to purge the temporary binary meshes cached by the Octave framework.
Frees up disk space after large batch generations.
"""
import shutil
from pathlib import Path

def main():
    root_dir = Path(__file__).resolve().parent.parent
    cache_dir = root_dir / "cached_meshes"

    print("\n--- Utility: Cache Cleaner ---")

    if not cache_dir.exists():
        print("[i] Cache directory does not exist. Nothing to clean.\n")
        return

    # Calculate size before deletion
    total_size_bytes = sum(f.stat().st_size for f in cache_dir.glob('**/*') if f.is_file())
    total_size_mb = total_size_bytes / (1024 * 1024)

    try:
        shutil.rmtree(cache_dir)
        print(f"[✓] Successfully purged {total_size_mb:.2f} MB of cached binary data.")
        print(f"[✓] Directory '{cache_dir.name}' removed.\n")
    except Exception as e:
        print(f"[x] Failed to delete cache directory: {e}\n")

if __name__ == "__main__":
    main()