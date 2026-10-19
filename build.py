"""
build.py
--------
Build helper for GRACE (PyQt6 version).

Usage:
    python build.py           # full clean build
    python build.py --no-upx  # skip UPX compression
    python build.py --no-clean  # skip cleaning previous build

Output:
    dist/GRACE/GRACE.exe      (Windows)
    dist/GRACE/GRACE           (Mac/Linux)
    dist/GRACE.app/            (Mac only)
"""

import os
import sys
import shutil
import subprocess
import argparse


def run(cmd):
    print(f"\n>>> {' '.join(str(c) for c in cmd)}\n")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(result.returncode)


def clean():
    for path in ("build", "dist", "__pycache__"):
        if os.path.exists(path):
            print(f"Removing {path}/")
            shutil.rmtree(path)
    for root, dirs, files in os.walk("."):
        for f in files:
            if f.endswith(".pyc"):
                os.remove(os.path.join(root, f))
        dirs[:] = [d for d in dirs if d not in ("dist", "build", ".git", "venv", ".venv")]


def check_pyinstaller():
    try:
        import PyInstaller
        print(f"PyInstaller {PyInstaller.__version__} found.")
    except ImportError:
        print("ERROR: PyInstaller not found. Run: pip install pyinstaller")
        sys.exit(1)


def check_upx():
    return shutil.which("upx") is not None


def build(use_upx=True):
    if not use_upx:
        _patch_upx(False)
    try:
        run([sys.executable, "-m", "PyInstaller", "GRACE.spec", "--noconfirm"])
    finally:
        if not use_upx:
            _patch_upx(True)


def _patch_upx(enabled):
    with open("GRACE.spec", "r") as f:
        c = f.read()
    if enabled:
        c = c.replace("upx=False", "upx=True")
    else:
        c = c.replace("upx=True", "upx=False")
    with open("GRACE.spec", "w") as f:
        f.write(c)


def report():
    dist_dir = os.path.join("dist", "GRACE")
    if not os.path.exists(dist_dir):
        print("Build output not found.")
        return
    total_bytes = sum(
        os.path.getsize(os.path.join(dp, f))
        for dp, _, files in os.walk(dist_dir)
        for f in files
    )
    exe_name = "GRACE.exe" if sys.platform == "win32" else "GRACE"
    print("\n" + "=" * 60)
    print("BUILD COMPLETE")
    print("=" * 60)
    print(f"  Output  : dist/GRACE/")
    print(f"  Size    : {total_bytes / (1024*1024):.1f} MB")
    print(f"  Exe     : dist/GRACE/{exe_name}")
    if sys.platform == "darwin":
        print(f"  Bundle  : dist/GRACE.app/")
    print()
    print("  Distribute the entire dist/GRACE/ folder.")
    print("  BLAST+ binaries from blast/<platform>/ are bundled —")
    print("  users do not need to install BLAST+ separately.")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(description="Build GRACE with PyInstaller")
    parser.add_argument("--no-upx",   action="store_true")
    parser.add_argument("--no-clean", action="store_true")
    args = parser.parse_args()

    print("GRACE — PyInstaller Build")
    print(f"Platform : {sys.platform}")
    print(f"Python   : {sys.version.split()[0]}")

    check_pyinstaller()

    upx_available = check_upx()
    use_upx = upx_available and not args.no_upx
    if not upx_available:
        print("UPX not found — building without compression.")

    if not args.no_clean:
        print("\nCleaning previous build...")
        clean()

    print("\nRunning PyInstaller...")
    build(use_upx)
    report()


if __name__ == "__main__":
    main()