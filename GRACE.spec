# GRACE.spec
# ----------
# PyInstaller build spec for GRACE (PyQt6 version) — Multi-platform build.
#
# Produces platform-specific builds:
#   Windows: dist/GRACE/ folder with GRACE.exe
#   macOS:   dist/GRACE.app bundle
#   Linux:   dist/GRACE/ folder with GRACE executable
#
# Build with:
#   python build.py
#
# Or directly:
#   pyinstaller GRACE.spec
#
# Before building, ensure BLAST+ binaries are in the correct platform directory:
#   Windows: blast/windows/blastn.exe, makeblastdb.exe, etc.
#   macOS:   blast/macos/blastn, makeblastdb, etc.
#   Linux:   blast/linux/blastn, makeblastdb, etc.

import sys
import os
from PyInstaller.utils.hooks import collect_all, collect_submodules, collect_data_files

# ---------------------------------------------------------------------------
# Platform detection
# ---------------------------------------------------------------------------

IS_WINDOWS = sys.platform.startswith("win")
IS_MAC     = sys.platform == "darwin"
IS_LINUX   = sys.platform.startswith("linux")

# ---------------------------------------------------------------------------
# Collect packages PyInstaller cannot fully analyse automatically
# ---------------------------------------------------------------------------

qt_datas,       qt_binaries,       qt_hiddenimports       = collect_all("PyQt6")
primer3_datas,  primer3_binaries,  primer3_hiddenimports  = collect_all("primer3")
psutil_datas,   psutil_binaries,   psutil_hiddenimports   = collect_all("psutil")

pandas_hiddenimports = collect_submodules("pandas")
numpy_hiddenimports  = collect_submodules("numpy")

# ---------------------------------------------------------------------------
# Application source files (platform-specific BLAST binaries)
# ---------------------------------------------------------------------------

# Determine which BLAST directory to include
if IS_WINDOWS:
    blast_dir = ("blast/windows", "blast/windows")
elif IS_MAC:
    blast_dir = ("blast/macos", "blast/macos")
elif IS_LINUX:
    blast_dir = ("blast/linux", "blast/linux")
else:
    blast_dir = None
    print(f"WARNING: Unknown platform {sys.platform}, BLAST binaries will not be included!")

app_datas = [
    ("main.py",          "."),
    ("app_state.py",     "."),
    ("ui",               "ui"),
    ("core",             "core"),
]

# Only bundle assets/ if it actually exists — PyInstaller hard-fails if a
# datas source path doesn't exist, and this folder is optional (fonts/
# icons). Without this check, a fresh git clone/CI checkout without an
# assets/ folder committed would break the build on every platform.
if os.path.isdir("assets"):
    app_datas.append(("assets", "assets"))

# Add BLAST binaries if platform was detected
if blast_dir:
    app_datas.append(blast_dir)

# ---------------------------------------------------------------------------
# Hidden imports
# ---------------------------------------------------------------------------

all_hidden = (
    qt_hiddenimports
    + primer3_hiddenimports
    + psutil_hiddenimports
    + pandas_hiddenimports
    + numpy_hiddenimports
    + [
        "multiprocessing",
        "multiprocessing.pool",
        "multiprocessing.managers",
        "multiprocessing.reduction",
        "multiprocessing.resource_tracker",
        "multiprocessing.popen_fork",
        "multiprocessing.popen_forkserver",
        "multiprocessing.popen_spawn_posix",
        "subprocess",
        "tempfile",
        "socket",
        "threading",
        "dataclasses",
        "re",
        "io",
        "json",
        "copy",
        "shutil",
        "itertools",
        "functools",
        "bisect",
        "collections",
    ]
)

# Add Windows-specific imports only on Windows
if IS_WINDOWS:
    all_hidden.append("multiprocessing.popen_spawn_win32")

# ---------------------------------------------------------------------------
# Combined datas and binaries
# ---------------------------------------------------------------------------

all_datas = (
    app_datas
    + qt_datas
    + primer3_datas
    + psutil_datas
)

all_binaries = (
    qt_binaries
    + primer3_binaries
    + psutil_binaries
)

# ---------------------------------------------------------------------------
# Platform-specific icon
# ---------------------------------------------------------------------------

if IS_WINDOWS:
    icon_path = "assets/icon.ico" if os.path.exists("assets/icon.ico") else None
elif IS_MAC:
    icon_path = "assets/icon.icns" if os.path.exists("assets/icon.icns") else None
elif IS_LINUX:
    icon_path = "assets/icon.png" if os.path.exists("assets/icon.png") else None
else:
    icon_path = None

# ---------------------------------------------------------------------------
# ANALYSIS
# ---------------------------------------------------------------------------

a = Analysis(
    ["main.py"],
    pathex=["."],
    binaries=all_binaries,
    datas=all_datas,
    hiddenimports=all_hidden,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        "tkinter",
        "matplotlib",
        "IPython",
        "jupyter",
        "notebook",
        "scipy",
        "sklearn",
        "tensorflow",
        "torch",
        "streamlit",
        "webview",
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="GRACE",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    windowed=True,
    icon=icon_path,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="GRACE",
)

# macOS app bundle (only built on macOS)
if IS_MAC:
    app = BUNDLE(
        coll,
        name="GRACE.app",
        icon=icon_path,
        bundle_identifier="com.grace.app",
        info_plist={
            "NSPrincipalClass":           "NSApplication",
            "NSHighResolutionCapable":    True,
            "CFBundleShortVersionString": "1.0.0",
        },
    )
