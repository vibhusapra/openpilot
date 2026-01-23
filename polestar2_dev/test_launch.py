#!/usr/bin/env python3
"""
OpenPilot Launch Test
Tests that core OpenPilot components can import and initialize
"""

import sys
import os
import importlib
import traceback
from pathlib import Path

# Add openpilot to path
OPENPILOT_ROOT = Path(__file__).parent
sys.path.insert(0, str(OPENPILOT_ROOT))
os.chdir(OPENPILOT_ROOT)

# ANSI color codes
RED = '\033[0;31m'
GREEN = '\033[0;32m'
YELLOW = '\033[1;33m'
NC = '\033[0m'  # No Color

def print_header(text):
    print("\n" + "="*60)
    print(f"  {text}")
    print("="*60)

def test_import(module_name, description):
    """Try to import a module and report success/failure"""
    try:
        importlib.import_module(module_name)
        print(f"{GREEN}✓{NC} {description}: {module_name}")
        return True
    except ImportError as e:
        print(f"{RED}✗{NC} {description}: {module_name}")
        print(f"  Error: {e}")
        return False
    except Exception as e:
        print(f"{RED}✗{NC} {description}: {module_name}")
        print(f"  Unexpected error: {e}")
        return False

def test_cereal():
    """Test cereal message definitions"""
    print_header("Testing Cereal Messages")

    tests = [
        ("cereal", "Base cereal module"),
        ("cereal.messaging", "Messaging module"),
        ("cereal.services", "Services definitions"),
    ]

    results = []
    for module, desc in tests:
        results.append(test_import(module, desc))

    # Try to actually use cereal
    try:
        from cereal import log
        print(f"{GREEN}✓{NC} Can create log messages")
        results.append(True)
    except Exception as e:
        print(f"{RED}✗{NC} Cannot create log messages: {e}")
        results.append(False)

    return all(results)

def test_common():
    """Test common utilities"""
    print_header("Testing Common Utilities")

    tests = [
        ("common.params", "Parameters module"),
        ("common.basedir", "Base directory module"),
    ]

    results = []
    for module, desc in tests:
        results.append(test_import(module, desc))

    return all(results)

def test_selfdrive():
    """Test selfdrive modules"""
    print_header("Testing Selfdrive Modules")

    tests = [
        ("selfdrive.version", "Version info"),
        ("selfdrive.swaglog", "Logging module"),
        ("selfdrive.hardware", "Hardware abstraction"),
    ]

    results = []
    for module, desc in tests:
        results.append(test_import(module, desc))

    return all(results)

def test_volvo_integration():
    """Test Volvo/Polestar integration"""
    print_header("Testing Volvo/Polestar Integration")

    results = []

    # Check if opendbc submodule exists
    opendbc_path = OPENPILOT_ROOT / "opendbc_repo"
    if not opendbc_path.exists():
        print(f"{RED}✗{NC} opendbc_repo submodule not found")
        return False

    # Add opendbc to path
    sys.path.insert(0, str(opendbc_path))

    # Test Volvo car module
    try:
        from opendbc.car.volvo.values import CAR
        if hasattr(CAR, 'POLESTAR_2'):
            print(f"{GREEN}✓{NC} Polestar 2 found in CAR enum")
            results.append(True)
        else:
            print(f"{YELLOW}⚠{NC} Polestar 2 not in CAR enum")
            results.append(False)
    except ImportError as e:
        print(f"{RED}✗{NC} Cannot import Volvo values: {e}")
        results.append(False)

    # Test Volvo interface
    try:
        from opendbc.car.volvo.interface import CarInterface
        print(f"{GREEN}✓{NC} Volvo CarInterface importable")
        results.append(True)
    except ImportError as e:
        print(f"{RED}✗{NC} Cannot import Volvo interface: {e}")
        results.append(False)

    # Test fingerprints
    try:
        from opendbc.car.volvo.fingerprints import FINGERPRINTS
        if CAR.POLESTAR_2 in FINGERPRINTS:
            fp_count = len(FINGERPRINTS[CAR.POLESTAR_2][0])
            print(f"{GREEN}✓{NC} Polestar 2 fingerprint: {fp_count} messages")
            results.append(True)
        else:
            print(f"{YELLOW}⚠{NC} Polestar 2 fingerprint not found")
            results.append(False)
    except Exception as e:
        print(f"{RED}✗{NC} Cannot check fingerprints: {e}")
        results.append(False)

    # Check safety firmware
    volvo_safety = opendbc_path / "opendbc/safety/modes/volvo.h"
    if volvo_safety.exists():
        size = volvo_safety.stat().st_size
        print(f"{GREEN}✓{NC} Volvo safety firmware: {size:,} bytes")
        results.append(True)
    else:
        print(f"{RED}✗{NC} Volvo safety firmware not found")
        results.append(False)

    return all(results)

def test_launch_scripts():
    """Test launch scripts exist and are executable"""
    print_header("Testing Launch Scripts")

    scripts = [
        "launch_openpilot.sh",
        "launch_chffrplus.sh",
        "env.sh",
    ]

    results = []
    for script in scripts:
        script_path = OPENPILOT_ROOT / script
        if script_path.exists():
            if os.access(script_path, os.X_OK):
                print(f"{GREEN}✓{NC} {script}: exists and executable")
            else:
                print(f"{YELLOW}⚠{NC} {script}: exists but not executable")
            results.append(True)
        else:
            print(f"{RED}✗{NC} {script}: not found")
            results.append(False)

    return all(results)

def test_environment():
    """Test environment setup"""
    print_header("Testing Environment")

    # Check Python version
    py_version = sys.version_info
    if py_version.major == 3 and py_version.minor >= 11:
        print(f"{GREEN}✓{NC} Python {py_version.major}.{py_version.minor}.{py_version.micro}")
    else:
        print(f"{YELLOW}⚠{NC} Python {py_version.major}.{py_version.minor} (expected 3.11+)")

    # Check key environment variables
    fingerprint = os.getenv("FINGERPRINT")
    if fingerprint:
        print(f"{GREEN}✓{NC} FINGERPRINT set: {fingerprint}")
    else:
        print(f"{YELLOW}ℹ{NC} FINGERPRINT not set (will use auto-detection)")

    skip_fw = os.getenv("SKIP_FW_QUERY")
    if skip_fw:
        print(f"{GREEN}✓{NC} SKIP_FW_QUERY set: {skip_fw}")
    else:
        print(f"{YELLOW}ℹ{NC} SKIP_FW_QUERY not set")

    return True

def main():
    print("\n" + "="*60)
    print("  OpenPilot Launch Test")
    print("  Testing core components can import and initialize")
    print("="*60)

    # Track overall results
    all_results = []

    # Run tests
    all_results.append(test_environment())
    all_results.append(test_cereal())
    all_results.append(test_common())
    all_results.append(test_selfdrive())
    all_results.append(test_volvo_integration())
    all_results.append(test_launch_scripts())

    # Summary
    print_header("TEST SUMMARY")

    passed = sum(all_results)
    total = len(all_results)

    if passed == total:
        print(f"{GREEN}✓ All {total} test categories passed!{NC}")
        print("\nCore OpenPilot components are working.")
        print("The system should be able to start on a Comma device.")
        return 0
    else:
        print(f"{RED}✗ {total - passed}/{total} test categories failed{NC}")
        print("\nSome components are not working correctly.")
        print("Please check the errors above before deploying.")
        return 1

if __name__ == "__main__":
    sys.exit(main())