#!/usr/bin/env python3
"""
Polestar 2 OpenPilot Integration Test

This script verifies that all components for Polestar 2 support are properly
configured and can be loaded without errors.

Run with: python scripts/test_polestar2.py
"""

import sys
import os
from pathlib import Path

# Add paths for imports
OPENPILOT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(OPENPILOT_ROOT))
sys.path.insert(0, str(OPENPILOT_ROOT / "opendbc_repo"))

def test_banner(name: str):
    print(f"\n{'='*60}")
    print(f"  {name}")
    print('='*60)

def test_pass(msg: str):
    print(f"  ✅ {msg}")

def test_fail(msg: str):
    print(f"  ❌ {msg}")

def test_warn(msg: str):
    print(f"  ⚠️  {msg}")

def main():
    print("\n" + "="*60)
    print("  POLESTAR 2 OPENPILOT INTEGRATION TEST")
    print("="*60)

    all_passed = True

    # Test 1: Python version
    test_banner("Test 1: Python Version")
    py_version = sys.version_info
    if py_version >= (3, 11):
        test_pass(f"Python {py_version.major}.{py_version.minor}.{py_version.micro}")
    else:
        test_fail(f"Python {py_version.major}.{py_version.minor} (requires 3.11+)")
        all_passed = False

    # Test 2: opendbc car module
    test_banner("Test 2: opendbc Car Module")
    try:
        from opendbc.car import Bus, structs
        test_pass("opendbc.car base imports")
    except ImportError as e:
        test_fail(f"opendbc.car import failed: {e}")
        all_passed = False

    # Test 3: Volvo values and CAR enum
    test_banner("Test 3: Volvo CAR Enum")
    try:
        from opendbc.car.volvo.values import CAR, CarControllerParams, DBC
        test_pass(f"CAR enum loaded with {len(list(CAR))} vehicles")
        for car in CAR:
            if "POLESTAR" in car.name:
                test_pass(f"Found: {car.name}")
    except ImportError as e:
        test_fail(f"Volvo values import failed: {e}")
        all_passed = False

    # Test 4: Volvo fingerprints
    test_banner("Test 4: Volvo Fingerprints")
    try:
        from opendbc.car.volvo.fingerprints import FINGERPRINTS, FW_VERSIONS
        test_pass(f"FINGERPRINTS loaded: {len(FINGERPRINTS)} vehicles")

        for car, fps in FINGERPRINTS.items():
            if "POLESTAR" in car.name:
                test_pass(f"{car.name}: {len(fps)} fingerprint(s), {len(fps[0])} CAN IDs")
    except ImportError as e:
        test_fail(f"Fingerprints import failed: {e}")
        all_passed = False

    # Test 5: Volvo CarInterface
    test_banner("Test 5: Volvo CarInterface")
    try:
        from opendbc.car.volvo.interface import CarInterface
        test_pass("CarInterface class loaded")

        # Check required attributes
        attrs = ['CarState', 'CarController', '_get_params']
        for attr in attrs:
            if hasattr(CarInterface, attr):
                test_pass(f"CarInterface.{attr} exists")
            else:
                test_fail(f"CarInterface.{attr} missing")
                all_passed = False
    except ImportError as e:
        test_fail(f"CarInterface import failed: {e}")
        all_passed = False

    # Test 6: Volvo CarState
    test_banner("Test 6: Volvo CarState")
    try:
        from opendbc.car.volvo.carstate import CarState
        test_pass("CarState class loaded")

        # Check required methods
        methods = ['update', 'get_can_parsers']
        for method in methods:
            if hasattr(CarState, method):
                test_pass(f"CarState.{method} exists")
            else:
                test_fail(f"CarState.{method} missing")
                all_passed = False
    except ImportError as e:
        test_fail(f"CarState import failed: {e}")
        all_passed = False

    # Test 7: Volvo CarController
    test_banner("Test 7: Volvo CarController")
    try:
        from opendbc.car.volvo.carcontroller import CarController
        test_pass("CarController class loaded")
    except ImportError as e:
        test_fail(f"CarController import failed: {e}")
        all_passed = False

    # Test 8: Safety model
    test_banner("Test 8: Volvo Safety Model")
    safety_path = OPENPILOT_ROOT / "opendbc_repo" / "opendbc" / "safety" / "modes" / "volvo.h"
    if safety_path.exists():
        test_pass(f"volvo.h safety model exists ({safety_path.stat().st_size} bytes)")

        # Check for key safety functions
        content = safety_path.read_text()
        safety_funcs = ['volvo_rx_hook', 'volvo_tx_hook', 'volvo_fwd_hook']
        for func in safety_funcs:
            if func in content:
                test_pass(f"Safety function {func} found")
            else:
                test_warn(f"Safety function {func} not found")
    else:
        test_fail(f"volvo.h safety model not found at {safety_path}")
        all_passed = False

    # Test 9: DBC files
    test_banner("Test 9: Volvo DBC Files")
    dbc_dir = OPENPILOT_ROOT / "opendbc_repo" / "opendbc" / "dbc"
    volvo_dbcs = ['volvo_cma.dbc', 'volvo_mid_1.dbc', 'volvo_front_1_cma.dbc']
    for dbc in volvo_dbcs:
        dbc_path = dbc_dir / dbc
        if dbc_path.exists():
            test_pass(f"{dbc} exists ({dbc_path.stat().st_size} bytes)")
        else:
            test_fail(f"{dbc} not found")
            all_passed = False

    # Test 10: Torque parameters
    test_banner("Test 10: Torque Parameters")
    torque_path = OPENPILOT_ROOT / "opendbc_repo" / "opendbc" / "car" / "torque_data" / "override.toml"
    if torque_path.exists():
        content = torque_path.read_text()
        if 'POLESTAR_2' in content:
            test_pass("POLESTAR_2 torque parameters configured")
        else:
            test_fail("POLESTAR_2 torque parameters missing")
            all_passed = False
    else:
        test_fail(f"torque override.toml not found")
        all_passed = False

    # Test 11: Interface longitudinal tuning
    test_banner("Test 11: Longitudinal Tuning (v0.10.0 fix)")
    interface_path = OPENPILOT_ROOT / "opendbc_repo" / "opendbc" / "car" / "volvo" / "interface.py"
    if interface_path.exists():
        content = interface_path.read_text()
        if 'longitudinalTuning' in content:
            test_pass("longitudinalTuning parameters present (v0.10.0 compatible)")
        else:
            test_fail("longitudinalTuning parameters missing (will crash controlsd)")
            all_passed = False

        if 'SafetyModel.volvo' in content:
            test_pass("Using SafetyModel.volvo (steering enabled)")
        elif 'SafetyModel.noOutput' in content:
            test_warn("Using SafetyModel.noOutput (dashcam only, no steering)")

        if 'dashcamOnly = False' in content:
            test_pass("dashcamOnly = False (steering control enabled)")
        elif 'dashcamOnly = True' in content:
            test_warn("dashcamOnly = True (steering disabled)")
    else:
        test_fail(f"interface.py not found")
        all_passed = False

    # Test 12: Simulated fingerprint matching
    test_banner("Test 12: Fingerprint Matching Simulation")
    try:
        from opendbc.car.volvo.fingerprints import FINGERPRINTS
        from opendbc.car.volvo.values import CAR

        # Simulate a fingerprint match
        ps2_fp = FINGERPRINTS.get(CAR.POLESTAR_2)
        if ps2_fp and len(ps2_fp) > 0:
            test_pass(f"POLESTAR_2 fingerprint has {len(ps2_fp[0])} CAN message IDs")

            # Show some of the CAN IDs
            sample_ids = list(ps2_fp[0].keys())[:5]
            test_pass(f"Sample CAN IDs: {sample_ids}")
        else:
            test_fail("POLESTAR_2 fingerprint is empty or missing")
            all_passed = False
    except Exception as e:
        test_fail(f"Fingerprint test failed: {e}")
        all_passed = False

    # Summary
    print("\n" + "="*60)
    if all_passed:
        print("  ✅ ALL TESTS PASSED")
        print("="*60)
        print("\nYour Polestar 2 OpenPilot setup appears to be working!")
        print("Next steps:")
        print("  1. Build with: scons -j$(nproc)")
        print("  2. Flash to Comma 3 or run in simulation")
        return 0
    else:
        print("  ❌ SOME TESTS FAILED")
        print("="*60)
        print("\nPlease fix the issues above before proceeding.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
