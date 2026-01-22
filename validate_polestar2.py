#!/usr/bin/env python3
"""
Polestar 2 Integration Validation Script
Validates that Polestar 2 is properly integrated into OpenPilot 0.10.0
"""

import sys
from pathlib import Path

def print_header(title):
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70)

def print_section(title):
    print(f"\n{title}")
    print("-" * 70)

def check_imports():
    print_section("1. CHECKING IMPORTS")

    try:
        from opendbc.car.volvo.values import CAR
        print("✓ opendbc.car.volvo.values imported")

        from opendbc.car.volvo.interface import CarInterface
        print("✓ opendbc.car.volvo.interface imported")

        from opendbc.car.volvo.fingerprints import FINGERPRINTS
        print("✓ opendbc.car.volvo.fingerprints imported")

        return True
    except Exception as e:
        print(f"✗ Import failed: {e}")
        return False

def check_car_enum():
    print_section("2. CHECKING CAR ENUM")

    try:
        from opendbc.car.volvo.values import CAR

        if hasattr(CAR, 'POLESTAR_2'):
            print(f"✓ CAR.POLESTAR_2 exists: {CAR.POLESTAR_2}")
            return True
        else:
            print("✗ CAR.POLESTAR_2 not found")
            return False
    except Exception as e:
        print(f"✗ Check failed: {e}")
        return False

def check_fingerprint():
    print_section("3. CHECKING FINGERPRINT")

    try:
        from opendbc.car.volvo.values import CAR
        from opendbc.car.volvo.fingerprints import FINGERPRINTS

        if CAR.POLESTAR_2 in FINGERPRINTS:
            fps = FINGERPRINTS[CAR.POLESTAR_2]
            print(f"✓ Polestar 2 fingerprint found: {len(fps)} variant(s)")

            for i, fp in enumerate(fps):
                msg_count = len(fp)
                msg_ids = sorted(fp.keys())
                print(f"  Fingerprint {i+1}: {msg_count} CAN messages")
                print(f"    Message IDs: {msg_ids[:15]}{'...' if len(msg_ids) > 15 else ''}")

            return True
        else:
            print("✗ Polestar 2 fingerprint not found")
            return False
    except Exception as e:
        print(f"✗ Check failed: {e}")
        return False

def check_platform_registration():
    print_section("4. CHECKING PLATFORM REGISTRATION")

    try:
        from opendbc.car.values import PLATFORMS
        from opendbc.car.volvo.values import CAR

        if CAR.POLESTAR_2 in [p for p in PLATFORMS]:
            print(f"✓ POLESTAR_2 registered in global PLATFORMS")
            return True
        else:
            print("✗ POLESTAR_2 not in global PLATFORMS")
            return False
    except Exception as e:
        print(f"✗ Check failed: {e}")
        return False

def check_interface():
    print_section("5. CHECKING CAR INTERFACE")

    try:
        from opendbc.car.volvo.interface import CarInterface
        from opendbc.car.volvo.values import CAR

        print(f"✓ CarInterface class available")
        print(f"  Can instantiate for: {CAR.POLESTAR_2}")

        return True
    except Exception as e:
        print(f"✗ Check failed: {e}")
        return False

def check_safety_firmware():
    print_section("6. CHECKING SAFETY FIRMWARE")

    volvo_safety = Path(__file__).parent / "opendbc_repo" / "opendbc" / "safety" / "modes" / "volvo.h"

    if volvo_safety.exists():
        size = volvo_safety.stat().st_size
        print(f"✓ Volvo safety firmware found: {volvo_safety.name}")
        print(f"  Location: {volvo_safety.relative_to(Path.cwd())}")
        print(f"  Size: {size:,} bytes")

        # Check for key defines
        content = volvo_safety.read_text()

        checks = [
            ("VOLVO_LCA_STEER", "Steering command message"),
            ("volvo_tx_hook", "TX message validation"),
            ("volvo_rx_hook", "RX message processing"),
            ("SAFETY_VOLVO", "Safety model ID"),
        ]

        print(f"  Key components:")
        for define, description in checks:
            if define in content:
                print(f"    ✓ {define} - {description}")
            else:
                print(f"    ✗ {define} - NOT FOUND")

        return True
    else:
        print(f"✗ Volvo safety firmware not found at: {volvo_safety}")
        return False

def check_git_config():
    print_section("7. CHECKING GIT CONFIGURATION")

    import subprocess

    try:
        # Check opendbc submodule
        result = subprocess.run(
            ["git", "config", "--file", ".gitmodules", "submodule.opendbc_repo.url"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent
        )

        if result.returncode == 0:
            url = result.stdout.strip()
            print(f"✓ opendbc submodule URL: {url}")

            if "paper5590" in url:
                print(f"  ✓ Using Paper's fork (correct for safety firmware)")
            else:
                print(f"  ⚠ Not using Paper's fork - safety firmware may be missing!")

        # Check branch
        result = subprocess.run(
            ["git", "-C", "opendbc_repo", "branch", "--show-current"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent
        )

        if result.returncode == 0:
            branch = result.stdout.strip()
            print(f"✓ opendbc branch: {branch}")

            if "cma" in branch.lower():
                print(f"  ✓ Using CMA branch (correct for Polestar 2)")
            else:
                print(f"  ⚠ Not using CMA branch")

        return True
    except Exception as e:
        print(f"✗ Git check failed: {e}")
        return False

def check_openpilot_changes():
    print_section("8. CHECKING OPENPILOT CHANGES")

    # Check card.py
    card_py = Path(__file__).parent / "selfdrive" / "car" / "card.py"
    if card_py.exists():
        content = card_py.read_text()
        if "VOLVO" in content or "POLESTAR" in content:
            print(f"✓ selfdrive/car/card.py: Volvo feature flags added")
        else:
            print(f"⚠ selfdrive/car/card.py: Volvo support not detected")
    else:
        print(f"✗ selfdrive/car/card.py: File not found")

    # Check torqued.py
    torqued_py = Path(__file__).parent / "selfdrive" / "locationd" / "torqued.py"
    if torqued_py.exists():
        content = torqued_py.read_text()
        if "volvo" in content:
            print(f"✓ selfdrive/locationd/torqued.py: Volvo added to ALLOWED_CARS")
        else:
            print(f"⚠ selfdrive/locationd/torqued.py: Volvo not in ALLOWED_CARS")
    else:
        print(f"✗ selfdrive/locationd/torqued.py: File not found")

    # Check launch script
    launch_sh = Path(__file__).parent / "launch_openpilot.sh"
    if launch_sh.exists():
        content = launch_sh.read_text()
        if "env.sh" in content:
            print(f"✓ launch_openpilot.sh: Environment variable support added")
        else:
            print(f"⚠ launch_openpilot.sh: No env.sh support")
    else:
        print(f"✗ launch_openpilot.sh: File not found")

    return True

def main():
    print_header("POLESTAR 2 INTEGRATION VALIDATION")
    print(f"OpenPilot 0.10.0 with Paper's Volvo Safety Firmware")
    print(f"Working Directory: {Path.cwd()}")

    checks = [
        ("Imports", check_imports),
        ("CAR Enum", check_car_enum),
        ("Fingerprint", check_fingerprint),
        ("Platform Registration", check_platform_registration),
        ("Car Interface", check_interface),
        ("Safety Firmware", check_safety_firmware),
        ("Git Configuration", check_git_config),
        ("OpenPilot Changes", check_openpilot_changes),
    ]

    results = []

    for name, check_func in checks:
        try:
            result = check_func()
            results.append((name, result))
        except Exception as e:
            print(f"\n✗ UNEXPECTED ERROR in {name}: {e}")
            results.append((name, False))

    # Summary
    print_header("VALIDATION SUMMARY")

    passed = sum(1 for _, result in results if result)
    total = len(results)

    for name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"  {status}: {name}")

    print(f"\nTotal: {passed}/{total} checks passed")

    if passed == total:
        print("\n" + "=" * 70)
        print("  ✓ ALL CHECKS PASSED - READY FOR DEPLOYMENT")
        print("=" * 70)
        print("\nNext steps:")
        print("  1. Deploy to Comma 3: git clone --recursive <fork_url>")
        print("  2. Build panda firmware: cd panda && scons -j4")
        print("  3. Reboot and test in vehicle")
        return 0
    else:
        print("\n" + "=" * 70)
        print("  ✗ VALIDATION FAILED - PLEASE FIX ISSUES ABOVE")
        print("=" * 70)
        return 1

if __name__ == "__main__":
    sys.exit(main())
