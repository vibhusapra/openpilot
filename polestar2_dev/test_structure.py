#!/usr/bin/env python3
"""
Simple structural test for Polestar 2 support without heavy dependencies
"""

import os
import ast
from pathlib import Path

def test_file_structure():
    """Test that all Volvo files exist"""
    print("\n" + "="*60)
    print("TEST 1: File Structure")
    print("="*60)

    base = Path(__file__).parent / "opendbc_repo/opendbc"

    # Check Volvo files
    volvo_files = [
        "car/volvo/__init__.py",
        "car/volvo/carcontroller.py",
        "car/volvo/carstate.py",
        "car/volvo/interface.py",
        "car/volvo/values.py",
        "car/volvo/fingerprints.py",
        "car/volvo/volvocan.py",
        "car/volvo/helpers.py",
    ]

    all_exist = True
    for file in volvo_files:
        path = base / file
        if path.exists():
            size = path.stat().st_size
            print(f"✓ {file} ({size} bytes)")
        else:
            print(f"❌ {file} NOT FOUND")
            all_exist = False

    # Check DBC files
    dbc_files = [
        "dbc/volvo_cma.dbc",
        "dbc/volvo_front_1_cma.dbc",
        "dbc/volvo_front_1_spa.dbc",
        "dbc/volvo_mid_1.dbc",
    ]

    for file in dbc_files:
        path = base / file
        if path.exists():
            size = path.stat().st_size
            print(f"✓ {file} ({size} bytes)")
        else:
            print(f"❌ {file} NOT FOUND")
            all_exist = False

    return all_exist

def test_volvo_registration():
    """Test that Volvo is registered in values.py"""
    print("\n" + "="*60)
    print("TEST 2: Volvo Registration")
    print("="*60)

    values_path = Path(__file__).parent / "opendbc_repo/opendbc/car/values.py"

    with open(values_path) as f:
        content = f.read()

    checks = [
        ("Volvo import", "from opendbc.car.volvo.values import CAR as VOLVO"),
        ("Volvo in Platform", "| VOLVO"),
    ]

    all_good = True
    for check_name, check_str in checks:
        if check_str in content:
            print(f"✓ {check_name}: Found")
        else:
            print(f"❌ {check_name}: NOT FOUND")
            all_good = False

    return all_good

def test_polestar2_definition():
    """Test that POLESTAR_2 is defined"""
    print("\n" + "="*60)
    print("TEST 3: Polestar 2 Definition")
    print("="*60)

    values_path = Path(__file__).parent / "opendbc_repo/opendbc/car/volvo/values.py"

    with open(values_path) as f:
        tree = ast.parse(f.read())

    # Look for CAR class
    polestar_found = False
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == "CAR":
            # Check for POLESTAR_2 attribute
            for item in node.body:
                if isinstance(item, ast.Assign):
                    for target in item.targets:
                        if isinstance(target, ast.Name) and target.id == "POLESTAR_2":
                            polestar_found = True
                            print("✓ POLESTAR_2 defined in CAR class")
                            break

    if not polestar_found:
        print("❌ POLESTAR_2 NOT defined in CAR class")

    # Check fingerprints
    fp_path = Path(__file__).parent / "opendbc_repo/opendbc/car/volvo/fingerprints.py"

    with open(fp_path) as f:
        content = f.read()

    if "CAR.POLESTAR_2" in content:
        print("✓ POLESTAR_2 fingerprint exists")
    else:
        print("❌ POLESTAR_2 fingerprint NOT found")
        polestar_found = False

    return polestar_found

def test_interface_params():
    """Test interface.py has correct function signature"""
    print("\n" + "="*60)
    print("TEST 4: Interface Parameters")
    print("="*60)

    interface_path = Path(__file__).parent / "opendbc_repo/opendbc/car/volvo/interface.py"

    with open(interface_path) as f:
        content = f.read()

    # Check for correct function signature
    if "def _get_params(ret: structs.CarParams, candidate, fingerprint, car_fw, alpha_long, is_release, docs)" in content:
        print("✓ _get_params has correct signature for 0.10.0")
    else:
        print("❌ _get_params signature incorrect")
        return False

    # Check safety model
    if "structs.CarParams.SafetyModel.noOutput" in content:
        print("✓ Uses noOutput safety model (correct)")
    elif "structs.CarParams.SafetyModel.volvo" in content:
        print("❌ References non-existent volvo safety model")
        return False

    # Check brand
    if "ret.brand = 'volvo'" in content:
        print("✓ Brand set to 'volvo'")
    else:
        print("❌ Brand not set correctly")
        return False

    return True

def test_car_specs():
    """Test that Polestar 2 has correct specs"""
    print("\n" + "="*60)
    print("TEST 5: Polestar 2 Specifications")
    print("="*60)

    values_path = Path(__file__).parent / "opendbc_repo/opendbc/car/volvo/values.py"

    with open(values_path) as f:
        content = f.read()

    # Check for Polestar 2 specs
    checks = [
        ("Mass", "mass=2123"),
        ("Wheelbase", "wheelbase=2.735"),
        ("Steer ratio", "steerRatio=15.8"),
        ("CMA platform", "VolvoCMAPlatformConfig"),
    ]

    all_good = True
    # Find the POLESTAR_2 section
    if "POLESTAR_2" in content:
        polestar_section = content[content.index("POLESTAR_2"):]
        # Get just the Polestar 2 definition (until next car or end of class)
        if "VOLVO_S60" in polestar_section:
            polestar_section = polestar_section[:polestar_section.index("VOLVO_S60")]

        for check_name, check_str in checks:
            if check_str in polestar_section:
                print(f"✓ {check_name}: {check_str}")
            else:
                print(f"❌ {check_name}: NOT FOUND")
                all_good = False
    else:
        print("❌ POLESTAR_2 definition not found")
        all_good = False

    return all_good

def main():
    """Run all tests"""
    print("\n" + "="*60)
    print("POLESTAR 2 STRUCTURAL VALIDATION")
    print("="*60)
    print("Testing OpenPilot 0.10.0 + Volvo/Polestar 2 support")

    tests = [
        ("File Structure", test_file_structure),
        ("Volvo Registration", test_volvo_registration),
        ("Polestar 2 Definition", test_polestar2_definition),
        ("Interface Parameters", test_interface_params),
        ("Car Specifications", test_car_specs),
    ]

    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
        except Exception as e:
            print(f"\n❌ Exception in {test_name}: {e}")
            import traceback
            traceback.print_exc()
            result = False
        results.append((test_name, result))

    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)

    for test_name, result in results:
        status = "✓ PASS" if result else "❌ FAIL"
        print(f"  {status}: {test_name}")

    passed = sum(1 for _, r in results if r)
    total = len(results)

    print(f"\nTotal: {passed}/{total} tests passed")

    if passed == total:
        print("\n🎉 ALL STRUCTURAL TESTS PASSED!")
        print("\nThe Polestar 2 integration is structurally ready.")
        print("\nNext steps:")
        print("1. Push to GitHub")
        print("2. Clone on Comma 3")
        print("3. Test with actual hardware")
        return 0
    else:
        print(f"\n⚠️ {total - passed} test(s) failed")
        print("Please review the errors above")
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(main())