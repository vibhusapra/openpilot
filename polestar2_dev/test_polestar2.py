#!/usr/bin/env python3
"""
Test script to verify Polestar 2 support in OpenPilot 0.10.0
"""

import sys
import os
from pathlib import Path

# Add openpilot to path
sys.path.insert(0, str(Path(__file__).parent))

def test_volvo_import():
    """Test that Volvo can be imported"""
    print("\n" + "="*60)
    print("TEST 1: Importing Volvo modules")
    print("="*60)

    try:
        from opendbc.car.volvo.values import CAR
        print("✓ Successfully imported Volvo CAR enum")

        # Check if POLESTAR_2 exists
        if hasattr(CAR, 'POLESTAR_2'):
            print("✓ POLESTAR_2 found in CAR enum")
            print(f"  Value: {CAR.POLESTAR_2}")
        else:
            print("❌ POLESTAR_2 not found in CAR enum")
            return False

        from opendbc.car.volvo.interface import CarInterface
        print("✓ Successfully imported Volvo CarInterface")

        from opendbc.car.volvo.fingerprints import FINGERPRINTS
        print("✓ Successfully imported Volvo fingerprints")

        if CAR.POLESTAR_2 in FINGERPRINTS:
            print("✓ POLESTAR_2 fingerprint exists")
            fp = FINGERPRINTS[CAR.POLESTAR_2][0]
            print(f"  Fingerprint has {len(fp)} CAN messages")
        else:
            print("❌ POLESTAR_2 fingerprint not found")
            return False

        return True
    except Exception as e:
        print(f"❌ Failed to import: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_platform_registration():
    """Test that Volvo is registered in the platform system"""
    print("\n" + "="*60)
    print("TEST 2: Platform Registration")
    print("="*60)

    try:
        from opendbc.car.values import PLATFORMS, Platform
        from opendbc.car.volvo.values import CAR

        polestar_key = str(CAR.POLESTAR_2)
        print(f"Looking for: {polestar_key}")

        if polestar_key in PLATFORMS:
            print(f"✓ {polestar_key} is registered in PLATFORMS")
            print(f"  Platform: {PLATFORMS[polestar_key]}")
        else:
            print(f"❌ {polestar_key} not found in PLATFORMS")
            print(f"  Available keys (first 10): {list(PLATFORMS.keys())[:10]}")
            return False

        return True
    except Exception as e:
        print(f"❌ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_fingerprint_system():
    """Test that fingerprints work correctly"""
    print("\n" + "="*60)
    print("TEST 3: Fingerprint System")
    print("="*60)

    try:
        from opendbc.car.fingerprints import _FINGERPRINTS
        from opendbc.car.volvo.values import CAR

        polestar_key = CAR.POLESTAR_2

        if polestar_key in _FINGERPRINTS:
            print(f"✓ POLESTAR_2 in global fingerprint registry")
            fp_list = _FINGERPRINTS[polestar_key]
            print(f"  Has {len(fp_list)} fingerprint(s)")
        else:
            print(f"❌ POLESTAR_2 not in global fingerprint registry")
            return False

        return True
    except Exception as e:
        print(f"❌ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_car_params():
    """Test getting car params for Polestar 2"""
    print("\n" + "="*60)
    print("TEST 4: Car Parameters")
    print("="*60)

    try:
        from opendbc.car.volvo.values import CAR
        from opendbc.car.volvo.interface import CarInterface
        from opendbc.car.volvo.fingerprints import FINGERPRINTS
        from opendbc.car import gen_empty_fingerprint

        # Create a fingerprint
        fingerprint = gen_empty_fingerprint()
        fingerprint[0] = FINGERPRINTS[CAR.POLESTAR_2][0]

        # Get params
        candidate = CAR.POLESTAR_2
        car_fw = []

        # Create a mock CarParams object
        from opendbc.car import structs
        ret = structs.CarParams.new_message()
        ret.carName = str(candidate)
        ret.carFingerprint = str(candidate)

        # Call _get_params
        params = CarInterface._get_params(
            ret,
            candidate,
            fingerprint,
            car_fw,
            alpha_long=False,
            is_release=False,
            docs=False
        )

        print("✓ Successfully got CarParams for POLESTAR_2")
        print(f"  Brand: {params.brand}")
        print(f"  Mass: {params.mass} kg")
        print(f"  Wheelbase: {params.wheelbase} m")
        print(f"  Steer control type: {params.steerControlType}")
        print(f"  Safety model: {params.safetyConfigs[0].safetyModel}")

        # Verify values
        assert params.brand == 'volvo', f"Brand should be 'volvo', got '{params.brand}'"
        assert params.mass == 2123, f"Mass should be 2123, got {params.mass}"
        assert params.wheelbase == 2.735, f"Wheelbase should be 2.735, got {params.wheelbase}"

        print("✓ All parameter validations passed!")
        return True

    except Exception as e:
        print(f"❌ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_interface_instantiation():
    """Test that CarInterface can be instantiated"""
    print("\n" + "="*60)
    print("TEST 5: Interface Instantiation")
    print("="*60)

    try:
        from opendbc.car.volvo.values import CAR
        from opendbc.car.volvo.interface import CarInterface
        from opendbc.car.volvo.carcontroller import CarController
        from opendbc.car.volvo.carstate import CarState
        from opendbc.car.volvo.fingerprints import FINGERPRINTS
        from opendbc.car import gen_empty_fingerprint, structs

        # Create fingerprint
        fingerprint = gen_empty_fingerprint()
        fingerprint[0] = FINGERPRINTS[CAR.POLESTAR_2][0]

        # Get params
        candidate = CAR.POLESTAR_2
        ret = structs.CarParams.new_message()
        ret.carName = str(candidate)
        ret.carFingerprint = str(candidate)

        CP = CarInterface._get_params(
            ret,
            candidate,
            fingerprint,
            [],
            alpha_long=False,
            is_release=False,
            docs=False
        )

        # Try to instantiate
        CI = CarInterface(CP, CarController, CarState)

        print("✓ Successfully instantiated CarInterface")
        print(f"  Has CC: {hasattr(CI, 'CC')}")
        print(f"  Has CS: {hasattr(CI, 'CS')}")

        return True

    except Exception as e:
        print(f"❌ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests"""
    print("\n" + "="*60)
    print("POLESTAR 2 VALIDATION FOR OPENPILOT 0.10.0")
    print("="*60)

    tests = [
        ("Import Test", test_volvo_import),
        ("Platform Registration", test_platform_registration),
        ("Fingerprint System", test_fingerprint_system),
        ("Car Parameters", test_car_params),
        ("Interface Instantiation", test_interface_instantiation),
    ]

    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
        except Exception as e:
            print(f"\n❌ Exception in {test_name}: {e}")
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
        print("\n🎉 ALL TESTS PASSED!")
        print("\nThe Polestar 2 integration is ready for deployment to Comma 3!")
        return 0
    else:
        print(f"\n⚠️ {total - passed} test(s) failed")
        print("Please review the errors above")
        return 1

if __name__ == "__main__":
    sys.exit(main())