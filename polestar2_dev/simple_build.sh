#!/bin/bash

# Simple build test - Direct and minimal

echo "=========================================="
echo "  Simple Build Test - Polestar 2"
echo "=========================================="
echo

echo "Working directory: $(pwd)"
echo "Python: $(which python3)"
echo "Scons: $(which scons)"
echo

echo "Cleaning..."
scons -c -s 2>/dev/null
rm -rf .sconsign.dblite cereal/gen 2>/dev/null
echo

echo "Building with verbose output..."
echo "=========================================="

# Direct build command with error output
scons -u -j1 2>&1 | head -100

echo
echo "=========================================="
echo "Checking build results..."

if [ -f "cereal/libcereal.a" ]; then
    echo "✓ cereal/libcereal.a built"
else
    echo "✗ cereal/libcereal.a missing"
fi

if [ -f "msgq_repo/libmsgq.a" ]; then
    echo "✓ msgq_repo/libmsgq.a built"
else
    echo "✗ msgq_repo/libmsgq.a missing"
fi

if [ -f "common/libcommon.a" ]; then
    echo "✓ common/libcommon.a built"
else
    echo "✗ common/libcommon.a missing"
fi

echo "=========================================="