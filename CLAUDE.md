# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## About openpilot

openpilot is an operating system for robotics, currently used as an advanced driver assistance system for 300+ supported cars. It's a distributed, real-time autonomous driving system built on a message-passing architecture with hardware-enforced safety.

## Build System

**Primary build tool**: SCons (not Make)

### Common build commands:

```bash
# Full build (C++, Cython, cereal)
scons -j$(nproc)

# Build specific component
scons -j$(nproc) selfdrive/controls/

# Clean build
rm -rf .sconsign.dblite && scons -c && scons -j$(nproc)

# Build with sanitizers
scons --asan       # Address sanitizer
scons --ubsan      # Undefined behavior sanitizer

# Generate compile_commands.json (for IDE support)
scons -j$(nproc)   # Already generated as part of build
```

**Build architecture detection**: Darwin (macOS), x86_64 (Linux PC), aarch64 (ARM64), larch64 (comma device)

## Testing

### Running tests:

```bash
# Run all tests
pytest

# Run specific test file
pytest selfdrive/controls/tests/test_longitudinal_mpc.py

# Run tests for a component
pytest selfdrive/controls/

# Run with specific options
pytest -n auto                    # Parallel execution
pytest -v                         # Verbose
pytest --durations=10            # Show 10 slowest tests
pytest -m "not slow"             # Skip slow tests
pytest -k "test_name_pattern"    # Run tests matching pattern
```

### Test file location:
Tests are co-located with source code in `test_*.py` files or `tests/` directories.

## Linting and Code Quality

```bash
# Run all linters (ruff, mypy, etc.)
pre-commit run --all-files

# Individual linters
ruff check .                  # Python linting
ruff format .                 # Python formatting
mypy selfdrive/               # Type checking
codespell                     # Spell checking

# Run specific pre-commit hook
pre-commit run ruff --all-files
```

**Code style**:
- Python: 2-space indentation, line length 160
- C++: clang-format (2-space indentation)
- Import rules: Use `openpilot.selfdrive`, `openpilot.common`, etc. (not bare `selfdrive` or `common`)

## Running openpilot

```bash
# Set up environment
./launch_env.sh

# Launch full openpilot stack (device only)
./launch_openpilot.sh

# Run individual processes for development
python selfdrive/controls/controlsd.py
python selfdrive/car/card.py

# Run with simulator
cd tools/sim && ./launch_openpilot.sh
```

## Architecture Overview

### Message-Passing System

openpilot uses a **cereal** (Cap'n Proto) message-passing architecture:

- **Message definitions**: `cereal/log.capnp`
- **Messaging library**: `msgq/` (zero-copy shared memory IPC)
- **Pattern**: Pub/Sub with `PubMaster` (publish) and `SubMaster` (subscribe)
- **Camera frames**: Separate `VisionIPC` for zero-copy shared memory

### Core Process Model

**manager.py** (`system/manager/manager.py`) supervises all processes defined in `system/manager/process_config.py`.

**Key processes and their frequencies**:

| Process | Frequency | Purpose |
|---------|-----------|---------|
| `pandad` | - | CAN bus interface, panda firmware communication |
| `card` | 100 Hz | Car interface (CAN parsing/writing) |
| `camerad` | 20 Hz | Camera capture (C++) |
| `modeld` | 20 Hz | Vision model inference (tinygrad) |
| `controlsd` | 100 Hz | Lateral/longitudinal control loops |
| `plannerd` | 20 Hz | Path and speed planning |
| `locationd` | - | Sensor fusion, localization (Kalman filter) |
| `calibrationd` | - | Device calibration |
| `selfdrived` | 10 Hz | State machine, event handling, alerts |
| `loggerd` | - | Data logging and video encoding |
| `radard` | 20 Hz | Radar fusion (when available) |

### Data Flow

**Sense-Plan-Act cycle** at 100 Hz:

```
CAN Bus → pandad → card (CarState) → selfdrived (events) → controlsd (actuators) → card → pandad → CAN Bus
                            ↓
Cameras → camerad → modeld (vision) → plannerd (path/speed) → controlsd
```

### Safety Model

- **Hardware-enforced safety**: panda firmware (STM32) runs safety models from opendbc
- **Safety code**: `panda/board/safety/` (MISRA C compliant)
- Prevents unsafe commands from reaching the car's CAN bus
- Cannot be bypassed by openpilot software

### Directory Structure

**selfdrive/**: Core autonomy stack
- `controls/`: Control loops (lateral/longitudinal control, planning)
  - `controlsd.py`: Main 100Hz control loop
  - `plannerd.py`: Longitudinal planning
  - `lib/latcontrol_*.py`: Lateral control algorithms (PID, torque, angle)
- `car/`: Car interface abstraction
  - `card.py`: Main car process
  - Car-specific code is in `opendbc/car/[manufacturer]/`
- `modeld/`: Vision model inference (GPU/NPU accelerated)
- `locationd/`: Localization and sensor fusion
- `ui/`: On-device UI

**system/**: System services
- `manager/`: Process supervisor and lifecycle
- `loggerd/`: High-performance C++ logging and video encoding
- `hardware/`: Hardware abstraction layer
- `camerad/`: Camera pipeline (C++)
- `updated/`: OTA update system
- `athena/`: Cloud connectivity

**common/**: Shared utilities
- `params.py/cc`: Persistent key-value storage
- `realtime.py`: Rate keeping, real-time scheduling
- `transformations/`: Coordinate frame math
- `swaglog.py`: Logging infrastructure

**opendbc/**: Vehicle interfaces (git submodule at `opendbc_repo/`)
- `opendbc/car/[manufacturer]/`: Per-manufacturer implementations
  - `interface.py`: CarInterface implementation
  - `carstate.py`: CAN message parsing → CarState
  - `carcontroller.py`: CarControl → CAN messages
  - `values.py`: Supported models, fingerprints
  - `radar_interface.py`: Radar parsing (if available)
- `opendbc/dbc/`: CAN database files

**panda/**: CAN hardware interface (git submodule)
- Firmware for comma's panda device (STM32-based)
- `board/safety/`: Safety model implementations
- Python client library

**cereal/**: Message definitions and code generation
- `log.capnp`: All message type definitions
- Generated code for Python and C++

**tools/**: Development tools
- `replay/`: Log replay for testing
- `sim/`: Simulator integration (MetaDrive)
- `cabana/`: CAN analysis tool
- `lib/`: Shared utilities (logreader, route handling)

### Car-Specific Code

All car-specific code lives in **opendbc**:

1. **Fingerprinting**: On startup, `card.py` reads CAN messages and matches them against fingerprints in `opendbc/car/[manufacturer]/values.py`
2. **Interface**: Each manufacturer implements `CarInterface` with methods for state parsing and control
3. **Safety model**: Compiled into panda firmware from `panda/board/safety/safety_[manufacturer].h`
4. **DBC files**: Define CAN message structure in `opendbc/dbc/[manufacturer]_[model].dbc`

To add a new car, create implementations in opendbc following the existing manufacturer patterns.

### Python and C++ Split

**Python**: Most business logic (controls, planning, car interfaces, state machines)
**C++**: Performance-critical code (logging, camera pipeline, messaging infrastructure)
**Cython**: Bridge for performance-critical Python modules (params, some transforms)

### Key Constants

- `DT_CTRL = 0.01` (100 Hz control loop)
- `DT_MDL = 0.05` (20 Hz model inference)
- Main branch: `master`
- Device target: comma 3X (Qualcomm Snapdragon with Adreno GPU)

## Development Tips

- **Logs for debugging**: openpilot logs all messages to `~/.comma/media/0/realdata/` on device (or `/tmp/` on PC)
- **Replay logs**: Use `tools/replay/replay` to replay drives for debugging
- **Parameter storage**: `Params()` from `common/params.py` provides persistent key-value storage
- **Real-time scheduling**: Use `Ratekeeper` from `common/realtime.py` to maintain loop timing
- **Message logging**: All cereal messages are automatically logged by `loggerd`

## Contributing

openpilot prioritizes: **safety > stability > quality > features** (in that order)

- Pull requests against `master` branch
- Must pass CI tests (GitHub Actions)
- Check `docs/CONTRIBUTING.md` for full guidelines
- Simple, well-tested PRs are most likely to be merged
- Large refactors or style changes are generally not accepted
