class LiveTestingManager:
  TESTING_FILE = "/data/openpilot/live_testing.txt"

  # Valid function names for dot-notation overrides
  VALID_FUNCTIONS = {'lca_5', 'lca', 'lca_2', 'lca_3', 'lca_4'}

  def __init__(self):
    pass  # No state needed - reads file fresh every call

  @staticmethod
  def _parse_value(value_str: str):
    """Auto-detect type: try bool, then int, then float, else keep as string."""
    if value_str in ('True', 'False'):
      return value_str == 'True'
    try:
      return int(value_str)
    except ValueError:
      try:
        return float(value_str)
      except ValueError:
        return value_str

  def load_config(self) -> dict | None:
    """
    Load live testing configuration from file.

    Supports two formats:
    - Global keys: lat_active=True, apply_angle=5.0
    - Function-scoped keys (dot notation): lca_5.lca_turn_bits=128

    Function-scoped keys are grouped into nested dicts with UPPERCASE field names
    (to match DBC signal names).

    Returns:
      - None if file doesn't exist, parse fails, or lat_active=False
      - dict with global keys and function-scoped nested dicts

    Example return values:
      - File doesn't exist → None
      - File contains lat_active=False → None
      - File contains: lat_active=True, lca_5.lca_turn_bits=128, lca.lca_steer=100
        → {'lat_active': True, 'lca_5': {'LCA_TURN_BITS': 128}, 'lca': {'LCA_STEER': 100}}
    """
    try:
      with open(self.TESTING_FILE, 'r') as f:
        lines = f.readlines()

      config = {}

      for line in lines:
        line = line.strip()

        if not line or line.startswith('#'):
          continue

        if '=' not in line:
          continue

        key, value = line.split('=', 1)
        key = key.strip()
        value = value.strip()

        # Strip inline comments
        if '#' in value:
          value = value.split('#', 1)[0].strip()

        # Check for dot notation (function-scoped override)
        if '.' in key:
          func_name, field_name = key.split('.', 1)
          if func_name not in self.VALID_FUNCTIONS:
            continue
          if func_name not in config:
            config[func_name] = {}
          # Convert field name to UPPERCASE (DBC signals are uppercase)
          config[func_name][field_name.upper()] = self._parse_value(value)
        else:
          # Global keys (lat_active, apply_angle)
          if key == 'lat_active':
            if value == 'True':
              config['lat_active'] = True
            elif value == 'False':
              config['lat_active'] = False
          elif key == 'apply_angle':
            try:
              config['apply_angle'] = float(value)
            except ValueError:
              continue

      # If lat_active is explicitly False, return None
      if 'lat_active' in config and config['lat_active'] is False:
        return None

      if not config:
        return None

      return config

    except FileNotFoundError:
      return None
    except Exception:
      return None
