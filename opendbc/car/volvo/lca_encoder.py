"""
LCA Target Angle Encoder – Asymmetric Left/Right Model
======================================================

Assumptions from empirical tests and desired behavior:

- SCALE is the same for both directions:

    SCALE = 0.05596 degrees per count

LEFT (angle > 0):
-----------------
- Byte 6 starts at 128, byte 7 goes 0→255
- When byte 7 exceeds 255, it wraps to 0 and byte 6 increments to 129, etc.

- Mapping using left_counts:
    left_counts = (byte6 - 128) * 256 + byte7
    angle_deg = left_counts * SCALE
    - So left_counts = round(angle_deg / SCALE)

- Near 0°, we treat (128, 0) as the canonical "neutral" command.

RIGHT (angle < 0, theoretical/ideal behavior):
----------------------------------------------
- At small right angles, zero point is:
    byte6 = 255, byte7 = 255  -> angle ≈ 0°

- Increasing right angle (more negative) should:
    - Decrease byte7 from 255 → 254 → ... → 0
    - Then wrap:
        (byte6, byte7) = (255,   0)  -> next step ->
        (byte6, byte7) = (254, 255)
      and continue with byte6 decreasing as needed.

- We model this as a simple "right counts" index:

    right_counts = (255 - byte6) * 256 + (255 - byte7)   # 0 at (255,255)

    angle_deg = - right_counts * SCALE

SPECIAL:
--------
- (byte6, byte7) = (186, 0) is the "inactive / no command" state and
  should NOT be used to request 0°; encode_inactive() returns this and
  decode() returns None for it.

- For angle >= 0, encoder uses LEFT encoding; for angle < 0, uses RIGHT encoding.
"""

SCALE = 0.05596  # degrees per count (symmetric magnitude)


class LCATargetAngleEncoder:
    """Encode/decode target steering angles to LCA_5 bytes 6+7."""

    @staticmethod
    def encode(target_angle_deg: float):
        """
        Encode target steering angle to (byte6, byte7) for ACTIVE LCA.

        Args:
            target_angle_deg: Target angle in degrees.
                              Positive = LEFT, Negative = RIGHT, 0 = straight.

        Returns:
            (byte6, byte7): Tuple of bytes for LCA_5 message.
        """

        # LEFT: byte6 starts at 128, byte7 0→255, then byte6 increments and byte7 wraps
        if target_angle_deg >= 0:
            # left_counts = (byte6 - 128) * 256 + byte7
            left_counts = int(round(target_angle_deg / SCALE))
            left_counts = max(0, min(0xFFFF, left_counts))  # Allow full range

            # Decompose left_counts into (byte6, byte7)
            k = left_counts // 256
            r = left_counts % 256
            byte6 = 128 + k
            byte7 = r

            # Avoid inactive pattern accidentally
            if byte6 == 186 and byte7 == 0:
                byte7 = 1
            return (byte6, byte7)

        # RIGHT: use right_counts with (255,255) as zero
        # angle_deg = - right_counts * SCALE
        # => right_counts = -angle_deg / SCALE
        right_counts = int(round(-target_angle_deg / SCALE))
        if right_counts < 0:
            right_counts = 0
        if right_counts > 0xFFFF:
            right_counts = 0xFFFF

        # Decompose right_counts into (byte6, byte7)
        # right_counts = (255 - byte6) * 256 + (255 - byte7)
        # Let k = 255 - byte6, r = 255 - byte7:
        #   right_counts = k * 256 + r
        k = right_counts // 256
        r = right_counts % 256

        byte6 = 255 - k
        byte7 = 255 - r

        # Avoid inactive pattern accidentally
        if byte6 == 186 and byte7 == 0:
            # Nudge one step further right
            right_counts = min(right_counts + 1, 0xFFFF)
            k = right_counts // 256
            r = right_counts % 256
            byte6 = 255 - k
            byte7 = 255 - r

        return (byte6, byte7)

    @staticmethod
    def encode_inactive():
        """
        Return the byte values for when LCA is inactive/standby.

        Returns:
            (186, 0): The "no command" state.
        """
        return (186, 0)

    @staticmethod
    def encode_lca_steer(target_angle_deg: float) -> int:
        """
        Encode target angle to 1-byte LCA_STEER value (0-255) for LCA message (0x58).
        Derived from the same encoding as LCA_5_STEER for consistency.

        Args:
            target_angle_deg: Target angle in degrees.
                              Positive = LEFT, Negative = RIGHT, 0 = straight.

        Returns:
            int: LCA_STEER value (byte7 from the 2-byte encoding)
        """
        byte6, byte7 = LCATargetAngleEncoder.encode(target_angle_deg)
        return byte7

    @staticmethod
    def decode(byte6: int, byte7: int):
        """
        Decode (byte6, byte7) to target steering angle.

        Args:
            byte6, byte7: LCA_5 message bytes.

        Returns:
            angle_deg: Target steering angle in degrees, or
                       None if this is the inactive/no-command state.
        """
        # Inactive / no command
        if byte6 == 186 and byte7 == 0:
            return None

        # LEFT side: byte6 >= 128 (starts at 128, increments as angle increases)
        if byte6 >= 128:
            # left_counts = (byte6 - 128) * 256 + byte7
            left_counts = (byte6 - 128) * 256 + byte7
            return SCALE * left_counts

        # RIGHT side region: byte6 < 128 (starts at 255, decrements as angle increases)
        # Note: byte6 values 255, 254, ... wrap around, so byte6 < 128 means we've
        # wrapped past the inactive zone. Use right_counts formula.
        right_counts = (255 - byte6) * 256 + (255 - byte7)
        return - right_counts * SCALE
