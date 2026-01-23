#pragma once

#include "opendbc/safety/declarations.h"

// Volvo CMA platform CAN message addresses
#define VOLVO_LCA_STEER           0x58U    // TX from VCU1 to PSCM, LCA steering command (0x58)
#define VOLVO_BUS1_SPEED          0x70U   // RX from BCM, vehicle speed (BUS1_SPEED)
#define VOLVO_LCA_2               0x69U   // RX from BCM, brake pedal, cruise state
#define VOLVO_SAS                 0x55U    // RX from SAS, steering angle sensor
#define VOLVO_PSCM                0x16U    // RX from PSCM, driver steering input
#define VOLVO_GEAR_POSITION       0x80U   // RX from transmission, gear position
#define VOLVO_ECM_1               0x250U   // RX from ECM, accelerator pedal position (0x250)
#define VOLVO_DRIVER_INPUT        0x15U
#define VOLVO_BUS1_CRUISE_CONTROL 0x340U   // RX from BCM, cruise control state (BUS1_CRUISE_CONTROL)
#define VOLVO_LCA_3               0x57U   // TX from VCU1 to PSCM
#define VOLVO_LCA_5               0x67U   // TX LCA_5 message (formerly SPEED_1, contains wheel speeds and LCA signals)
#define VOLVO_SPEED               0x60U   // RX/TX SPEED message
#define VOLVO_SPEED_2             0x68U   // RX
#define VOLVO_0x1a                0x1aU   // RX
#define VOLVO_EGSM                0x45U   // RX from EGSM
#define VOLVO_PSCM_RELATED        0x17U   // RX from PSCM, related messages
#define VOLVO_LCA_4               0x90U   // TX LCA_4 message (PA status spoofing)
#define VOLVO_LCA_6               0x97U   // TX LCA_6 message
#define VOLVO_LCA_7               0x92U   // TX LCA_7 message


// CAN bus definitions for Volvo CMA platform
// Using same naming as carstate.py for consistency: main, pt, party
#define VOLVO_MAIN_BUS    0U  // Bus.main - VCU1 car side
#define VOLVO_PT_BUS      1U  // Bus.pt - VCU1 ECM side (where ECM is)
#define VOLVO_PARTY_BUS   2U  // Bus.party - VCU PSCM/BCM2 side (BCM2, SAS, EGSM, PSCM, where LCA is sent to)

static void volvo_rx_hook(const CANPacket_t *msg) {
  // Basic vehicle state monitoring - very relaxed implementation

  // Main bus (bus 0) messages
  if (msg->bus == VOLVO_MAIN_BUS) {
    // Update brake pedal and cruise state from BCM2
    if (msg->addr == VOLVO_LCA_2) {
      // DBC: SG_ BRAKE_PEDAL_PRESSED_A : 47|1@0+ (-1,1) - inverted in DBC, so we invert raw bit
      // DBC: SG_ BRAKE_PEDAL_PRESSED_B : 46|1@0+ (1,0) - not inverted
      //bool brake_a = !((msg->data[5] >> 7) & 1U); // Raw bit, active low (DBC inverts it)
      bool brake_b = (msg->data[5] >> 6) & 1U; // Raw bit, active high
      //brake_pressed = brake_a || brake_b;
      brake_pressed = brake_b;
    }
  }

  // PT bus (bus 1) messages
  if (msg->bus == VOLVO_PT_BUS) {
    if (msg->addr == VOLVO_ECM_1) {
      // Gas pedal position - ACCELERATOR_PEDAL_POS
      // DBC: SG_ ACCELERATOR_PEDAL_POS : 31|8@0+ (1,0) [0|255]
      // carstate.py: > 20+1 (20 baseline + 1 tolerance)
      uint8_t gas_pedal_position = msg->data[3];
      gas_pressed = gas_pedal_position > 20+1; // Match carstate.py tolerance
    }

    // Update vehicle speed from BUS1_SPEED
    if (msg->addr == VOLVO_BUS1_SPEED) {
      // Signal: BUS1_SPEED (0.015625 m/s per bit)
      // DBC: SG_ BUS1_SPEED : 23|16@0+ (0.015625,0) [0|65535] "m/s" XXX
      uint16_t speed_raw = ((msg->data[2] & 0xFFU) << 8) | msg->data[3];
      vehicle_moving = (speed_raw * 0.01886) > 0.1; // > 0.1 m/s
      UPDATE_VEHICLE_SPEED(speed_raw * 0.01886);
    }

    if (msg->addr == VOLVO_BUS1_CRUISE_CONTROL) {
      // SG_ CRUISE_CONTROL_ENABLED : 56|1@0+ (1,0) [0|1] "" XXX
      // SG_ CRUISE_CONTROL_ENABLED_IDLE_TRAFFIC : 57|1@0+ (1,0) [0|1] "" XXX
      bool cruise_enabled = ((msg->data[7] & 1U) || (msg->data[7] & 2U));
      pcm_cruise_check(cruise_enabled);
    }
  }

  // Party bus (bus 2) messages - BCM2, SAS, PSCM, EGSM
  if (msg->bus == VOLVO_PARTY_BUS) {

    // Update steering angle from SAS - Deprecated: Use PSCM instead
    /*if (msg->addr == VOLVO_SAS) {
      // DBC: SG_ SAS_ANGLE_SENSOR : 6|15@0- (-0.05596,0)
      // carstate.py uses SAS (not PSCM) for steering angle
      // Bit position 6, 15 bits, signed, little endian
      int angle_raw = ((msg->data[0] & 0x7FU) << 8) | msg->data[1];
      if (msg->data[0] & 0x80U) {
        angle_raw = -angle_raw;
      }
      update_sample(&angle_meas, angle_raw);
    }*/

    // Update driver steering input from PSCM
    if (msg->addr == VOLVO_PSCM) {
      // DBC: SG_ DRIVER_INPUT_DEVIATION : 47|8@0- (1,0)
      // carstate.py uses abs(DRIVER_INPUT_DEVIATION) for steering torque and pressed detection
      // Bit position 47, 8 bits, signed
      int driver_input = msg->data[5];
      update_sample(&torque_driver, driver_input);
      // The following might not be right, probably needs to be re-made, tests don't catch any inconsistencies with CarState (no checks are made)
      /*uint8_t byte0_masked = msg->data[0] & 0x7F;   // keep only lower 7 bits
      uint16_t combined = ((uint16_t)byte0_masked << 8) | msg->data[1];
      update_sample(&angle_meas, -(int16_t)combined); // openpilot excepts a right turn to be negative*/
    }
  }
}

static bool volvo_tx_hook(const CANPacket_t *msg) {
  bool tx = true;

  // Very relaxed safety policy - only basic frame ID checks
  if (msg->addr == VOLVO_LCA_STEER) {
    // LCA message flows: VCU1 (main bus) -> PSCM (party bus)
    // We're acting as VCU1, so we send LCA message to party bus (bus 2)
    if (msg->bus != VOLVO_PARTY_BUS) {
      tx = false;  // Wrong bus
    }
  }

  if (msg->addr == VOLVO_PSCM) {
    // PSCM message: we relay from party bus (bus 2) to main bus (bus 0)
    // So we TX on main bus (bus 0)
    if (msg->bus != VOLVO_MAIN_BUS) {
      tx = false;  // Wrong bus
    }
  }

  if (msg->addr == VOLVO_DRIVER_INPUT) {
    // Driver input message: we relay from party bus (bus 2) to main bus (bus 0)
    // So we TX on main bus (bus 0)
    if (msg->bus != VOLVO_MAIN_BUS) {
      tx = false;  // Wrong bus
    }
  }

  if (msg->addr == VOLVO_SAS) {
    // SAS message: we relay from party bus (bus 2) to main bus (bus 0)
    // So we TX on main bus (bus 0)
    if (msg->bus != VOLVO_MAIN_BUS) {
      tx = false;  // Wrong bus
    }
  }

  if (msg->addr == VOLVO_LCA_2) {
    // LCA_2 -> PSCM
    if (msg->bus != VOLVO_PARTY_BUS) {
      tx = false;  // Wrong bus
    }
  }

  return tx;
}

static safety_config volvo_init(uint16_t param) {
  SAFETY_UNUSED(param);

  // Define allowed TX messages - very permissive
  static const CanMsg VOLVO_TX_MSGS[] = {
    {VOLVO_LCA_STEER, VOLVO_PARTY_BUS, 8, .check_relay = true},  // LCA steering command to party bus
    {VOLVO_PSCM, VOLVO_MAIN_BUS, 8, .check_relay = true},  // PSCM message sent to main bus (relay from party bus)
    {VOLVO_LCA_3, VOLVO_PARTY_BUS, 8, .check_relay = true},  // LCA_3 message sent to party bus
    {VOLVO_LCA_2, VOLVO_PARTY_BUS, 8, .check_relay = true},  // LCA_2 message sent to party bus (spoof PILOT_ASSIST_ENGAGED for PSCM)
    {VOLVO_LCA_4, VOLVO_PARTY_BUS, 8, .check_relay = true},  // LCA_4 message sent to party bus (spoof LCA_ENABLE for PA state)
    {VOLVO_LCA_5, VOLVO_PARTY_BUS, 8, .check_relay = true},  // LCA_5 message sent to party bus (wheel speeds + LCA signals)
    {VOLVO_LCA_6, VOLVO_PARTY_BUS, 8, .check_relay = true},  // LCA_6 message sent to party bus
    {VOLVO_LCA_7, VOLVO_PARTY_BUS, 8, .check_relay = true},  // LCA_7 message sent to party bus
    //{VOLVO_SPEED, VOLVO_PARTY_BUS, 8, .check_relay = true},  // SPEED message sent to main bus
    //{VOLVO_SPEED_2, VOLVO_PARTY_BUS, 8, .check_relay = true},  // SPEED_2 message sent to main bus
    //{VOLVO_0x1a, VOLVO_PARTY_BUS, 8, .check_relay = true},  // 0x1a message sent to main bus
    //{VOLVO_GEAR_POSITION, VOLVO_PARTY_BUS, 8, .check_relay = true},  // GEAR_POSITION message sent from main to party bus
    //{VOLVO_EGSM, VOLVO_MAIN_BUS, 8, .check_relay = true},  // EGSM message sent from party to main bus
    {VOLVO_PSCM_RELATED, VOLVO_MAIN_BUS, 8, .check_relay = true},  // PSCM_RELATED message sent to party bus
  };

  // Define RX checks - include all messages present in route
  static RxCheck volvo_rx_checks[] = {
    {.msg = {{VOLVO_GEAR_POSITION, VOLVO_MAIN_BUS, 8, 40U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},  // Using bus 2
    {.msg = {{VOLVO_BUS1_SPEED, VOLVO_PT_BUS, 8, 100U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_LCA_2, VOLVO_MAIN_BUS, 8, 50U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_LCA_4, VOLVO_MAIN_BUS, 8, 29U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_LCA_6, VOLVO_MAIN_BUS, 8, 25U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_LCA_7, VOLVO_MAIN_BUS, 8, 29U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_SAS, VOLVO_PARTY_BUS, 8, 100U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_PSCM, VOLVO_PARTY_BUS, 8, 100U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_DRIVER_INPUT, VOLVO_PARTY_BUS, 8, 100U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_ECM_1, VOLVO_PT_BUS, 8, 17U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_LCA_STEER, VOLVO_MAIN_BUS, 8, 100U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_BUS1_CRUISE_CONTROL, VOLVO_PT_BUS, 8, 10U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_LCA_3, VOLVO_MAIN_BUS, 8, 67U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_LCA_5, VOLVO_MAIN_BUS, 8, 50U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_SPEED, VOLVO_MAIN_BUS, 8, 50U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_SPEED_2, VOLVO_MAIN_BUS, 8, 50U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_EGSM, VOLVO_PARTY_BUS, 8, 100U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
    {.msg = {{VOLVO_PSCM_RELATED, VOLVO_PARTY_BUS, 8, 100U, .ignore_checksum = true, .ignore_counter = true, .ignore_quality_flag = true}, { 0 }, { 0 }}},
  };

  return BUILD_SAFETY_CFG(volvo_rx_checks, VOLVO_TX_MSGS);
}

const safety_hooks volvo_hooks = {
  .init = volvo_init,
  .rx = volvo_rx_hook,
  .tx = volvo_tx_hook,
  // No custom fwd hook - stock LCA always blocked by .check_relay = true
};
