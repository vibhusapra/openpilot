#!/usr/bin/env python3
# pylint: skip-file

import argparse
import bisect
import select
from tqdm import tqdm
import sys
import termios
import time
import tty
from collections import defaultdict

import cereal.messaging as messaging
from tools.lib.framereader import FrameReader
from tools.lib.logreader import LogReader
from tools.lib.route import Route

from cereal.visionipc.visionipc_pyx import VisionIpcServer, VisionStreamType  # pylint: disable=no-name-in-module, import-error
from common.transformations.camera import eon_f_frame_size, tici_f_frame_size

IGNORE = ['initData', 'sentinel']

def input_ready():
  return select.select([sys.stdin], [], [], 0) == ([sys.stdin], [], [])


def get_frame(fr, id):
  img = fr.get(id, pix_fmt="rgb24")
  bts = None
  if img is not None:
    img = img[0][:, :, ::-1]  # Convert RGB to BGR, which is what the camera outputs
    img = img.flatten()
    bts = img.tobytes()

  return bts

def replay(route, loop):
  r = Route(route)

  fr = FrameReader(r.camera_paths()[0], readahead=True)
  img = fr.get(0, pix_fmt="rgb24")[0]
  vipc_server = VisionIpcServer("camerad")
  vipc_server.create_buffers(VisionStreamType.VISION_STREAM_RGB_BACK, 4, True, img.shape[1], img.shape[0])
  vipc_server.start_listener()

  socks = {}

  for segment in range(12, 25):
    lr = LogReader(r.log_paths()[segment])
    fr = FrameReader(r.camera_paths()[segment], readahead=True)

    # Build mapping from frameId to segmentId from roadEncodeIdx, type == fullHEVC
    msgs = [m for m in lr if m.which() not in IGNORE]
    msgs = sorted(msgs, key=lambda m: m.logMonoTime)
    frame_idx = {m.roadEncodeIdx.frameId: m.roadEncodeIdx.segmentId for m in msgs if m.which() == 'roadEncodeIdx' and m.roadEncodeIdx.type == 'fullHEVC'}

    frames = {}
    for j in tqdm(frame_idx.values()):
      frames[j] = get_frame(fr, j)

    lag = 0

    for i in tqdm(range(len(msgs) - 1)):
      msg = msgs[i].as_builder()
      next_msg = msgs[i + 1]

      start_time = time.time()
      w = msg.which()

      if w == 'roadCameraState':
        try:
            bts = frames.get(frame_idx[msg.roadCameraState.frameId], None)
            # bts = get_frame(fr, frame_idx[msg.roadCameraState.frameId])
            if bts is not None:
              vipc_server.send(VisionStreamType.VISION_STREAM_RGB_BACK, bts, i, 0, 0)
        except (KeyError, ValueError):
          pass

      if w not in socks:
        socks[w] = messaging.pub_sock(w)

      try:
        if socks[w]:
          socks[w].send(msg.to_bytes())
      except messaging.messaging_pyx.MultiplePublishersError:
        socks[w] = None

      lag += (next_msg.logMonoTime - msg.logMonoTime) / 1e9
      lag -= time.time() - start_time

      dt = max(lag, 0)
      lag -= dt
      time.sleep(dt)

      if lag < -1.0 and i % 1000 == 0:
        print(f"{-lag:.2f} s behind")


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--loop", action='store_true')
  parser.add_argument("route", type=str, help="route name + segment number for offline usage")
  args = parser.parse_args()

  orig_settings = termios.tcgetattr(sys.stdin)
  tty.setcbreak(sys.stdin)

  try:
    replay(args.route, args.loop)
    termios.tcsetattr(sys.stdin, termios.TCSADRAIN, orig_settings)
  except Exception:
    termios.tcsetattr(sys.stdin, termios.TCSADRAIN, orig_settings)
    raise
