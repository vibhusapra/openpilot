#!/usr/bin/env bash

if [ -f env.sh ]; then
  source env.sh
fi

exec ./launch_chffrplus.sh
