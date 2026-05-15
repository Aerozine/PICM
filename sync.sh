#!/bin/bash
while true; do
    echo "=== Syncing $(date) ==="
    rsync -Pavz PostPro/video/*.mp4 root@10.0.0.1:/var/www/pi/video/
    echo "=== Done. Next sync in 15min ==="
    sleep 900
done
