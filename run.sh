for d in ./short-cadence/Kepler/*; do
  if [ -d "$d" ]; then
    for d in "$d"; do
      if [ -f "$d" ]; then
        ./reader "$d"
      fi
      break
    done
  fi
done