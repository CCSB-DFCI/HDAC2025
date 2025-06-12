for d in */ ; do
  if [ -f "$d/quant.sf" ]; then
    echo "===== $d/quant.sf ====="
    head -1 "$d/quant.sf"
    sort -k4,4nr "$d/quant.sf" | head -6
    echo
  fi
done
