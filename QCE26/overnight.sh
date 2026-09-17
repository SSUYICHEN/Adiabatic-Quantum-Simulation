cat > overnight.sh <<'EOF'
#!/bin/bash
run() {
    local name="$1"; shift
    echo "=== $name START $(date '+%F %T') ===" | tee -a overnight.log
    ( ulimit -v 230686720; /usr/bin/time -v "$@" ) > "logs/${name}.out" 2> "logs/${name}.err"
    echo "=== $name END   $(date '+%F %T') exit=$? ===" | tee -a overnight.log
}

mkdir -p logs

run N6_L40   python aqs_ddsim.py polarization --N 6 --steps 40  --shots 1024
run N8_L40   python aqs_ddsim.py polarization --N 8 --steps 40  --shots 1024


echo "ALL DONE $(date '+%F %T')" | tee -a overnight.log
