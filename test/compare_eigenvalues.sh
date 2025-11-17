#!/bin/bash

echo "=========================================="
echo "Eigenvalue Comparison: θ=9° vs θ=10°"
echo "=========================================="
echo ""

for L in 0 1 2 3; do
    echo "L=$L channel:"
    echo "  θ=9°:"
    grep -A 20 "L =           $L" output_theta9.txt | grep "λ\[  1\]" | head -1
    grep -A 20 "L =           $L" output_theta9.txt | grep "λ\[  2\]" | head -1
    grep -A 20 "L =           $L" output_theta9.txt | grep "λ\[  3\]" | head -1
    echo ""
    echo "  θ=10°:"
    grep -A 20 "L =           $L" output_theta10.txt | grep "λ\[  1\]" | head -1
    grep -A 20 "L =           $L" output_theta10.txt | grep "λ\[  2\]" | head -1
    grep -A 20 "L =           $L" output_theta10.txt | grep "λ\[  3\]" | head -1
    echo ""
    echo "----------------------------------------"
done
