#!/bin/bash

# A script to test the python ports of the perl scripts.

# Exit on error
set -e

# Activate python environment
source .venv/bin/activate
export PYTHONPATH=.

# Define directories
INPUT_DIR="test_data/input"
PERL_OUTPUT_DIR="test_data/perl_output"
PYTHON_OUTPUT_DIR="test_data/python_output"

# Clean up previous results
rm -f $PERL_OUTPUT_DIR/*
rm -f $PYTHON_OUTPUT_DIR/*

# A function to run and compare a script
test_script() {
    local script_name=$1
    shift
    local args="$@"

    echo "Testing $script_name..."

    # Run Perl script
    perl "scripts/${script_name}.pl" $args > "$PERL_OUTPUT_DIR/${script_name}.out"

    # Run Python script
    python "scripts/${script_name}.py" $args > "$PYTHON_OUTPUT_DIR/${script_name}.out"

    # Compare outputs
    diff -wB <(sort "$PERL_OUTPUT_DIR/${script_name}.out") <(sort "$PYTHON_OUTPUT_DIR/${script_name}.out")
    echo "$script_name: OK"
    echo "--------------------"
}

test_script_stdin() {
    local script_name=$1
    local input_file=$2
    shift 2
    local args="$@"

    echo "Testing $script_name..."

    # Run Perl script
    cat "$input_file" | perl "scripts/${script_name}.pl" $args > "$PERL_OUTPUT_DIR/${script_name}.out"

    # Run Python script
    cat "$input_file" | python "scripts/${script_name}.py" $args > "$PYTHON_OUTPUT_DIR/${script_name}.out"

    # Compare outputs
    diff -wB <(sort "$PERL_OUTPUT_DIR/${script_name}.out") <(sort "$PYTHON_OUTPUT_DIR/${script_name}.out")
    echo "$script_name: OK"
    echo "--------------------"
}


# --- Tests for scripts with text output ---

test_script_stdin clstr2txt "$INPUT_DIR/test.clstr"
test_script_stdin clstr_renumber "$INPUT_DIR/test.clstr"
test_script_stdin clstr_rep "$INPUT_DIR/test.clstr"
test_script_stdin clstr_quality_eval "$INPUT_DIR/test.clstr"
test_script_stdin clstr_quality_eval_by_link "$INPUT_DIR/test.clstr"

test_script clstr_reduce "$INPUT_DIR/test.clstr" "1,2-3" 2
test_script clstr_reps_faa_rev "$INPUT_DIR/test.clstr" "$INPUT_DIR/test.fasta" 2

# --- Tests for scripts with more complex I/O ---

# clstr_list and clstr_list_sort produce binary files, which are hard to diff directly.
# We will test them by their effect: create a store file, sort it, and then we need a way to view the content.
# For now, we will just run them to see if they crash.

echo "Testing clstr_list..."
perl scripts/clstr_list.pl "$INPUT_DIR/test.clstr" "$PERL_OUTPUT_DIR/clstr_list.store"
python scripts/clstr_list.py "$INPUT_DIR/test.clstr" "$PYTHON_OUTPUT_DIR/clstr_list.store"
echo "clstr_list: OK (ran without crash)"
echo "--------------------"

echo "Testing clstr_list_sort..."
perl scripts/clstr_list_sort.pl "$PERL_OUTPUT_DIR/clstr_list.store" "$PERL_OUTPUT_DIR/clstr_list_sort_no.store" "no"
python scripts/clstr_list_sort.py "$PYTHON_OUTPUT_DIR/clstr_list.store" "$PYTHON_OUTPUT_DIR/clstr_list_sort_no.store" "no"
echo "clstr_list_sort (sort by no): OK (ran without crash)"

perl scripts/clstr_list_sort.pl "$PERL_OUTPUT_DIR/clstr_list.store" "$PERL_OUTPUT_DIR/clstr_list_sort_len.store" "len"
python scripts/clstr_list_sort.py "$PYTHON_OUTPUT_DIR/clstr_list.store" "$PYTHON_OUTPUT_DIR/clstr_list_sort_len.store" "len"
echo "clstr_list_sort (sort by len): OK (ran without crash)"

perl scripts/clstr_list_sort.pl "$PERL_OUTPUT_DIR/clstr_list.store" "$PERL_OUTPUT_DIR/clstr_list_sort_des.store" "des"
python scripts/clstr_list_sort.py "$PYTHON_OUTPUT_DIR/clstr_list.store" "$PYTHON_OUTPUT_DIR/clstr_list_sort_des.store" "des"
echo "clstr_list_sort (sort by des): OK (ran without crash)"
echo "--------------------"

# clstr_merge_noorder requires multiple clstr files. Let's create them.
echo "Creating data for clstr_merge_noorder"
cat > test_data/input/merge_master.clstr << EOL
>Cluster 0
0	100aa, >seq1... *
1	100aa, >seq2... at 100.00%
EOL

cat > test_data/input/merge_slave1.clstr << EOL
>Cluster 5
0	100aa, >seq1... *
1	90aa, >seq3... at 90.00%
EOL

test_script clstr_merge_noorder "test_data/input/merge_master.clstr" "test_data/input/merge_slave1.clstr"

echo "All tests passed!"
