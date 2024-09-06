#!/usr/bin/env bash

set -e

# Source message handler functions
source ./ci/scripts/handle_messages.sh

# Function to check dependencies for Python and Fortran linters
function check_dependencies {
    if ! command -v git &> /dev/null; then
        error "git is not installed. Please install it and try again."
        exit 1
    fi

    if ! command -v fortran-linter &> /dev/null; then
        error "fortran-linter is not installed. Please install it and try again."
        exit 1
    fi

    if ! command -v fprettify &> /dev/null; then
        error "fprettify is not installed. Please install it and try again."
        exit 1
    fi

    if ! command -v flake8 &> /dev/null; then
        error "flake8 is not installed. Please install it and try again."
        exit 1
    fi

    if ! command -v pylint &> /dev/null; then
        error "pylint is not installed. Please install it and try again."
        exit 1
    fi
}

# Run Flake8 linter on Python files
function run_flake8 {
    info "Running Flake8 linter"
    
    readarray -t python_files < <(git ls-files '*.py' ':!:tests')
    if [ ${#python_files[@]} -eq 0 ]; then
        info "No Python files to lint"
        return
    fi

    flake8 "${python_files[@]}"
    if [ $? -ne 0 ]; then
        error "Flake8 did not pass"
        exit 2
    fi
    info "Flake8 passed"
}

# Run pylint on non-test Python code
function run_pylint {
    info "Performing pylint, skipping tests"
    msg_format='{C}:{line:3d},{column:2d}: {msg} ({msg_id}: {symbol})'

    readarray -t python_files < <(git ls-files '*.py' ':!:tests/**/*')
    if [ ${#python_files[@]} -eq 0 ]; then
        info "No Python files to lint"
        return
    fi

    pylint "${python_files[@]}" \
           --ignore-paths='.*/tests/.*' \
           --msg-template="$msg_format" \
           --output-format=colorized --rcfile=./.pylintrc
    if [ $? -ne 0 ]; then
        error "Pylint did not pass"
        exit 2
    fi
    info "Pylint passed"
}

# Function to get all Fortran files
get_fortran_files() {
    find . -type f \( -name "*.f90" -o -name "*.f" \) -not -path "./.sandbox/*"
}

# Run Fortran linter
function run_fortran_linter {
    info "Running Fortran linter"
    
    readarray -t fortran_files < <(get_fortran_files)
    if [ ${#fortran_files[@]} -eq 0 ]; then
        info "No Fortran files to lint"
        return
    fi

    for file in "${fortran_files[@]}"; do
        if [ -f "$file" ]; then
            echo "Linting file: $file"
            if ! fortran-linter "$file"; then
                error "Fortran linter did not pass for file: $file"
                exit 2
            fi
        else
            warning "File not found: $file"
        fi
    done
    info "Fortran linter passed"
}

# Run fprettify for Fortran formatting
function run_fprettify {
    info "Running fprettify for Fortran formatting"

    readarray -t fortran_files < <(get_fortran_files)
    if [ ${#fortran_files[@]} -eq 0 ]; then
        info "No Fortran files to format"
        return
    fi

    for file in "${fortran_files[@]}"; do
        if [ -f "$file" ]; then
            echo "Formatting file: $file"
            if ! fprettify "$file"; then
                warning "fprettify failed for file: $file"
            fi
        else
            warning "File not found: $file"
        fi
    done
    info "fprettify completed"
}

# Check dependencies first
check_dependencies

# Run linting and formatting for Python and Fortran
run_flake8
run_pylint
run_fortran_linter
run_fprettify
