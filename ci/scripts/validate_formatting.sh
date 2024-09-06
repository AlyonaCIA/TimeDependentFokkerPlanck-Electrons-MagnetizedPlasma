#!/usr/bin/env bash

set -e

# Source message handler functions
source ci/scripts/handle_messages.sh

# --- Fortran Section ---

# Function to get all Fortran files
get_fortran_files() {
    find . -type f \( -name "*.f90" -o -name "*.f" \) -not -path "./.sandbox/*"
}

# Function to check formatting with fprettify
check_fprettify() {
    info "Checking formatting with fprettify for Fortran files"
    
    readarray -t fortran_files < <(get_fortran_files)
    if [ ${#fortran_files[@]} -eq 0 ]; then
        info "No Fortran files to check"
        return
    fi

    echo "Processing Fortran files: ${fortran_files[@]}"
    for file in "${fortran_files[@]}"; do
        if ! diff -q <(fprettify --silent "$file") "$file"; then
            error "fprettify did not pass for file: $file"
            exit 2
        fi
    done
    info "fprettify passed for all Fortran files!"
}

# Function to check linting with fortran-linter
check_fortran_linter() {
    info "Checking linting with fortran-linter"

    readarray -t fortran_files < <(get_fortran_files)
    if [ ${#fortran_files[@]} -eq 0 ]; then
        info "No Fortran files to lint"
        return
    fi

    echo "Processing Fortran files: ${fortran_files[@]}"
    for file in "${fortran_files[@]}"; do
        if ! fortran-linter "$file"; then
            error "fortran-linter did not pass for file: $file"
            exit 2
        fi
    done
    info "fortran-linter passed for all Fortran files!"
}

# --- Python Section ---

# Function to get all Python files
get_python_files() {
    find . -type f -name "*.py" -not -path "./.sandbox/*"
}

# Function to check formatting with autopep8
check_autopep8() {
    info "Checking formatting with autopep8 for Python files"

    readarray -t python_files < <(get_python_files)
    if [ ${#python_files[@]} -eq 0 ]; then
        info "No Python files to check"
        return
    fi

    echo "Processing Python files: ${python_files[@]}"
    if autopep8 --diff "${python_files[@]}" --exit-code; then
        info "Autopep8 passed!"
    else
        error "Autopep8 did not pass"
        exit 2
    fi
}

# Function to check doc formats with docformatter
check_docformatter() {
    info "Checking doc formats with docformatter"

    readarray -t python_files < <(get_python_files)
    if [ ${#python_files[@]} -eq 0 ]; then
        info "No Python files to check"
        return
    fi

    echo "Processing Python files: ${python_files[@]}"
    if docformatter --wrap-summaries 88 --wrap-descriptions 88 --check "${python_files[@]}"; then
        info "Docformatter passed!"
    else
        error "Docformatter did not pass"
        exit 2
    fi
}

# Function to check for unused imports and variables with autoflake
check_autoflake() {
    info "Checking for unused imports and variables with autoflake"

    readarray -t python_files < <(get_python_files)
    if [ ${#python_files[@]} -eq 0 ]; then
        info "No Python files to check"
        return
    fi

    echo "Processing Python files: ${python_files[@]}"
    if autoflake --check --remove-all-unused-imports --ignore-init-module-imports "${python_files[@]}"; then
        info "Autoflake passed!"
    else
        error "Autoflake did not pass"
        exit 2
    fi
}

# --- Main Script ---

# Fortran Validation
check_fprettify
check_fortran_linter

# Python Validation
check_autopep8
check_docformatter
check_autoflake
