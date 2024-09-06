#!/usr/bin/env bash

set -e

# Source message handler functions
source ci/scripts/handle_messages.sh

# --- Fortran Section ---

# Function to get all Fortran files
get_fortran_files() {
    find . -type f \( -name "*.f90" -o -name "*.f" \) -not -path "./.sandbox/*"
}

# Function to format Fortran files with fprettify
format_fortran_files() {
    info "Formatting Fortran files with fprettify"
    fortran_files=$(get_fortran_files)
    if [ -z "$fortran_files" ]; then
        info "No Fortran files to format"
        return
    fi

    # Apply fprettify to all Fortran files
    for file in $fortran_files; do
        info "Formatting $file"
        fprettify --silent "$file"
    done
    info "Fortran formatting completed"
}

# --- Python Section ---

# Function to get all Python files
get_python_files() {
    find . -type f -name "*.py" -not -path "./.sandbox/*"
}

# Function to format Python files
format_python_files() {
    info "Formatting Python files"

    readarray -t python_files < <(get_python_files)
    if [ ${#python_files[@]} -eq 0 ]; then
        info "No Python files to format"
        return
    fi

    # Sorting imports with isort
    info "Sorting imports with isort"
    isort "${python_files[@]}"

    # Formatting Python files with autopep8
    info "Performing autopep8 on files"
    autopep8 --in-place "${python_files[@]}"

    # Performing selective aggressive autopep8
    info "Performing selective aggressive autopep8"
    autopep8 --global-config ci/.aggressive.pep8 --aggressive --in-place "${python_files[@]}"

    # Formatting docstrings with docformatter
    info "Performing docformatter"
    docformatter --wrap-summaries 88 --wrap-descriptions 88 --in-place "${python_files[@]}"

    # Removing unused imports and variables with autoflake
    info "Removing unused imports and variables with autoflake"
    autoflake --in-place --ignore-init-module-imports --remove-all-unused-imports --remove-unused-variables "${python_files[@]}"

    info "Python formatting completed"
}

# --- Main Script ---

# Format Fortran files
format_fortran_files

# Format Python files
format_python_files

info "All formatting completed"
